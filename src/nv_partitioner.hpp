#pragma once

#include <vector>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <random>
#include <functional>

#include "partitioner.hpp"
#include "graph.hpp"
#include "util.hpp"
#include "min_heap.hpp"

/* Neighbor-expansion based Vertex partitioner */
template <typename GraphViewT>
class NvPartitioner : public Partitioner
{
  private:
    const vid_t num_vertices;
    const size_t num_edges;
    vid_t num_assigned;

    using AdjColorList = CSRGraph<uint24_t>::AdjWeightList;
    CSRGraph<uint24_t> graph;
    // the color of each edge in the adj list, determined by the src vertex
    std::vector<AdjColorList> adj_colors;

    const int p;
    std::vector<std::vector<vid_t>> partitions;
    std::vector<int> assignments;

    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<vid_t> dis;

    using assignment_func_t = std::function<int(const AdjColorList &)>;
    assignment_func_t assignment_f;

    int random_assignment(const AdjColorList &color_list) {
        return color_list[dis(gen) % color_list.size()].v;
    }

    int min_first_assignment(const AdjColorList &color_list) {
        size_t part_id = static_cast<size_t>(-1);
        size_t min_size = num_vertices;
        for (const auto color : color_list) {
            if (this->partitions[color.v].size() < min_size) {
                part_id = color.v;
                min_size = partitions[color.v].size();
            }
        }
        CHECK_NE(part_id, static_cast<size_t>(-1));
        return part_id;
    }

    int balanced_assignment(const AdjColorList &color_list) {
        // TODO: a combination of random and balanced assignments
        return 0;
    }

    void select_assignment_func(NE_POLICY nep) {
        switch(nep) {
            case NE_RANDOM:
                this->assignment_f = [this](const AdjColorList &list) {
                    return random_assignment(list);
                };
                LOG(INFO) << "neighbor expansion policy: random";
                break;
            case NE_MIN_FIRST:
                this->assignment_f = [this](const AdjColorList &list) {
                    return min_first_assignment(list);
                };
                LOG(INFO) << "neighbor expansion policy: min cluster first";
                break;
            case NE_BALANCED:
                this->assignment_f = [this](const AdjColorList &list) {
                    return min_first_assignment(list);
                };
                LOG(INFO) << "neighbor expansion policy: balanced";
                break;
            default:
                LOG(ERROR) << "invalid NE_POLICY value: " << nep;
        }
    }

    int new_assignment(const AdjColorList &color_list) {
        return this->assignment_f(color_list);
    }

  public:
    // init partitioner with no assigned vertices
    // the initial vertex in each partition is randomly picked
    NvPartitioner(vid_t nv, size_t ne, const GraphViewT &edges, int p, NE_POLICY nep)
    : num_vertices(nv), num_edges(ne), num_assigned(0), graph(nv, ne, edges)
    , p(p), partitions(p), assignments(nv, -1), rd(), gen(rd())
    {
        LOG(INFO) << "Initializing the vertex partitioner: p = " << p;
        CHECK_LT(p, 1<<24) << "Can't handle partition number more than " << (1<<24);

        for (vid_t vid = 0; vid < nv; ++vid) {
            adj_colors.push_back(graph.NWeights(vid));
            adj_colors[vid].clear(); // set length to 0
        }

        dis.param(std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));
        // sample p vertices as the starting point
        std::unordered_set<vid_t> sampled;
        while (sampled.size() < p)
            sampled.insert(dis(gen));
        LOG(INFO) << "Sample " << sampled.size() << " vertices as partitioning seeds";

        size_t part_id = 0;
        for (const auto vid : sampled) {
            assignments[vid] = part_id;
            partitions[part_id++].push_back(vid);
        }
        num_assigned += p;

        select_assignment_func(nep);
    }

    // init partitioner with some "seed" vertices already assigned to partitions
    NvPartitioner(vid_t nv, size_t ne, const GraphViewT &edges, int p,
        const std::vector<std::vector<vid_t>> &seeds, NE_POLICY nep) 
    : num_vertices(nv), num_edges(ne), num_assigned(0), graph(nv, ne, edges)
    , p(p), partitions(seeds), assignments(nv, -1), rd(), gen(rd())
    {
        LOG(INFO) << "Initializing the vertex partitioner: p = " << p;
        CHECK_LT(p, 1<<24) << "Can't handle partition number more than " << (1<<24);
        CHECK_EQ(p, seeds.size()) << "Seeds should be partitioned into " << p << " parts";

        for (vid_t vid = 0; vid < nv; ++vid) {
            adj_colors.push_back(graph.NWeights(vid));
            adj_colors[vid].clear(); // set length to 0
        }

        dis.param(std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));

        size_t part_id = 0;
        vid_t num_seeds = 0;
        for (const auto &part : partitions) {
            for (const vid_t vid : part)
                assignments[vid] = part_id;
            part_id++;
            num_seeds += part.size();
        }
        num_assigned += num_seeds;

        select_assignment_func(nep);
    }

    void split();
    void split_queue();

    // get partition assignments (0, ..., p-1) of each vertex from partitions vector
    std::vector<int> get_assignments() {
        return assignments;
    }
    std::vector<std::vector<vid_t>> get_vertex_sets() {
        return partitions;
    }
};

template <typename GraphViewT>
void NvPartitioner<GraphViewT>::split()
{
    Timer compute_timer;
    compute_timer.start();

    std::vector<vid_t> next_frontiers;
    std::vector<std::vector<vid_t>> frontiers = this->partitions;
    vid_t num_frontiers = 1;
    int rounds = 0;

    // keep expansion until no new nodes could be added
    while(num_frontiers != 0) {
        num_frontiers = 0;
        int part_id = 0;
        for (const auto &part : frontiers) {
            // mark each outgoing edge of the partition
            for (const vid_t src : part) {
                for (const vid_t dst : graph.N(src)) {
                    if (this->assignments[dst] == -1) {
                        this->adj_colors[dst].push_back(part_id);
                        next_frontiers.push_back(dst);
                    }
                }
            }
            ++part_id;
        }
        // generate the frontiers for the next round
        for (auto &part : frontiers) part.clear();
        for (const vid_t v : next_frontiers) {
            // there might be duplicate nodes in next_frontiers
            if (this->assignments[v] == -1) {
                int assign = new_assignment(this->adj_colors[v]);
                this->assignments[v] = assign;
                this->partitions[assign].push_back(v);
                frontiers[assign].push_back(v);
                num_frontiers++;
            }
        }
        num_assigned += num_frontiers;
        next_frontiers.clear();
        LOG(INFO) << "Neighbor expansion Round " << rounds++ << ": "
            << num_frontiers << " vertices got assigned";
    }

    compute_timer.stop();

    if (num_assigned < num_vertices)
        LOG(INFO) << num_assigned << "/" << num_vertices << " vertices got assigned";
    rep (i, p)
        DLOG(INFO) << "vertices in partition " << i << ": " << this->partitions[i].size();
    std::vector<size_t> counts;
    std::transform(this->partitions.begin(), this->partitions.end(),
        std::back_inserter(counts), [](std::vector<vid_t> &part) { return part.size(); });
    size_t max_vertices = *std::max_element(counts.begin(), counts.end());
    size_t min_vertices = *std::min_element(counts.begin(), counts.end());
    LOG(INFO) << "vertex balance (max/avg): " << (double)max_vertices / ((double)num_vertices/p);
    LOG(INFO) << "vertex balance (max/min): " << (double)max_vertices / (double)min_vertices;
    LOG(INFO) << "time used for partitioning: " << compute_timer.get_time();
}

template <typename GraphViewT>
void NvPartitioner<GraphViewT>::split_queue()
{
    Timer compute_timer;
    compute_timer.start();

    std::unordered_set<vid_t> frontiers;
    std::vector<bool> node_added(this->p, false);
    size_t node_added_count = 0;
    MinHeap<int, vid_t> frontier_queue;

    frontiers.reserve(std::min(num_assigned * (vid_t)num_edges / num_vertices,
                               num_vertices - num_assigned));
    int part_id = 0;
    for (const auto &part : this->partitions) {
        // mark each outgoing edge of the partition
        for (const vid_t src : part) {
            for (const vid_t dst : graph.N(src)) {
                if (this->assignments[dst] == -1) {
                    this->adj_colors[dst].push_back(part_id);
                    frontiers.insert(dst);
                }
            }
        }
        ++part_id;
    }
    frontier_queue.reserve(num_vertices);
    for (const auto node : frontiers)
        frontier_queue.insert(this->adj_colors[node].size(), node);
    frontiers.clear();


    while (frontier_queue.size() != 0) {
        // extract p least-connected nodes in the frontier,
        // assign them to smallest adjacent clusters, and update the frontier
        vid_t node;
        int conn;
        while (frontier_queue.get_min(conn, node)) {
            frontier_queue.remove(node);
            int assign = new_assignment(this->adj_colors[node]);
            if (!node_added[assign]) {
                this->assignments[node] = assign;
                this->partitions[assign].push_back(node);
                node_added[assign] = true;
                ++node_added_count;
                // update neighbors in frontier_queue
                for (const auto v : graph.N(node)) {
                    if (this->assignments[v] == -1) {
                        this->adj_colors[v].push_back(assign);
                        if (frontier_queue.contains(v)) {
                            frontier_queue.increase_key(v);
                        } else {
                            frontiers.insert(v);
                        }
                    }
                }
            } else {
                frontiers.insert(node);
            }
            if (node_added_count == this->p)
                break;
        }
        this->num_assigned += node_added_count;
        // more vertices to add into queue
        for (const auto node : frontiers)
            frontier_queue.insert(this->adj_colors[node].size(), node);
        frontiers.clear();

        LOG(INFO) << node_added_count << " vertices assigned";
        LOG(INFO) << "frontier queue size: " << frontier_queue.size();
        std::fill(node_added.begin(), node_added.end(), false);
        node_added_count = 0;
    }

    compute_timer.stop();

    if (num_assigned < num_vertices)
        LOG(INFO) << num_assigned << "/" << num_vertices << " vertices got assigned overall";
    rep (i, p)
        DLOG(INFO) << "vertices in partition " << i << ": " << this->partitions[i].size();
    std::vector<size_t> counts;
    std::transform(this->partitions.begin(), this->partitions.end(),
        std::back_inserter(counts), [](std::vector<vid_t> &part) { return part.size(); });
    size_t max_vertices = *std::max_element(counts.begin(), counts.end());
    size_t min_vertices = *std::min_element(counts.begin(), counts.end());
    LOG(INFO) << "vertex balance (max/avg): " << (double)max_vertices / ((double)num_vertices/p);
    LOG(INFO) << "vertex balance (max/min): " << (double)max_vertices / (double)min_vertices;
    LOG(INFO) << "time used for partitioning: " << compute_timer.get_time();
}