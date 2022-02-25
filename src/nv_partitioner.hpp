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
    std::vector<std::unordered_set<int>> adj_clusters;

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
    template <typename GraphViewT>
    NvPartitioner(vid_t nv, size_t ne, const GraphViewT &edges, int p, NE_POLICY nep)
    : num_vertices(nv), num_edges(ne), num_assigned(0), graph(nv, ne, edges)
    , p(p), partitions(p), assignments(nv, -1), rd(), gen(rd())
    {
        LOG(INFO) << "Initializing the vertex partitioner: p = " << p;
        CHECK_LT(p, 1<<24) << "Can't handle partition number more than " << (1<<24);

        for (vid_t vid = 0; vid < nv; ++vid) {
            adj_colors.push_back(graph.NWeights(vid));
            adj_colors[vid].clear(); // set length to 0
            adj_clusters.push_back(std::unordered_set<int>(std::min(p, (int)graph.num_neighbors(vid))));
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
    template <typename GraphViewT>
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
            adj_clusters.push_back(std::unordered_set<int>(std::min(p, (int)graph.num_neighbors(vid))));
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

    void partitioner_report();
    void split();
    void split_queue();
    void split_balanced();
    void split_balanced_v2();

    // get partition assignments (0, ..., p-1) of each vertex from partitions vector
    std::vector<int> get_assignments() {
        return assignments;
    }
    std::vector<std::vector<vid_t>> get_vertex_sets() {
        return partitions;
    }
};
