#include <unordered_map>
#include "nv_partitioner.hpp"

void NvPartitioner::partitioner_report()
{
    rep (i, p)
        DLOG(INFO) << "vertices in partition " << i << ": " << this->partitions[i].size();
    std::vector<size_t> counts;
    std::transform(this->partitions.begin(), this->partitions.end(),
        std::back_inserter(counts), [](std::vector<vid_t> &part) { return part.size(); });
    size_t max_vertices = *std::max_element(counts.begin(), counts.end());
    size_t min_vertices = *std::min_element(counts.begin(), counts.end());
    double avg_vertices = (double)std::accumulate(counts.begin(), counts.end(), 0) / p;
    LOG(INFO) << "vertex balance (max/avg): " << (double)max_vertices / avg_vertices;
    LOG(INFO) << "vertex balance (max/min): " << (double)max_vertices / (double)min_vertices;
    LOG(INFO) << "vertex balance (avg/min): " << avg_vertices / (double)min_vertices;
    size_t outliers = 0;
    for (auto c : counts) {
        outliers += static_cast<size_t>(c > 2 * avg_vertices);
    }
    LOG(INFO) << "number of outlier partitions: " << outliers;
}

void NvPartitioner::split()
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
    LOG(INFO) << "time used for partitioning: " << compute_timer.get_time();
    partitioner_report();
}

// a slow but guaranteed balanced partitioner
// it adds only 1 vertex to each cluster each round
void NvPartitioner::split_queue()
{
    Timer compute_timer;
    compute_timer.start();

    std::unordered_set<vid_t> frontiers;
    std::vector<bool> node_added(this->p, false);
    int node_added_count = 0;
    MinHeap<int, vid_t> frontier_queue;

    frontiers.reserve(std::min(num_assigned * (vid_t)num_edges / num_vertices,
                               num_vertices - num_assigned));
    int part_id = 0;
    for (const auto &part : this->partitions) {
        // mark each outgoing edge of the partition
        for (const vid_t src : part) {
            for (const vid_t dst : graph.N(src)) {
                if (this->assignments[dst] == -1) {
                    this->adj_clusters[dst].insert(part_id);
                    frontiers.insert(dst);
                }
            }
        }
        ++part_id;
    }
    frontier_queue.reserve(num_vertices);
    for (const auto node : frontiers)
        frontier_queue.insert(this->adj_clusters[node].size(), node);
    frontiers.clear();


    while (frontier_queue.size() != 0) {
        // extract p least-connected nodes in the frontier,
        // assign them to smallest adjacent clusters, and update the frontier
        vid_t node;
        int conn;
        while (frontier_queue.get_min(conn, node)) {
            frontier_queue.remove(node);
            // decide the cluster that could be assign node to
            int cluster_min = -1;
            size_t cluster_minsize = -1;
            for (auto cluster : this->adj_clusters[node]) {
                size_t cluster_size = this->partitions[cluster].size();
                if (!node_added[cluster] && cluster_size < cluster_minsize) {
                    cluster_min = cluster;
                    cluster_minsize = cluster_size;
                }
            }
            if (cluster_min != -1) {
                node_added[cluster_min] = true;
                this->assignments[node] = cluster_min;
                this->partitions[cluster_min].push_back(node);
                ++node_added_count;
                // update neighbors in frontier_queue
                for (const auto v : graph.N(node)) {
                    if (this->assignments[v] == -1) {
                        auto p = this->adj_clusters[v].insert(cluster_min);
                        if (p.second && frontier_queue.contains(v)) {
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
            frontier_queue.insert(this->adj_clusters[node].size(), node);
        frontiers.clear();

        LOG(INFO) << node_added_count << " vertices assigned";
        LOG(INFO) << "frontier queue size: " << frontier_queue.size();
        std::fill(node_added.begin(), node_added.end(), false);
        node_added_count = 0;
    }

    compute_timer.stop();
    if (num_assigned < num_vertices)
        LOG(INFO) << num_assigned << "/" << num_vertices << " vertices got assigned";
    LOG(INFO) << "time used for partitioning: " << compute_timer.get_time();
    partitioner_report();
}

void NvPartitioner::split_balanced() {
    Timer compute_timer;
    compute_timer.start();

    std::vector<std::vector<vid_t>> frontiers = this->partitions;
    std::unordered_map<vid_t, int> candidates(this->p);
    std::unordered_set<vid_t> new_candidates(this->p);
    std::vector<size_t> expansion_counts(this->p);
    vid_t num_frontiers = 1;
    int rounds = 0;
    size_t min_thresold = this->num_vertices / this->p * 0.01;
    // keep expansion until no new nodes could be added
    while(num_frontiers != 0) {
        num_frontiers = 0;
        int part_id = 0;
        for (const auto &part : frontiers) {
            // mark each outgoing edge of the partition
            for (const vid_t src : part) {
                for (const vid_t dst : graph.N(src)) {
                    if (this->assignments[dst] == -1) {
                        // update color list of unassigned nodes
                        this->adj_colors[dst].push_back(part_id);
                        new_candidates.insert(dst);
                        // add newly added nodes
                        if (candidates.find(dst) == candidates.end()) {
                            candidates.insert({dst, 0});
                        }
                    }
                }
            }
            ++part_id;
        }

        // decide assignments for newly added nodes
        for (const auto v : new_candidates) {
            int part = new_assignment(this->adj_colors[v]);
            candidates[v] = part;
        }
        // find min non-zero expansion elements in all clusters
        for (const auto &p : candidates)
            expansion_counts[p.second]++;
        size_t min_count = -1;
        for (auto count : expansion_counts)
            if (count != 0 && count < min_count)
                min_count = count;
        if (min_count == (size_t)-1) {
            // all elements in delta_assignments are zero
            break;
        }
        min_count = (min_count < min_thresold) ? min_thresold : min_count;

        // cluster expansion
        for (auto &part : frontiers) part.clear();
        for (auto it = candidates.cbegin(); it != candidates.cend();) {
            // randomly select vertices to add into the next frontier
            // each cluster gets min_count new nodes expectedly
            vid_t v = it->first;
            int part = it->second;
            if (dis(gen) % expansion_counts[part] < min_count) {
                frontiers[part].push_back(v);
                this->assignments[v] = part;
                this->partitions[part].push_back(v);
                it = candidates.erase(it);
                num_frontiers++;
            } else {
                ++it;
            }
        }

        num_assigned += num_frontiers;
        std::fill(expansion_counts.begin(), expansion_counts.end(), 0);
        new_candidates.clear();
        LOG(INFO) << "NEV Round " << rounds++ << ": " << num_frontiers
            << " vertices got assigned, " << num_assigned << " overall";
    }

    compute_timer.stop();

    if (num_assigned < num_vertices)
        LOG(INFO) << num_assigned << "/" << num_vertices << " vertices got assigned";
    LOG(INFO) << "time used for partitioning: " << compute_timer.get_time();
    partitioner_report();
}

void NvPartitioner::split_balanced_v2()
{
    Timer compute_timer;
    compute_timer.start();

    std::vector<std::vector<vid_t>> frontiers = this->partitions;
    std::unordered_set<vid_t> new_candidates(this->p);
    std::vector<float> adj_prob;
    vid_t num_frontiers = -1;
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
                        // update color list of unassigned nodes
                        // this->adj_colors[dst].push_back(part_id);
                        this->adj_clusters[dst].insert(part_id);
                        new_candidates.insert(dst);
                    }
                }
            }
            ++part_id;
        }
        for (auto &part : frontiers) part.clear();

        // initialize adj_prob matrix (dense format)
        // if we are to use CSR, we have to do CSR <-> CSC transformation repetitively
        // TODO: use CSR & CSC as index to the dense format matrix
        size_t row_size = new_candidates.size();
        adj_prob.resize(this->p * row_size);
        adj_prob.clear();
        for (auto iter = new_candidates.begin(); iter != new_candidates.end(); ++iter) {
            const auto row = std::distance(new_candidates.begin(), iter);
            for (const auto part : this->adj_clusters[*iter]) {
                adj_prob[this->p * row + part] = 1;
            }
        }
        float target_col_sum = (float)row_size / this->p;
        std::vector<float> col_sum(row_size);

        // 3 iterations for now
        for (int it = 0; it < 3; ++it) {
            col_sum.clear();
            // column normalize
            for (auto i = 0; i < row_size; ++i) {
                std::transform(col_sum.begin(), col_sum.end(), &adj_prob[this->p * i],
                    col_sum.begin(), std::plus<float>{});
            }
            for (auto i = 0; i < row_size; ++i) {
                std::transform(col_sum.begin(), col_sum.end(), &adj_prob[this->p * i],
                    &adj_prob[this->p * i], [target_col_sum](float x, float y)->float { return y/x*target_col_sum;});
            }
            // row normalize
            for (auto i = 0; i < row_size; ++i) {
                auto sum = std::accumulate(&adj_prob[this->p * i], &adj_prob[this->p * (i+1)], 0);
                std::transform(&adj_prob[this->p*i], &adj_prob[this->p*(i+1)], &adj_prob[this->p*i],
                    [sum](float prob)->float { return prob/sum; });
            }
        }

        // determine assignments based on adj_prob
        for (auto iter = new_candidates.begin(); iter != new_candidates.end(); ++iter) {
            const auto row = std::distance(new_candidates.begin(), iter);
            vid_t v = *iter;
            float rand = dis(gen) / (float)this->num_vertices;
            // TODO: int assign = ...
            int assign = new_assignment(this->adj_colors[v]);
            this->assignments[v] = assign;
            this->partitions[assign].push_back(v);
            frontiers[assign].push_back(v);
            num_frontiers++;
        }

        new_candidates.clear();
        num_assigned += num_frontiers;
        LOG(INFO) << "Neighbor expansion Round " << rounds++ << ": "
            << num_frontiers << " vertices got assigned";
    }

    compute_timer.stop();
    LOG(INFO) << num_assigned << "/" << num_vertices << " vertices got assigned";
    LOG(INFO) << "time used for partitioning: " << compute_timer.get_time();
    partitioner_report();
}
