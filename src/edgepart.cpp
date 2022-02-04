#include <vector>

#include "edgepart.hpp"
#include "graph_view.hpp"
#include "ne_partitioner.hpp"
#include "nv_partitioner.hpp"

std::vector<std::vector<vid_t>> EdgePart_PartGraphCOO(vid_t nv, size_t ne, vid_t *row, vid_t *col, vid_t npart)
{
    NePartitioner<GraphViewRawCOO> partitioner(nv, ne, GraphViewRawCOO(row, col, ne), npart);
    partitioner.split();
    return partitioner.get_vertex_sets();
}

std::vector<std::vector<vid_t>> VertexPart_PartGraphCOOFromSeeds(
    vid_t nv, size_t ne, vid_t *row, vid_t *col, vid_t npart,
    const std::vector<vid_t> &seed_vertices, const std::vector<vid_t> &seed_assignments,
    NE_POLICY nes)
{
    std::vector<std::vector<vid_t>> seeds(npart);
    for (size_t i = 0; i < seed_vertices.size(); ++i) {
        seeds[seed_assignments[i]].push_back(seed_vertices[i]);
    }
    NvPartitioner<GraphViewRawCOO> partitioner(nv, ne, GraphViewRawCOO(row, col, ne), npart, seeds, nes);

    if (nes == NE_BALANCED)
        partitioner.split_queue();
    else
        partitioner.split();
    return partitioner.get_vertex_sets();
}
