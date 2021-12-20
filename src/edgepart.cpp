#include <string>

#include "edgepart.hpp"
#include "graph_view.hpp"
#include "ne_partitioner.hpp"

std::vector<std::vector<vid_t>> EdgePart_PartGraphCOO(vid_t nv, size_t ne, vid_t *row, vid_t *col, vid_t npart)
{
    NePartitioner<GraphViewRawCOO> partitioner(nv, ne, GraphViewRawCOO(row, col, ne), npart);
    partitioner.split();
    return partitioner.get_vertex_sets();
}