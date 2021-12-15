#include "graph.hpp"
#include "conversions.hpp"
#include "util.hpp"

/*
void graph_t::build(const std::vector<edge_t> &edges)
{
    if (edges.size() > nedges)
        neighbors = (uint40_t *)realloc(neighbors, sizeof(uint40_t) * edges.size());
    CHECK(neighbors) << "allocation failed";
    nedges = edges.size();

    std::vector<size_t> count(num_vertices, 0);
    for (size_t i = 0; i < nedges; i++)
        count[edges[i].first]++;

    vdata[0] = adjlist_t(neighbors);
    for (vid_t v = 1; v < num_vertices; v++) {
        count[v] += count[v-1];
        vdata[v] = adjlist_t(neighbors + count[v-1]);
    }
    // vdata[from] is an adjlist storing eid's of the array "edges"
    for (size_t i = 0; i < edges.size(); i++)
        vdata[edges[i].first].push_back(i);
}

void graph_t::build_reverse(const std::vector<edge_t> &edges)
{
    if (edges.size() > nedges)
        neighbors = (uint40_t *)realloc(neighbors, sizeof(uint40_t) * edges.size());
    CHECK(neighbors) << "allocation failed";
    nedges = edges.size();

    std::vector<size_t> count(num_vertices, 0);
    for (size_t i = 0; i < nedges; i++)
        count[edges[i].second]++;

    vdata[0] = adjlist_t(neighbors);
    for (vid_t v = 1; v < num_vertices; v++) {
        count[v] += count[v-1];
        vdata[v] = adjlist_t(neighbors + count[v-1]);
    }
    for (size_t i = 0; i < edges.size(); i++)
        vdata[edges[i].second].push_back(i);
}
*/

void FileGraphEdgeList::load()
{
    // converter outputs 2 files
    // - vid is consecutive, starting from 0
    // - edgelist in binary: [num_v : uint32_t, num_e : uint64_t, [(from, to : uint32_t), ...] ]
    // - degree in binary: [deg : uint32_t, ...]
    Timer convert_timer;
    convert_timer.start();
    Converter *converter = new Converter(basefilename);
    convert(basefilename, converter);
    delete converter;
    convert_timer.stop();
    LOG(INFO) << "convert time: " << convert_timer.get_time();

    Timer read_timer;
    read_timer.start();
    LOG(INFO) << "loading...";

    std::ifstream fin(binedgelist_name(basefilename),
                    std::ios::binary | std::ios::ate);
    auto filesize = fin.tellg();
    LOG(INFO) << "file size: " << filesize;
    fin.seekg(0, std::ios::beg);

    fin.read((char *)&num_vertices, sizeof(num_vertices));
    fin.read((char *)&num_edges, sizeof(num_edges));

    LOG(INFO) << "num_vertices: " << num_vertices
            << ", num_edges: " << num_edges;
    CHECK_EQ(sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t), filesize);

    edges.resize(num_edges);
    fin.read((char *)edges.data(), sizeof(edge_t) * num_edges);

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)degrees.data(), num_vertices * sizeof(vid_t));
    degree_file.close();
    read_timer.stop();
    LOG(INFO) << "load time: " << read_timer.get_time();
}

FileGraphEdgeList::GraphViewT FileGraphEdgeList::get_view()
{
    return GraphViewT(this->edges);
}

void FileGraphCOO::load()
{
    // converter outputs 2 files
    // - vid is consecutive, starting from 0
    // - edgelist in binary: [num_v : uint32_t, num_e : uint64_t, [(from, to : uint32_t), ...] ]
    // - degree in binary: [deg : uint32_t, ...]
    Timer convert_timer;
    convert_timer.start();
    Converter *converter = new Converter(basefilename);
    convert(basefilename, converter);
    delete converter;
    convert_timer.stop();
    LOG(INFO) << "convert time: " << convert_timer.get_time();

    Timer read_timer;
    read_timer.start();
    LOG(INFO) << "loading...";

    std::ifstream fin(binedgelist_name(basefilename),
                    std::ios::binary | std::ios::ate);
    auto filesize = fin.tellg();
    LOG(INFO) << "file size: " << filesize;
    fin.seekg(0, std::ios::beg);

    fin.read((char *)&num_vertices, sizeof(num_vertices));
    fin.read((char *)&num_edges, sizeof(num_edges));

    LOG(INFO) << "num_vertices: " << num_vertices
            << ", num_edges: " << num_edges;
    CHECK_EQ(sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t), filesize);

    std::vector<edge_t> edges(num_edges);
    fin.read((char *)edges.data(), sizeof(edge_t) * num_edges);

    row.reserve(num_edges);
    col.reserve(num_edges);
    for (const auto &edge: edges) {
        row.push_back(edge.first);
        col.push_back(edge.second);
    }

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)degrees.data(), num_vertices * sizeof(vid_t));
    degree_file.close();
    read_timer.stop();
    LOG(INFO) << "load time: " << read_timer.get_time();
}

FileGraphCOO::GraphViewT FileGraphCOO::get_view()
{
    return GraphViewT(this->row, this->col);
}
