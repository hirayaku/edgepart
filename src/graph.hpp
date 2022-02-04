#pragma once

#include <vector>
#include <memory>
#include <algorithm>

#include "graph_view.hpp"
#include "util.hpp"

// TODO: don't use 40bits for EID - could be insufficient soon in the future
struct uint40_t {
        uint64_t v:40;
} __attribute__((packed));

struct uint24_t {
        uint32_t v:24;
} __attribute__((packed));

template <typename T>
class adjlist_t
{
  private:
    T *adj;
    vid_t len;

  public:
    adjlist_t() : adj(NULL), len(0) {}
    adjlist_t(T *adj, vid_t len = 0) : adj(adj), len(len) {}
    T *begin() { return adj; }
    T *end() { return adj + len; }
    T *begin() const { return adj; }
    T *end() const { return adj + len; }
    void increment() { len++; }
    void push_back(size_t data) { adj[len++].v = data; }
    size_t size() const { return len; }
    T &operator[](size_t idx) { return adj[idx]; };
    const T &operator[](size_t idx) const { return adj[idx]; };
    T &back() { return adj[len - 1]; };
    const T &back() const { return adj[len - 1]; };
    void pop_back() { len--; }
    void clear() { len = 0; }
};

class graph_t
{
  private:
    using adjlist = adjlist_t<uint40_t>;
    vid_t num_vertices;
    size_t nedges;
    uint40_t *neighbors;
    std::vector<adjlist> vdata;

  public:
    graph_t() : num_vertices(0), nedges(0), neighbors(NULL) {}

    ~graph_t()
    {
        if (neighbors)
            free(neighbors);
    }

    void resize(vid_t _num_vertices)
    {
        num_vertices = _num_vertices;
        vdata.resize(num_vertices);
    }

    size_t num_edges() const { return nedges; }

    template <typename GraphViewT>
    void build(const GraphViewT &edges);

    template <typename GraphViewT>
    void build_reverse(const GraphViewT &edges);

    adjlist &operator[](size_t idx) { return vdata[idx]; };
    const adjlist &operator[](size_t idx) const { return vdata[idx]; };
};

template <typename GraphViewT>
void graph_t::build(const GraphViewT &edges)
{
    if (edges.size() > nedges)
        neighbors = (uint40_t *)realloc(neighbors, sizeof(uint40_t) * edges.size());
    CHECK(neighbors) << "allocation failed";
    nedges = edges.size();

    std::vector<size_t> count(num_vertices, 0);
    for (size_t i = 0; i < nedges; i++)
        count[edges[i].first]++;

    vdata[0] = adjlist(neighbors);
    for (vid_t v = 1; v < num_vertices; v++) {
        count[v] += count[v-1];
        vdata[v] = adjlist(neighbors + count[v-1]);
    }
    // vdata[from] is an adjlist storing eid's of the array "edges"
    for (size_t i = 0; i < edges.size(); i++)
        vdata[edges[i].first].push_back(i);
}

template <typename GraphViewT>
void graph_t::build_reverse(const GraphViewT &edges)
{
    if (edges.size() > nedges)
        neighbors = (uint40_t *)realloc(neighbors, sizeof(uint40_t) * edges.size());
    CHECK(neighbors) << "allocation failed";
    nedges = edges.size();

    std::vector<size_t> count(num_vertices, 0);
    for (size_t i = 0; i < nedges; i++) {
        CHECK_LT(edges[i].second, num_vertices);
        count[edges[i].second]++;
    }

    vdata[0] = adjlist(neighbors);
    for (vid_t v = 1; v < num_vertices; v++) {
        count[v] += count[v-1];
        vdata[v] = adjlist(neighbors + count[v-1]);
    }
    for (size_t i = 0; i < edges.size(); i++)
        vdata[edges[i].second].push_back(i);
};

template <typename EDATA_T>
class CSRGraph {
  private:
    vid_t nvertices;
    size_t nedges;
    std::vector<size_t> edge_offsets;
    std::vector<vid_t> edge_vids;
    std::vector<EDATA_T> edge_data;

  public:
    using AdjList = adjlist_t<vid_t>;
    using AdjWeightList = adjlist_t<EDATA_T>;

    CSRGraph(): nvertices(0), nedges(0), edge_offsets(0), edge_vids(0), edge_data(0) {}

    template <typename GraphViewT>
    CSRGraph(vid_t nv, size_t ne, const GraphViewT &edges, const EDATA_T *data = nullptr)
    : nvertices(nv), nedges(ne), edge_offsets(nv+1), edge_vids(ne), edge_data(ne)
    {
      std::vector<size_t> count(nv, 0);
      for (size_t i = 0; i < nedges; i++)
          count[edges[i].first]++;

      edge_offsets[0] = 0;
      for (vid_t i = 1; i <= nv; ++i)
        edge_offsets[i] = edge_offsets[i-1] + count[i-1];

      count.clear();
      count = edge_offsets;
      for (size_t i = 0; i < nedges; ++i) {
        auto edge = edges[i];
        edge_vids[count[edge.first]++] = edge.second;
      }

      if (data) {
        edge_data.assign(data, data+ne);
      }
    }
    ~CSRGraph() {}

    vid_t num_vertices() const { return nvertices; }
    size_t num_edges() const { return nedges; }
    vid_t num_neighbors(vid_t src) const {
      return static_cast<vid_t>(edge_offsets[src+1] - edge_offsets[src]);
    }
    AdjList N(vid_t src) {
      return AdjList(edge_vids.data() + edge_offsets[src], num_neighbors(src));
    }
    AdjWeightList NWeights(vid_t src) {
      return AdjWeightList(edge_data.data() + edge_offsets[src], num_neighbors(src));
    }
};

// abstract class for graph objects backed up by a file
struct FileGraph
{ 
  const std::string basefilename;
  FileGraph(std::string basefilename): basefilename(basefilename) {}
  virtual ~FileGraph() {}
  virtual void load() = 0; // must specify how to load data from the file
};

struct FileGraphEdgeList : public FileGraph
{
  using GraphViewT = GraphViewEdgeList;
  vid_t num_vertices;
  size_t num_edges;
  std::vector<edge_t> edges;
  std::vector<vid_t> degrees;

  FileGraphEdgeList(std::string basefilename): FileGraph(basefilename) {}
  virtual void load();
  GraphViewT get_view();
};

struct FileGraphCOO : public FileGraph
{
  using GraphViewT = GraphViewCOO;
  vid_t num_vertices;
  size_t num_edges;
  std::vector<vid_t> row;
  std::vector<vid_t> col;
  std::vector<vid_t> degrees;

  FileGraphCOO(std::string basefilename): FileGraph(basefilename) {}
  virtual void load();
  GraphViewT get_view();
};

// struct FileGraphMemmap: public FileGraph
// ...

template <typename FILE_GRAPH>
FileGraph *makeFileGraph(std::string basefilename)
{
  FileGraph *p = new FILE_GRAPH(basefilename);
  p->load();
  return p;
}
