#pragma once

#include <vector>
#include <memory>
#include <parallel/algorithm>

#include "graph_view.hpp"
#include "util.hpp"

// TODO: don't use 40bits - could be insufficient soon in the future
struct uint40_t {
        uint64_t v:40;
} __attribute__((packed));

class adjlist_t
{
  private:
    uint40_t *adj;
    vid_t len;

  public:
    adjlist_t() : adj(NULL), len(0) {}
    adjlist_t(uint40_t *adj, vid_t len = 0) : adj(adj), len(len) {}
    uint40_t *begin() { return adj; }
    uint40_t *end() { return adj + len; }
    void increment() { len++; }
    void push_back(size_t data) { adj[len++].v = data; }
    size_t size() const { return len; }
    uint40_t &operator[](size_t idx) { return adj[idx]; };
    const uint40_t &operator[](size_t idx) const { return adj[idx]; };
    uint40_t &back() { return adj[len - 1]; };
    const uint40_t &back() const { return adj[len - 1]; };
    void pop_back() { len--; }
    void clear() { len = 0; }
};

class graph_t
{
  private:
    vid_t num_vertices;
    size_t nedges;
    uint40_t *neighbors;
    std::vector<adjlist_t> vdata;

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

    adjlist_t &operator[](size_t idx) { return vdata[idx]; };
    const adjlist_t &operator[](size_t idx) const { return vdata[idx]; };
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

    vdata[0] = adjlist_t(neighbors);
    for (vid_t v = 1; v < num_vertices; v++) {
        count[v] += count[v-1];
        vdata[v] = adjlist_t(neighbors + count[v-1]);
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

    vdata[0] = adjlist_t(neighbors);
    for (vid_t v = 1; v < num_vertices; v++) {
        count[v] += count[v-1];
        vdata[v] = adjlist_t(neighbors + count[v-1]);
    }
    for (size_t i = 0; i < edges.size(); i++)
        vdata[edges[i].second].push_back(i);
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
