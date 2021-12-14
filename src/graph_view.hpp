#include <vector>
#include <memory>
#include "dense_bitset.hpp"

typedef std::vector<edge_t> GraphEdgeList;

class GraphViewEdgeList
{
  private:
    const std::vector<edge_t> &edges;
    std::unique_ptr<dense_bitset> valid_bits;

  public:
    struct edge_ref {
      const vid_t &first;
      const vid_t &second;
      edge_ref(GraphViewEdgeList &edgelist, size_t idx)
      : first(edgelist.edges[idx].first)
      , second(edgelist.edges[idx].second)
      , valid_bits(*edgelist.valid_bits)
      , idx(idx) {}
      bool valid() const { return valid_bits.get(idx); }
      void remove() { valid_bits.clear_bit_unsync(idx); }

      private:
        dense_bitset &valid_bits;
        size_t idx;
    };

    struct edge_cref {
      const vid_t &first;
      const vid_t &second;
      edge_cref(const GraphViewEdgeList &edgelist, size_t idx)
      : first(edgelist.edges[idx].first)
      , second(edgelist.edges[idx].second)
      {}
    };

    GraphViewEdgeList() = delete;
    GraphViewEdgeList(GraphViewEdgeList &) = delete;
    GraphViewEdgeList(GraphViewEdgeList &&that)
    : edges(that.edges), valid_bits(std::move(that.valid_bits))
    {}
    GraphViewEdgeList(const std::vector<edge_t> &edges)
    : edges(edges), valid_bits(new dense_bitset(edges.size()))
    { valid_bits->fill(); }

    size_t size() const { return edges.size(); }
    edge_ref operator[](size_t idx) { return edge_ref(*this, idx); }
    edge_cref operator[](size_t idx) const { return edge_cref(*this, idx); }
};

class GraphViewCOO
{
  private:
    const std::vector<vid_t> &src;
    const std::vector<vid_t> &dst;
    std::unique_ptr<dense_bitset> valid_bits;
  
  public:
    struct edge_ref {
      const vid_t &first;
      const vid_t &second;
      edge_ref(GraphViewCOO &edgelist, size_t idx)
      : first(edgelist.src[idx])
      , second(edgelist.dst[idx])
      , valid_bits(*edgelist.valid_bits)
      , idx(idx) {}
      bool valid() const { return valid_bits.get(idx); }
      void remove() { valid_bits.clear_bit_unsync(idx); }

      private:
        dense_bitset &valid_bits;
        size_t idx;
    };

    struct edge_cref {
      const vid_t &first;
      const vid_t &second;
      edge_cref(const GraphViewCOO &edgelist, size_t idx)
      : first(edgelist.src[idx]), second(edgelist.dst[idx])
      {}
    };

    GraphViewCOO() = delete;
    GraphViewCOO(GraphViewCOO &) = delete;
    GraphViewCOO(GraphViewCOO &&that)
    : src(that.src), dst(that.dst), valid_bits(std::move(that.valid_bits))
    {}
    GraphViewCOO(const std::vector<vid_t> &src, const std::vector<vid_t> &dst)
    : src(src), dst(dst), valid_bits(new dense_bitset(src.size()))
    {
      CHECK_EQ(src.size(), dst.size());
      valid_bits->fill();
    }
    size_t size() const { return src.size(); }
    edge_ref operator[](size_t idx) { return edge_ref(*this, idx); }
    edge_cref operator[](size_t idx) const { return edge_cref(*this, idx); }
};

// using COO with raw pointers
class GraphViewRawCOO
{
  private:
    const vid_t *&src;
    const vid_t *&dst;
    const size_t array_len;
    std::unique_ptr<dense_bitset> valid_bits;
  
  public:
    struct edge_ref {
      const vid_t &first;
      const vid_t &second;
      edge_ref(GraphViewRawCOO &edgelist, size_t idx)
      : first(edgelist.src[idx])
      , second(edgelist.dst[idx])
      , valid_bits(*edgelist.valid_bits)
      , idx(idx) {}
      bool valid() const { return valid_bits.get(idx); }
      void remove() { valid_bits.clear_bit_unsync(idx); }

      private:
        dense_bitset &valid_bits;
        size_t idx;
    };

    struct edge_cref {
      const vid_t &first;
      const vid_t &second;
      edge_cref(const GraphViewRawCOO &edgelist, size_t idx)
      : first(edgelist.src[idx]), second(edgelist.dst[idx])
      {}
    };

    GraphViewRawCOO() = delete;
    GraphViewRawCOO(GraphViewRawCOO &) = delete;
    GraphViewRawCOO(GraphViewRawCOO &&that)
    : src(that.src), dst(that.dst), array_len(that.array_len)
    , valid_bits(std::move(that.valid_bits))
    {}
    GraphViewRawCOO(const vid_t *src, const vid_t *dst, const size_t array_len)
    : src(src), dst(dst), array_len(array_len)
    , valid_bits(new dense_bitset(array_len))
    { valid_bits->fill(); }
    size_t size() const { return array_len; }
    edge_ref operator[](size_t idx) { return edge_ref(*this, idx); }
    edge_cref operator[](size_t idx) const { return edge_cref(*this, idx); }
};
