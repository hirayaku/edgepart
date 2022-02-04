#pragma once

#include <stddef.h>
#include <stdint.h>
#include <vector>

typedef int64_t vid_t;

enum NE_POLICY {
    NE_RANDOM = 0,
    NE_MIN_FIRST = 1,
    NE_BALANCED = 2
};

/*
 * For the following APIs, the input COO graph should be bidirected
 * (or undirected with each edge repeated once)
 * - vertex id starts consecutively from 0 to num_vertices - 1
 * - edges are specified by vid_t *row and vid_t *col
 * - memory referred by vid_t *row, vid_t *col is not managed by this function
 */

std::vector<std::vector<vid_t>> EdgePart_PartGraphCOO(
    vid_t nv, size_t ne, vid_t *row, vid_t *col, vid_t npart);

std::vector<std::vector<vid_t>> VertexPart_PartGraphCOOFromSeeds(
    vid_t nv, size_t ne, vid_t *row, vid_t *col, vid_t npart,
    const std::vector<vid_t> &, const std::vector<vid_t> &, NE_POLICY);