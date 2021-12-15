#pragma once

#include <stddef.h>
#include <stdint.h>
#include <vector>

typedef int32_t vid_t;

/*
 * The input COO graph should be bidirected (or undirected but also unsymmetrized):
 * - vertex id starts consecutively from 0 to num_vertices - 1
 * - edges are specified by vid_t *row and vid_t *col
 * 
 * vid_t *row, vid_t *col is not managed by this function
 * initially, part_vertices and vertex_counts should be empty
 */
std::vector<std::vector<vid_t>> EdgePart_PartGraphCOO(vid_t nv, size_t ne, vid_t *row, vid_t *col, vid_t npart);
