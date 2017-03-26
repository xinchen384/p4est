/*
 * Tests union and intersection operation over two 2D rectangle images.
 * Images are in "C header format" (see rect_horizontal1.h, rect_vertical1.h) for details.
 *
 *  Created on: Dec 6, 2016
 *  Author: James Fox
 */

#include <p4est_vtk.h>
#include "rect_horizontal1.h"
#include "rect_vertical1.h"
#include <stdio.h>
#include <p4est_bits.h>

#define P4EST_STEP1_PATTERN_LEVEL 5
#define P4EST_STEP1_PATTERN_LENGTH (1 << P4EST_STEP1_PATTERN_LEVEL) // image dimension: 2^5 = 32
static const int plv = P4EST_STEP1_PATTERN_LEVEL;
static const int ple = P4EST_STEP1_PATTERN_LENGTH;

/**
 * Flags used to "mark" quadrants during refinement process. Current union and intersection
 * implementation assume the following: quadrant value > 0 => quadrant of "marked" pixels,
 * quadrant value = 0 => quadrant of "unmarked" pixels (not pixels of interest).
 */
static const int UNMARKED = 0;
static const int MARKED = 1;

/**
 * User-implemented callback function, for refinement of horizontal image input (represented as
 * a single-tree forest).
 *
 * Here we assume grayscale image, where p[0] = p[1] = p[2].
 */
static int refine_fn_h(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  int                 tilelen;
  int                 offsi, offsj;
  int                 i, j;
  const char         *d;
  unsigned char       p[3];

  P4EST_ASSERT (quadrant->level <= plv); // shouldn't exceed max refinement level
  	  	  	  	  	  	  	  	  	  	 // which corresponds to individual pixel
  tilelen = 1 << (plv - quadrant->level);       /* Pixel size of quadrant */
  offsi = quadrant->x / P4EST_QUADRANT_LEN (plv);       /* Pixel x offset */
  offsj = quadrant->y / P4EST_QUADRANT_LEN (plv);       /* Pixel y offset */
  P4EST_ASSERT (offsi >= 0 && offsj >= 0);
  int markedCount = 0;
  for (j = 0; j < tilelen; ++j) {
    P4EST_ASSERT (offsj + j < ple);
    for (i = 0; i < tilelen; ++i) {
      P4EST_ASSERT (offsi + i < ple);
      d = recth_header_data + 4 * (ple * (ple - 1 - (offsj + j)) + (offsi + i));
      RECTH_HEADER_PIXEL (d, p);
      // pixel value 0 is most black, pixel value 255 is most black.
      // we assume black pixels are the pixels of interest.
      if (p[0] < 128) {
        markedCount += 1;
      }
    }
  }

  /**
   * After counting pixels in this quadrant/tile, need to make decision about whether to
   * refine further or not.
   */
  if (!markedCount) { // base case 1: no marked pixels in this tile
	  quadrant->p.user_int = UNMARKED;
	  return 0;
  } else if (tilelen*tilelen == 1 && markedCount == 1) { // base case 2: this single pixel is marked
	  quadrant->p.user_int = MARKED;
	  return 0;
  } else {
	  return 1; // not a base case; refine further
  }
}

/**
 * User-implemented callback function, for refinement of vertical image input
 *
 * Implementation and logic is identical to refine_fn_h, except for some names.
 * See above for detailed comments.
 */
static int refine_fn_v(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  int                 tilelen;
  int                 offsi, offsj;
  int                 i, j;
  const char         *d;
  unsigned char       p[3];

  P4EST_ASSERT (which_tree == 0);
  P4EST_ASSERT (quadrant->level <= plv);

  tilelen = 1 << (plv - quadrant->level);
  offsi = quadrant->x / P4EST_QUADRANT_LEN (plv);
  offsj = quadrant->y / P4EST_QUADRANT_LEN (plv);
  P4EST_ASSERT (offsi >= 0 && offsj >= 0);
  int markedCount = 0;
  for (j = 0; j < tilelen; ++j) {
    P4EST_ASSERT (offsj + j < ple);
    for (i = 0; i < tilelen; ++i) {
      P4EST_ASSERT (offsi + i < ple);
      d = rectv_header_data + 4 * (ple * (ple - 1 - (offsj + j)) + (offsi + i));
      RECTV_HEADER_PIXEL (d, p);
      if (p[0] < 128) {
        markedCount += 1;
      }
    }
  }

  if (!markedCount) {
	  quadrant->p.user_int = UNMARKED;
	  return 0;
  } else if (tilelen*tilelen == 1 && markedCount == 1) {
	  quadrant->p.user_int = MARKED;
	  return 0;
  } else {
	  return 1;
  }
}



int main (int argc, char **argv)
{
	int                 mpiret;
	int                 recursive;
	sc_MPI_Comm         mpicomm;
	p4est_t            *p4est_h, *p4est_v, *p4est_u, *p4est_i;
	p4est_connectivity_t *conn_h, *conn_v, *conn_u, *conn_i;
        int si ;

	/* Initialize MPI; "dummy" in single-process case */
	mpiret = sc_MPI_Init (&argc, &argv);
	SC_CHECK_MPI (mpiret);
	mpicomm = sc_MPI_COMM_WORLD;

	/* Create four forests: two input forests, one for union, one for intersection */
	conn_h = p4est_connectivity_new_unitsquare();
	conn_v = p4est_connectivity_new_unitsquare();
	conn_u = p4est_connectivity_new_unitsquare();
	conn_i = p4est_connectivity_new_unitsquare();
	p4est_h = p4est_new (mpicomm, conn_h, 0, NULL, NULL); // horizontal rectangle input
	p4est_v = p4est_new (mpicomm, conn_v, 0, NULL, NULL); // vertical rectangle input
	p4est_u = p4est_new (mpicomm, conn_u, 0, NULL, NULL); // rectangle union output
	p4est_i = p4est_new (mpicomm, conn_i, 0, NULL, NULL); // rectangle intersection output

	/* Refine two input forests (recursively) */
	recursive = 1;
	p4est_refine (p4est_h, recursive, refine_fn_h, NULL);
	p4est_refine (p4est_v, recursive, refine_fn_v, NULL);

	/* Union operation */
	p4est_union(p4est_h, p4est_v, p4est_u);

	// print union results
	p4est_tree_t *tree = p4est_tree_array_index (p4est_u->trees, 0); // we only need 1 tree for simple examples
	sc_array_t *quadrants = &(tree->quadrants);
	int n_quads = quadrants->elem_count;
	for ( si = 0; si < n_quads; si++) {
		p4est_quadrant_t *q = p4est_quadrant_array_index(quadrants, si);
		int value = q->p.user_int;
		if (value > 0) {
			printf("(%d, %d), lvl %d, ", q->x / P4EST_QUADRANT_LEN(plv), q->y / P4EST_QUADRANT_LEN(5), q->level);
			printf("value %d\n", value);
		}
	}

	/* Intersection operation */
	p4est_intersection(p4est_h, p4est_v, p4est_i);

	// print intersection results
	tree = p4est_tree_array_index (p4est_i->trees, 0);
	quadrants = &(tree->quadrants);
	n_quads = quadrants->elem_count;
	for (si = 0; si < n_quads; si++) {
		p4est_quadrant_t *q = p4est_quadrant_array_index(quadrants, si);
		int value = q->p.user_int;
		if (value > 0) {
			printf("(%d, %d), lvl %d, ", q->x / P4EST_QUADRANT_LEN(plv), q->y / P4EST_QUADRANT_LEN(5), q->level);
			printf("value %d\n", value);
		}
	}

	/* Write the forests to disk for visualization. Can open with ParaView */
	p4est_vtk_write_file (p4est_h, NULL, P4EST_STRING "_setop2_rectangle_h");
	p4est_vtk_write_file (p4est_v, NULL, P4EST_STRING "_setop2_rectangle_v");
	p4est_vtk_write_file (p4est_u, NULL, P4EST_STRING "_setop2_rectangle_union");
	p4est_vtk_write_file (p4est_i, NULL, P4EST_STRING "_setop2_rectangle_intersection");

	/* Destroy the p4est and the connectivity structure. */
	p4est_destroy (p4est_h);
	p4est_destroy (p4est_v);
	p4est_destroy (p4est_u);
	p4est_destroy (p4est_i);
	p4est_connectivity_destroy (conn_h);
	p4est_connectivity_destroy (conn_v);
	p4est_connectivity_destroy (conn_u);
	p4est_connectivity_destroy (conn_i);

	// MPI clean-up
	sc_finalize ();
	mpiret = sc_MPI_Finalize ();
	SC_CHECK_MPI (mpiret);
	return 0;
}

