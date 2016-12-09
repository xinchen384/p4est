/*
 * p4est_setop2.c
 *
 *  Created on: Dec 6, 2016
 *      Author: foxjas09
 */

#include <p4est_vtk.h>
#include "rect_horizontal1.h"
#include "rect_vertical1.h"
#include <stdio.h>
#include <p4est_bits.h>

#define P4EST_STEP1_PATTERN_LEVEL 5
#define P4EST_STEP1_PATTERN_LENGTH (1 << P4EST_STEP1_PATTERN_LEVEL)
static const int plv = P4EST_STEP1_PATTERN_LEVEL;
static const int ple = P4EST_STEP1_PATTERN_LENGTH;
static const int WHITE = 0;
static const int BLACK = 1;

static int refine_fn_h(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  int                 tilelen;
  int                 offsi, offsj;
  int                 i, j;
  const char         *d;
  unsigned char       p[3];

  P4EST_ASSERT (which_tree == 0);
  P4EST_ASSERT (quadrant->level <= plv);

  tilelen = 1 << (plv - quadrant->level);       /* Pixel size of quadrant */
  offsi = quadrant->x / P4EST_QUADRANT_LEN (plv);       /* Pixel x offset */
  offsj = quadrant->y / P4EST_QUADRANT_LEN (plv);       /* Pixel y offset */
  P4EST_ASSERT (offsi >= 0 && offsj >= 0);
  int blackCount = 0;
  for (j = 0; j < tilelen; ++j) {
    P4EST_ASSERT (offsj + j < ple);
    for (i = 0; i < tilelen; ++i) {
      P4EST_ASSERT (offsi + i < ple);
      d = recth_header_data + 4 * (ple * (ple - 1 - (offsj + j)) + (offsi + i)); // TODO: print out?
      RECTH_HEADER_PIXEL (d, p);
      P4EST_ASSERT (p[0] == p[1] && p[1] == p[2]);      /* Grayscale image */
      if (p[0] < 128) {
        blackCount += 1;
      }
    }
  }

  if (!blackCount) {
	  quadrant->p.user_int = WHITE;
	  return 0;
  } else if (blackCount == 1) {
	  quadrant->p.user_int = BLACK;
	  return 0;
  } else {
	  return 1; // "gray" case; further refining needed
  }
}

static int refine_fn_v(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  int                 tilelen;
  int                 offsi, offsj;
  int                 i, j;
  const char         *d;
  unsigned char       p[3];

  P4EST_ASSERT (which_tree == 0);
  P4EST_ASSERT (quadrant->level <= plv);

  tilelen = 1 << (plv - quadrant->level);       /* Pixel size of quadrant */
  offsi = quadrant->x / P4EST_QUADRANT_LEN (plv);       /* Pixel x offset */
  offsj = quadrant->y / P4EST_QUADRANT_LEN (plv);       /* Pixel y offset */
  P4EST_ASSERT (offsi >= 0 && offsj >= 0);
  int blackCount = 0;
  for (j = 0; j < tilelen; ++j) {
    P4EST_ASSERT (offsj + j < ple);
    for (i = 0; i < tilelen; ++i) {
      P4EST_ASSERT (offsi + i < ple);
      d = rectv_header_data + 4 * (ple * (ple - 1 - (offsj + j)) + (offsi + i)); // TODO: print out?
      RECTV_HEADER_PIXEL (d, p);
      P4EST_ASSERT (p[0] == p[1] && p[1] == p[2]);      /* Grayscale image */
      if (p[0] < 128) {
        blackCount += 1;
      }
    }
  }

  if (!blackCount) {
	  quadrant->p.user_int = WHITE;
	  return 0;
  } else if (blackCount == 1) {
	  quadrant->p.user_int = BLACK;
	  return 0;
  } else {
	  return 1; // "gray" case; further refining needed
  }
}



int main (int argc, char **argv)
{
	int                 mpiret;
	int                 recursive, partforcoarsen, balance;
	sc_MPI_Comm         mpicomm;
	p4est_t            *p4est_h, *p4est_v, *p4est_u, *p4est_i;
	p4est_connectivity_t *conn_h, *conn_v, *conn_u, *conn_i;

	/* Initialize MPI; see sc_mpi.h. */
	mpiret = sc_MPI_Init (&argc, &argv);
	SC_CHECK_MPI (mpiret);
	mpicomm = sc_MPI_COMM_WORLD;

	/* Create single-quadtree forests */
	conn_h = p4est_connectivity_new_unitsquare();
	conn_v = p4est_connectivity_new_unitsquare();
	conn_u = p4est_connectivity_new_unitsquare();
	conn_i = p4est_connectivity_new_unitsquare();
	p4est_h = p4est_new (mpicomm, conn_h, 0, NULL, NULL);
	p4est_v = p4est_new (mpicomm, conn_v, 0, NULL, NULL);
	p4est_u = p4est_new (mpicomm, conn_u, 0, NULL, NULL);
	p4est_i = p4est_new (mpicomm, conn_i, 0, NULL, NULL);

	/* Refine the forest recursively in parallel. */
	recursive = 1;
	p4est_refine (p4est_h, recursive, refine_fn_h, NULL);
	p4est_refine (p4est_v, recursive, refine_fn_v, NULL);

	/* Union operation */
	p4est_union(p4est_h, p4est_v, p4est_u);
	// print union results
	p4est_tree_t *tree = p4est_tree_array_index (p4est_u->trees, 0);
	sc_array_t *quadrants = &(tree->quadrants);
	int n_quads = quadrants->elem_count;
	for (int si = 0; si < n_quads; si++) {
		p4est_quadrant_t *q = p4est_quadrant_array_index(quadrants, si);
		int value = q->p.user_int;
		if (value > 0) {
			printf("(%d, %d), lvl %d, ", q->x / P4EST_QUADRANT_LEN(plv), q->y / P4EST_QUADRANT_LEN(5), q->level);
			printf("value %d\n", value);
		}
	}

	/* Intersection operation */
	p4est_intersection(p4est_h, p4est_v, p4est_i);
	tree = p4est_tree_array_index (p4est_i->trees, 0);
	quadrants = &(tree->quadrants);
	n_quads = quadrants->elem_count;
	// print intersection results
	for (int si = 0; si < n_quads; si++) {
		p4est_quadrant_t *q = p4est_quadrant_array_index(quadrants, si);
		int value = q->p.user_int;
		if (value > 0) {
			printf("(%d, %d), lvl %d, ", q->x / P4EST_QUADRANT_LEN(plv), q->y / P4EST_QUADRANT_LEN(5), q->level);
			printf("value %d\n", value);
		}
	}

	/* Write the forest to disk for visualization */
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

	sc_finalize ();
	mpiret = sc_MPI_Finalize ();
	SC_CHECK_MPI (mpiret);
	return 0;
}

