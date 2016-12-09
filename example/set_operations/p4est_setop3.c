/*
 * p4est_setop3.c
 *
 *  Created on: Dec 8, 2016
 *      Author: foxjas09
 */

#include <p4est_vtk.h>
#include "cat1.h"
#include "cat2.h"
#include <stdio.h>
#include <p4est_bits.h>

#define P4EST_STEP1_PATTERN_LEVEL 7
#define P4EST_STEP1_PATTERN_LENGTH (1 << P4EST_STEP1_PATTERN_LEVEL)
static const int plv = P4EST_STEP1_PATTERN_LEVEL;
static const int ple = P4EST_STEP1_PATTERN_LENGTH;
static const int WHITE = 0;
static const int BLACK = 1;

static int refine_fn1(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
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
      d = cat1_header_data + 4 * (ple * (ple - 1 - (offsj + j)) + (offsi + i)); // TODO: print out?
      CAT1_HEADER_PIXEL (d, p);
      P4EST_ASSERT (p[0] == p[1] && p[1] == p[2]);
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

static int refine_fn2(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
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
      d = cat2_header_data + 4 * (ple * (ple - 1 - (offsj + j)) + (offsi + i)); // TODO: print out?
      CAT2_HEADER_PIXEL (d, p);
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
	p4est_t            *p4est1, *p4est2, *p4est_u, *p4est_i;
	p4est_connectivity_t *conn1, *conn2, *conn_u, *conn_i;

	/* Initialize MPI; see sc_mpi.h. */
	mpiret = sc_MPI_Init (&argc, &argv);
	SC_CHECK_MPI (mpiret);
	mpicomm = sc_MPI_COMM_WORLD;

	/* Create single-quadtree forests */
	conn1 = p4est_connectivity_new_unitsquare();
	conn2 = p4est_connectivity_new_unitsquare();
	conn_u = p4est_connectivity_new_unitsquare();
	conn_i = p4est_connectivity_new_unitsquare();
	p4est1 = p4est_new (mpicomm, conn1, 0, NULL, NULL);
	p4est2 = p4est_new (mpicomm, conn2, 0, NULL, NULL);
	p4est_u = p4est_new (mpicomm, conn_u, 0, NULL, NULL);
	p4est_i = p4est_new (mpicomm, conn_i, 0, NULL, NULL);

	/* Refine the forest recursively in parallel. */
	recursive = 1;
	p4est_refine (p4est1, recursive, refine_fn1, NULL);
	p4est_refine (p4est2, recursive, refine_fn2, NULL);

	/* Union operation */
	p4est_union(p4est1, p4est2, p4est_u);

	/* Intersection operation */
	p4est_intersection(p4est1, p4est2, p4est_i);

	/* Write the forest to disk for visualization */
	p4est_vtk_write_file (p4est1, NULL, P4EST_STRING "_setop3_cat1");
	p4est_vtk_write_file (p4est2, NULL, P4EST_STRING "_setop3_cat2");
	p4est_vtk_write_file (p4est_u, NULL, P4EST_STRING "_setop3_cats_union");
	p4est_vtk_write_file (p4est_i, NULL, P4EST_STRING "_setop3_cats_intersection");

	/* Destroy the p4est and the connectivity structure. */
	p4est_destroy (p4est1);
	p4est_destroy (p4est2);
	p4est_destroy (p4est_u);
	p4est_destroy (p4est_i);
	p4est_connectivity_destroy (conn1);
	p4est_connectivity_destroy (conn2);
	p4est_connectivity_destroy (conn_u);
	p4est_connectivity_destroy (conn_i);

	sc_finalize ();
	mpiret = sc_MPI_Finalize ();
	SC_CHECK_MPI (mpiret);
	return 0;
}

