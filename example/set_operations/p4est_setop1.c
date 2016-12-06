/** \file p4est_setop1.c\
 *
 */

#include <p4est_vtk.h>
#include "rect_horizontal1.h"
#include "rect_vertical1.h"
#include "rectangle_union.h"
#include <stdio.h>

#define P4EST_STEP1_PATTERN_LEVEL 5
/** The dimension of the image data. */
#define P4EST_STEP1_PATTERN_LENGTH (1 << P4EST_STEP1_PATTERN_LEVEL)
static const int plv = P4EST_STEP1_PATTERN_LEVEL;
static const int ple = P4EST_STEP1_PATTERN_LENGTH;

static int refine_fn_h(p4est_t * p4est, p4est_topidx_t which_tree,
		p4est_quadrant_t * quadrant)
{
  int                 tilelen;
  int                 offsi, offsj;
  int                 i, j;
  const char         *d;
  unsigned char       p[3];

  P4EST_ASSERT (which_tree == 0);
  /* We do not want to refine deeper than a given maximum level. */
  if (quadrant->level > plv) {
    return 0;
  }

  tilelen = 1 << (plv - quadrant->level);       /* Pixel size of quadrant */
  // P4EST_QUADRANT_LEN (plv) would always return 1?
  offsi = quadrant->x / P4EST_QUADRANT_LEN (plv);       /* Pixel x offset */
  offsj = quadrant->y / P4EST_QUADRANT_LEN (plv);       /* Pixel y offset */
  P4EST_ASSERT (offsi >= 0 && offsj >= 0);
  for (j = 0; j < tilelen; ++j) {
    P4EST_ASSERT (offsj + j < ple);
    for (i = 0; i < tilelen; ++i) {
      P4EST_ASSERT (offsi + i < ple);
      d = recth_header_data + 4 * (ple * (ple - 1 - (offsj + j)) + (offsi + i)); // TODO: print out?
      RECTH_HEADER_PIXEL (d, p);
      P4EST_ASSERT (p[0] == p[1] && p[1] == p[2]);      /* Grayscale image */
      if (p[0] < 128) {
        return 1;
      }
    }
  }
  return 0;
}

static int refine_fn_v(p4est_t * p4est, p4est_topidx_t which_tree,
		p4est_quadrant_t * quadrant)
{
  int                 tilelen;
  int                 offsi, offsj;
  int                 i, j;
  const char         *d;
  unsigned char       p[3];

  P4EST_ASSERT (which_tree == 0);
  if (quadrant->level > plv) {
    return 0;
  }

  tilelen = 1 << (plv - quadrant->level);       /* Pixel size of quadrant */
  offsi = quadrant->x / P4EST_QUADRANT_LEN (plv);       /* Pixel x offset */
  offsj = quadrant->y / P4EST_QUADRANT_LEN (plv);       /* Pixel y offset */
//  printf("%d, %d, %d\n", offsi, offsj, quadrant->level);
  P4EST_ASSERT (offsi >= 0 && offsj >= 0);
  for (j = 0; j < tilelen; ++j) {
    P4EST_ASSERT (offsj + j < ple);
    for (i = 0; i < tilelen; ++i) {
      P4EST_ASSERT (offsi + i < ple);
      d = rectv_header_data + 4 * (ple * (ple - 1 - (offsj + j)) + (offsi + i));
      RECTV_HEADER_PIXEL (d, p);
      P4EST_ASSERT (p[0] == p[1] && p[1] == p[2]);      /* Grayscale image */
      if (p[0] < 128) {
        return 1;
      }
    }
  }
  return 0;
}

static int refine_fn_u(p4est_t * p4est, p4est_topidx_t which_tree,
		p4est_quadrant_t * quadrant)
{
  int                 tilelen;
  int                 offsi, offsj;
  int                 i, j;
  const char         *d;
  unsigned char       p[3];

  P4EST_ASSERT (which_tree == 0);
  if (quadrant->level > plv) {
    return 0;
  }

  tilelen = 1 << (plv - quadrant->level);       /* Pixel size of quadrant */
  offsi = quadrant->x / P4EST_QUADRANT_LEN (plv);       /* Pixel x offset */
  offsj = quadrant->y / P4EST_QUADRANT_LEN (plv);       /* Pixel y offset */
  P4EST_ASSERT (offsi >= 0 && offsj >= 0);
  for (j = 0; j < tilelen; ++j) {
    P4EST_ASSERT (offsj + j < ple);
    for (i = 0; i < tilelen; ++i) {
      P4EST_ASSERT (offsi + i < ple);
      d = rect_union_header_data + 4 * (ple * (ple - 1 - (offsj + j)) + (offsi + i));
      RECT_UNION_HEADER_PIXEL (d, p);
      P4EST_ASSERT (p[0] == p[1] && p[1] == p[2]);      /* Grayscale image */
      if (p[0] < 128) {
        return 1;
      }
    }
  }
  return 0;
}


int main (int argc, char **argv)
{
  int                 mpiret;
  int                 recursive, partforcoarsen, balance;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est_h, *p4est_v, *p4est_u;
  p4est_connectivity_t *conn_h, *conn_v, *conn_u;

  /* Initialize MPI; see sc_mpi.h. */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  /* Create single-quadtree forests */
  conn_h = p4est_connectivity_new_unitsquare();
  conn_v = p4est_connectivity_new_unitsquare();
  conn_u = p4est_connectivity_new_unitsquare();
  p4est_h = p4est_new (mpicomm, conn_h, 0, NULL, NULL);
  p4est_v = p4est_new (mpicomm, conn_v, 0, NULL, NULL);
  p4est_u = p4est_new (mpicomm, conn_u, 0, NULL, NULL);

  /* Refine the forest recursively in parallel.
   * The P4EST_ASSERT macro only activates with --enable-debug.
   * We check against the data dimensions in example/steps/hw32.h. */
  recursive = 1;
  p4est_refine (p4est_h, recursive, refine_fn_h, NULL);
  p4est_refine (p4est_v, recursive, refine_fn_v, NULL);
  p4est_refine (p4est_u, recursive, refine_fn_u, NULL);

  /* Write the forest to disk for visualization */
  p4est_vtk_write_file (p4est_h, NULL, P4EST_STRING "_setop1_rectangle_h");
  p4est_vtk_write_file (p4est_v, NULL, P4EST_STRING "_setop1_rectangle_v");
  p4est_vtk_write_file (p4est_u, NULL, P4EST_STRING "_setop1_rectangle_union");

  /* Destroy the p4est and the connectivity structure. */
  p4est_destroy (p4est_h);
  p4est_destroy (p4est_v);
  p4est_destroy (p4est_u);
  p4est_connectivity_destroy (conn_h);
  p4est_connectivity_destroy (conn_v);
  p4est_connectivity_destroy (conn_u);

  /* Verify that allocations internal to p4est and sc do not leak memory.
   * This should be called if sc_init () has been called earlier. */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
