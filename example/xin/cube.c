/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*
 * Usage: p4est_timings <configuration> <level>
 *        possible configurations:
 *        o unit      Refinement on the unit square.
 *        o periodic  Refinement on the unit square with periodic b.c.
 *        o three     Refinement on a forest with three trees.
 *        o moebius   Refinement on a 5-tree Moebius band.
 *        o star      Refinement on a 6-tree star shaped domain.
 *
 * Usage: p8est_timings <configuration> <level>
 *        possible configurations:
 *        o unit      Refinement on the unit cube.
 *        o periodic  Refinement on the unit cube with all-periodic b.c.
 *        o rotwrap   Refinement on the unit cube with weird periodic b.c.
 *        o twocubes  Refinement on a forest with two trees.
 *        o rotcubes  Refinement on a forest with six rotated trees.
 *        o shell     Refinement on a 24-tree spherical shell.
 */

#ifndef P4_TO_P8
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_nodes.h>
#include <p4est_vtk.h>
#include <p4est_lnodes.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_nodes.h>
#include <p8est_vtk.h>
#include <p8est_lnodes.h>
#endif

#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

#include "api.h"
#include <sys/time.h>

typedef struct
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;



static int
refine_fractal (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  int                 qid;

  if ((int) q->level >= refine_level) {
    return 0;
  }
  if ((int) q->level < refine_level - level_shift) {
    return 1;
  }

  qid = p4est_quadrant_child_id (q);
  return (qid == 0 || qid == 3
#ifdef P4_TO_P8
          || qid == 5 || qid == 6
#endif
    );
}

int
coarsen_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * q[])
{

  int i, n; 
  int offsj, offsi, offsk;
  int unit_len, tilelen;

  n = 4;
#ifdef P4_TO_P8
  n = 8;
#endif
  for (i=0; i<n; i++){
    if (q[i] == NULL) {
      return 0;
    }
  }

  if ((int) q[0]->level < min_level) {
    return 0;
  }
  /*
  unit_len = P4EST_QUADRANT_LEN (refine_level);
  tilelen = 1 << (refine_level - q[0]->level);
  offsi = (q[0]->x)/unit_len ;       
  offsj = (q[0]->y)/unit_len ;       
#ifdef P4_TO_P8
  offsk = (q[0]->z)/unit_len ;       
#endif
#ifdef P4_TO_P8
          printf ("coarsen ****  (x: %d y: %d z: %d) level %d tilelen: %d\n", offsi, offsj, offsk, q[0]->level, tilelen);
#endif
  */ 

  if (  q[0]->p.user_int == UNMARKED && q[1]->p.user_int == UNMARKED  
     && q[2]->p.user_int == UNMARKED && q[3]->p.user_int == UNMARKED
#ifdef P4_TO_P8
     && q[4]->p.user_int == UNMARKED && q[5]->p.user_int == UNMARKED
     && q[6]->p.user_int == UNMARKED && q[7]->p.user_int == UNMARKED 
#endif
  ){
    
    //printf("### coarsen function level: %d data: %d, %d, %d, %d\n ", 
    //  q[0]->level, q[0]->p.user_int, q[1]->p.user_int, q[2]->p.user_int, q[3]->p.user_int);
    //printf("following quadrants are going to be coarsened!!! return 1\n");
    return 1;
  }

#ifdef P4_TO_P8
  //printf ("coarsen **** (x: %d y: %d z: %d) level %d tilelen: %d\n", offsi, offsj, offsk, q[0]->level, tilelen);
#endif
  //printf("one quadrant is 1, return 0 !!! \n");
  return 0;
}

int p4est_diff_aafn (p4est_t *p4est_in1, p4est_t *p4est_in2, p4est_t *p4est,
                                       p4est_topidx_t which_tree,
                                       p4est_quadrant_t *quadrant)
{
        p4est_tree_t *tree1, *tree2;
        sc_array_t *quads1, *quads2;
        int q1_index, q2_index;
        p4est_quadrant_t *q1, *q2;
        int q1_match = 0;
        int tilelen, unit_len;

        if ( quadrant->level > refine_level ) return 0; 

        tilelen = 1 << (refine_level - quadrant->level);
        unit_len = P4EST_QUADRANT_LEN (refine_level);
        //printf ("**** in main set operation x: %d y: %d level %d tilelen: %d \n", quadrant->x/unit_len, quadrant->y/unit_len, quadrant->level, tilelen);
#ifdef P4_TO_P8
        //printf ("z: %d \n", quadrant->z/unit_len);
#endif
        // 1st input tree
        tree1 = p4est_tree_array_index(p4est_in1->trees, 0);
        quads1 = &(tree1->quadrants);
        q1_index = sc_array_bsearch(quads1, quadrant, p4est_quadrant_compare);
        if (q1_index != -1) { // quadrant match in 1st input tree
                //printf("found on in 1st input tree (((())))))\n");
                q1_match = 1;
                q1 = p4est_quadrant_array_index(quads1, q1_index);
                if (!q1->p.user_int) { // quadrant is "unmarked"; we're done
                        quadrant->p.user_int = UNMARKED;
                        return 0;
                }
        }
        // 2nd input tree
        tree2 = p4est_tree_array_index(p4est_in2->trees, 0);
        quads2 = &(tree2->quadrants);
        q2_index = sc_array_bsearch(quads2, quadrant, p4est_quadrant_compare);
        if (q2_index == -1 && !q1_match) { // no match on either input
                if (quadrant->p.user_int == GRAY) { // "subtree" case; refine
                        return -1;
                } else { // no match and NOT a subtree case; refine
                        return 1;
                }
        } else if (q1_match && q2_index == -1) { // "marked" match on 1st input, only
                if (quadrant->p.user_int == GRAY) { // in "subtree" case; done
                        quadrant->p.user_int = MARKED;  // tree2 is 0 >> 1 & 0 => 1 
                        return 0;
                } else { // entering "subtree" case; refine
                        printf(" tree1 is marked, but tree2 needs to be refined!!! enter sub tree !!!!!!!\n ");
                        return -1;
                }
        } else { // match in 2nd input tree
                q2 = p4est_quadrant_array_index(quads2, q2_index);
                if (q2->p.user_int) { // quadrant is "marked"; we're done
                        quadrant->p.user_int = UNMARKED;  // * & 1 => 0
                        return 0;
                } else if (q1_match && !q2->p.user_int) { // corresponding matches in both inputs, so this quadrant is "marked"
                        quadrant->p.user_int = MARKED;
                        return 0;
                } else if (!q1_match && !q2->p.user_int) { // no match in 1st input; matched "marked" in 2nd
                        if (quadrant->p.user_int == GRAY) { // in "subtree" case
                                printf(" tree2 is unmarked, tree1 had been refined!!! \n ");
                                quadrant->p.user_int = MARKED;
                                return 0;
                        } else { // entering "subtree" case; refine
                                return -1;
                        }
                } else { // shouldn't enter this clause
                        printf("Unexpected case triggered in p4est_diff_refine_fn\n");
                        return 0;
                }
        }
}

int
refine_sphere1 (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  int                 qid;
  int                 tilelen;
  int                 i, j, k;
  int                 markedCount;
  double                 offsi, offsj, offsk;
  double                 tx, ty, tz, temp;
  p4est_qcoord_t      unit_len, current_len;

  if ((int) q->level > refine_level) {
    return 0;
  }
  if ((int) q->level < min_level) {
    return 1;
  }
  tilelen = 1 << (refine_level - q->level);       /* Pixel size of quadrant */
  unit_len = P4EST_QUADRANT_LEN (refine_level);
  current_len = P4EST_QUADRANT_LEN (q->level); 
  offsi = (double)(q->x) / unit_len + 0.5;       /* Pixel x offset */
  offsj = (double)(q->y) / unit_len + 0.5;       /* Pixel y offset */
#ifdef P4_TO_P8
  offsk = (double)(q->z) / unit_len + 0.5;       /* Pixel z offset */
#endif
  
  P4EST_ASSERT (offsi >= 0 && offsj >= 0 && offsk >= 0);
  markedCount = 0; 
#ifdef P4_TO_P8
  for (k = 0; k < tilelen; ++k) {
#endif
    for (j = 0; j < tilelen; ++j) {
      for (i = 0; i < tilelen; ++i) {
        tx = offsi + i;
        ty = offsj + j;
#ifdef P4_TO_P8
        tz = offsk + k; 
#endif
        //tx++; ty++; tz++;
        temp = (tx-center)*(tx-center) + (ty-center)*(ty-center)
#ifdef P4_TO_P8
             + (tz-center_z)*(tz-center_z)
#endif
;
#ifdef P4_TO_P8
        //printf(" tx %f ty %f tz %f center %d radix %d \n ", tx, ty, tz, center, radix1);
#endif
        if ( temp <= radix1*radix1 ){
          markedCount += 1;
          //return 1;
        }
      }
    }
#ifdef P4_TO_P8
  }
#endif
#ifdef P4_TO_P8
#endif
  //printf ("x: %d y: %d z: %d markedCount: %d level %d tilelen: %d\n", q->x, q->y, q->z, markedCount, q->level, tilelen);
  //printf ("quadrant len %d tilelen*unit_len: %d\n", current_len, tilelen*unit_len);
  if (!markedCount) {
          q->p.user_int = UNMARKED;
          return 0;
  } else if (tilelen*tilelen*tilelen == 1 && markedCount == 1) {
  //printf("*************** offset offx: %f offy: %f  len: %d count %d temp %f\n", offsi-0.5, offsj-0.5, tilelen, markedCount, temp );
          q->p.user_int = MARKED;
          return 0;
  } else {
          //q->p.user_int = UNMARKED;
          return 1;
  }
  return 0;
}


int
refine_sphere2 (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  int                 qid;
  int                 tilelen;
  double                 offsi, offsj, offsk;
  int                 i, j, k;
  double                  tx, ty, tz, temp;
  int                 markedCount;
  p4est_qcoord_t      unit_len, offset;

  if ((int) q->level > refine_level) {
    return 0;
  }
  if ((int) q->level < min_level) {
    return 1;
  }

  tilelen = 1 << (refine_level - q->level);       /* Pixel size of quadrant */
  unit_len = P4EST_QUADRANT_LEN (refine_level);
  offsi = (double)(q->x) / unit_len + 0.5;       /* Pixel x offset */
  offsj = (double)(q->y) / unit_len + 0.5;       /* Pixel y offset */
#ifdef P4_TO_P8
  offsk = (double)(q->z) / unit_len + 0.5;       /* Pixel z offset */
#endif

  P4EST_ASSERT (offsi >= 0 && offsj >= 0 && offsk >= 0);
  markedCount = 0; 
#ifdef P4_TO_P8
  for (k = 0; k < tilelen; ++k) {
#endif
    for (j = 0; j < tilelen; ++j) {
      for (i = 0; i < tilelen; ++i) {
        tx = offsi + i;
        ty = offsj + j;
#ifdef P4_TO_P8
        tz = offsk + k; 
#endif
        //tx++; ty++; tz++;
        temp = (tx-center)*(tx-center) + (ty-center)*(ty-center)
#ifdef P4_TO_P8
             + (tz-center_z)*(tz-center_z)
#endif
;
         if ( temp <= radix2*radix2 ){
          markedCount += 1;
          //return 1;
        }
      }
    }
#ifdef P4_TO_P8
  }
#endif
  if (!markedCount) {
          q->p.user_int = UNMARKED;
          return 0;
  } else if (tilelen*tilelen*tilelen == 1 && markedCount == 1) {
  //printf("*************** offset offx: %f offy: %f  len: %d count %d temp %f\n", offsi-0.5, offsj-0.5, tilelen, markedCount, temp );
          q->p.user_int = MARKED;
          return 0;
  } else {
          //q->p.user_int = UNMARKED;
          return 1;
  }
  return 0;
}


int
refine_sphere3 (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  int                 qid;
  int                 tilelen;
  int                 i, j, k;
  int                 markedCount;
  double                 offsi, offsj, offsk;
  double                 tx, ty, tz, temp;
  p4est_qcoord_t      unit_len, current_len;

  if ((int) q->level > refine_level) {
    return 0;
  }
  if ((int) q->level < min_level) {
    return 1;
  }
  tilelen = 1 << (refine_level - q->level);       /* Pixel size of quadrant */
  unit_len = P4EST_QUADRANT_LEN (refine_level);
  current_len = P4EST_QUADRANT_LEN (q->level); 
  offsi = (double)(q->x) / unit_len + 0.5;       /* Pixel x offset */
  offsj = (double)(q->y) / unit_len + 0.5;       /* Pixel y offset */
#ifdef P4_TO_P8
  offsk = (double)(q->z) / unit_len + 0.5;       /* Pixel z offset */
#endif
  
  P4EST_ASSERT (offsi >= 0 && offsj >= 0 && offsk >= 0);
  markedCount = 0; 
#ifdef P4_TO_P8
  for (k = 0; k < tilelen; ++k) {
#endif
    for (j = 0; j < tilelen; ++j) {
      for (i = 0; i < tilelen; ++i) {
        tx = offsi + i;
        ty = offsj + j;
#ifdef P4_TO_P8
        tz = offsk + k; 
#endif
        //tx++; ty++; tz++;
        temp = (tx-center)*(tx-center) + (ty-center)*(ty-center)
#ifdef P4_TO_P8
             + (tz-center)*(tz-center)
#endif
;
        //if ( temp >= radix1*radix1 && temp <= radix2*radix2 ){
#ifdef P4_TO_P8
        //printf(" tx %f ty %f tz %f center %d radix %d \n ", tx, ty, tz, center, radix1);
#endif
        if ( temp <= radix2*radix2 ){
          markedCount += 1;
          //return 1;
        }
      }
    }
#ifdef P4_TO_P8
  }
#endif
#ifdef P4_TO_P8
#endif
  //printf ("x: %d y: %d z: %d markedCount: %d level %d tilelen: %d\n", q->x, q->y, q->z, markedCount, q->level, tilelen);
  //printf ("quadrant len %d tilelen*unit_len: %d\n", current_len, tilelen*unit_len);
  if (!markedCount) {
          q->p.user_int = UNMARKED;
          return 0;
  } else if (tilelen*tilelen*tilelen == 1 && markedCount == 1) {
  //printf("*************** offset offx: %f offy: %f  len: %d count %d temp %f\n", offsi-0.5, offsj-0.5, tilelen, markedCount, temp );
          q->p.user_int = MARKED;
          return 0;
  } else {
          //q->p.user_int = UNMARKED;
          return 1;
  }
  return 0;
}


int
main (int argc, char **argv)
{
  int                 i;
  int                 mpiret;
  int                 wrongusage;
  unsigned            crc, gcrc;
  const char         *config_name;
  const char         *load_name;
  p4est_locidx_t     *quadrant_counts;
  p4est_gloidx_t      count_refined, count_balanced;
  p4est_gloidx_t      prev_quadrant, next_quadrant;
  p4est_gloidx_t      global_shipped;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est1;
  p4est_t            *p4est2;
  p4est_t            *p4est_out;
  p4est_nodes_t      *nodes = NULL;
  p4est_ghost_t      *ghost;
  p4est_lnodes_t     *lnodes;
  sc_flopinfo_t       fi, snapshot;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  sc_options_t       *opt;
  int                 overlap;
  int                 subtree;
  int                 borders;
  int                 max_ranges;
  int                 use_ranges, use_ranges_notify, use_balance_verify;
  int                 oldschool, generate;
  int                 first_argc;
  int                 test_multiple_orders;
  int                 skip_nodes, skip_lnodes;
  int                 repartition_lnodes;

  //*********
  struct timeval t1, t2, t3, t4, t5;
  double elapsedTime;
  //*********

  /* initialize MPI and p4est internals */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpi->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
#ifndef P4EST_ENABLE_DEBUG
  sc_set_log_defaults (NULL, NULL, SC_LP_STATISTICS);
#endif
  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  P4EST_GLOBAL_PRODUCTIONF ("Size of %dtant: %lld bytes\n", P4EST_DIM,
                            (long long) sizeof (p4est_quadrant_t));


  /* get command line argument: maximum refinement level */
  //level_shift = 2;
  //refine_level = 5;

  /* print general setup information */
  //P4EST_GLOBAL_STATISTICSF
  //  ("Processors %d configuration %s level %d shift %d\n", mpi->mpisize,
  //   config_name, refine_level, level_shift);
  printf (" @@@@@@@@@@ test coordinates eighth %d %d %d \n", eighth, center, center_z);

  /* start overall timing */
  mpiret = sc_MPI_Barrier (mpi->mpicomm);
  //SC_CHECK_MPI (mpiret);
  //sc_flops_start (&fi);

  /* create connectivity and forest structures */
#ifndef P4_TO_P8
    
      connectivity = p4est_connectivity_new_unitsquare ();
#else
    
      printf(" created unit cube from connectivity!!! \n");
      connectivity = p8est_connectivity_new_unitcube ();
      //connectivity = p8est_connectivity_new_twocubes ();
#endif

    
  p4est1 = p4est_new_ext (mpi->mpicomm, connectivity,
                             0, min_level, 1, 0, NULL, NULL);
  
  p4est2 = p4est_copy(p4est1, 1); 
  p4est_out = p4est_copy(p4est1, 1); 

  //p4est_addtree( p4est, p4est1, 1 );

  max_ranges = -1;
  p4est1->inspect = P4EST_ALLOC_ZERO (p4est_inspect_t, 1);
  p4est1->inspect->use_balance_ranges = use_ranges;
  p4est1->inspect->use_balance_ranges_notify = use_ranges_notify;
  p4est1->inspect->use_balance_verify = use_balance_verify;
  p4est1->inspect->balance_max_ranges = max_ranges;
  P4EST_GLOBAL_STATISTICSF
    ("Balance: new overlap %d new subtree %d borders %d\n", overlap,
     (overlap && subtree), (overlap && borders));
  quadrant_counts = P4EST_ALLOC (p4est_locidx_t, p4est1->mpisize);

#ifdef P4EST_TIMINGS_VTK
  //p4est_vtk_write_file (p4est1, NULL, "tree1_mynew");
  //p4est_vtk_write_file (p4est2, NULL, "tree2_mynew");
#endif

  /* time refine */
  sc_flops_snap (&fi, &snapshot);
  p4est_refine (p4est1, 1, refine_sphere1, NULL);
  p4est_refine (p4est2, 1, refine_sphere2, NULL);
  //p4est_refine (p4est2, 1, refine_bowl, NULL);
  //p4est_refine (p4est1, 1, refine_fn1, NULL);
  //p4est_refine (p4est2, 1, refine_fn2, NULL);
  sc_flops_shot (&fi, &snapshot);

gettimeofday(&t1, NULL);
  //p4est_intersection(p4est1, p4est2, p4est_out);
  //p4est_union(p4est1, p4est2, p4est_out);
  p4est_diff(p4est2, p4est1, p4est_out, p4est_diff_aafn, coarsen_fn);
gettimeofday(&t2, NULL);


  p4est_remove(p4est_out);
  //p4est_remove(p4est1);
  //p4est_remove(p4est2);
  //sc_stats_set1 (&stats[TIMINGS_REFINE], snapshot.iwtime, "Refine");
  printf("finish operations here !!!\n");

#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est1, NULL, "tree1_refined");
  p4est_vtk_write_file (p4est2, NULL, "tree2_refined");
  p4est_vtk_write_file (p4est_out, NULL, "diff_refined");
#endif
  count_refined = p4est1->global_num_quadrants;

  //p4est_coarsen (p4est1, 1, coarsen_fn, NULL);
  //p4est_vtk_write_file (p4est1, NULL, "tree1_coarsed");

  /* time balance */
  sc_flops_snap (&fi, &snapshot);
  //p4est_balance (p4est1, P4EST_CONNECT_FULL, NULL);
  sc_flops_shot (&fi, &snapshot);
#ifdef P4EST_TIMINGS_VTK
  //p4est_vtk_write_file (p4est1, NULL, "timings_balanced");
#endif

  /* calculate and print timings */
  //sc_stats_compute (mpi->mpicomm, TIMINGS_NUM_STATS, stats);
  //sc_stats_print (p4est_package_id, SC_LP_STATISTICS,
  //                TIMINGS_NUM_STATS, stats, 1, 1);

  /* destroy the p4est and its connectivity structure */
  P4EST_FREE (quadrant_counts);
  P4EST_FREE (p4est1->inspect);
  p4est_destroy (p4est1);
  p4est_destroy (p4est2);
  p4est_destroy (p4est_out);
  p4est_connectivity_destroy (connectivity);

  elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
  printf(" operation p4est elapsedTime %f  ms.\n", elapsedTime);

  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
