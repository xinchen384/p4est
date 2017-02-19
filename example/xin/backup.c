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

 #define P4EST_TIMINGS_VTK 
/*
*/

typedef enum
{
  P4EST_CONFIG_NULL,
  P4EST_CONFIG_UNIT,
  P4EST_CONFIG_PERIODIC,
#ifndef P4_TO_P8
  P4EST_CONFIG_THREE,
  P4EST_CONFIG_MOEBIUS,
  P4EST_CONFIG_STAR
#else
  P4EST_CONFIG_ROTWRAP,
  P4EST_CONFIG_TWOCUBES,
  P4EST_CONFIG_ROTCUBES,
  P4EST_CONFIG_SHELL
#endif
}
timings_config_t;

/*
TIMINGS_BALANCE_A,
  TIMINGS_BALANCE_COMM,
  TIMINGS_BALANCE_B,
  TIMINGS_BALANCE_A_COUNT_IN,
  TIMINGS_BALANCE_A_COUNT_OUT,
  TIMINGS_BALANCE_COMM_SENT,
  TIMINGS_BALANCE_COMM_NZPEERS,
  TIMINGS_BALANCE_B_COUNT_IN,
  TIMINGS_BALANCE_B_COUNT_OUT,
  TIMINGS_BALANCE_RANGES,
  TIMINGS_BALANCE_NOTIFY,
  TIMINGS_BALANCE_NOTIFY_ALLGATHER,
  TIMINGS_BALANCE_A_ZERO_SENDS,
  TIMINGS_BALANCE_A_ZERO_RECEIVES,
  TIMINGS_BALANCE_B_ZERO_SENDS,
  TIMINGS_BALANCE_B_ZERO_RECEIVES,

TIMINGS_REBALANCE_A,
  TIMINGS_REBALANCE_COMM,
  TIMINGS_REBALANCE_B,
  TIMINGS_REBALANCE_A_COUNT_IN,
  TIMINGS_REBALANCE_A_COUNT_OUT,
  TIMINGS_REBALANCE_COMM_SENT,
  TIMINGS_REBALANCE_COMM_NZPEERS,
  TIMINGS_REBALANCE_B_COUNT_IN,
  TIMINGS_REBALANCE_B_COUNT_OUT,

*/
enum
{
  TIMINGS_REFINE,
  TIMINGS_BALANCE,
    TIMINGS_REBALANCE,
    TIMINGS_PARTITION,
  TIMINGS_GHOSTS,
  TIMINGS_NODES,
  TIMINGS_TRILINEAR_OBSOLETE,
  TIMINGS_REPARTITION,
  TIMINGS_LNODES,
  //TIMINGS_LNODES3,
  //TIMINGS_LNODES7,
  TIMINGS_NUM_STATS
};

typedef struct
{
  timings_config_t    config;
  int                 mpisize;
  int                 level;
  unsigned            checksum;
}
timings_regression_t;

typedef struct
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

#define P4EST_REFINE_LEVEL 6 
/** The dimension of the image data. */
#define P4EST_REFINE_LENGTH (1 << P4EST_REFINE_LEVEL)
static const int          refine_level = P4EST_REFINE_LEVEL;
static const int          level_shift = 2;
static const int          min_level = 2;
static const int          plen = P4EST_REFINE_LENGTH;

static const p4est_qcoord_t eighth = P4EST_QUADRANT_LEN (3);

static const p4est_qcoord_t center = P4EST_REFINE_LENGTH/2;
static const p4est_qcoord_t center_z = P4EST_REFINE_LENGTH;
static const p4est_qcoord_t radix1 = P4EST_REFINE_LENGTH/4;
static const p4est_qcoord_t radix2 = P4EST_REFINE_LENGTH*3/8;

static const p4est_qcoord_t p1 = P4EST_REFINE_LENGTH/4;
static const p4est_qcoord_t p2 = P4EST_REFINE_LENGTH*3/4;

static const int UNMARKED = 0;
static const int MARKED = 1;
static const int GRAY = -1;


static int
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

int in_the_mid(p4est_qcoord_t v1, p4est_qcoord_t v2, p4est_qcoord_t x){
  p4est_qcoord_t s1, s2;
  if ( v1 >= v2 ){
    s1 = v2;
    v2 = v1;
  } else{
    s1 = v1;
    s2 = v2;
  }
  if ( s1 <= x && x <= s2 ){
    return 1;
  } else{
    return 0;
  }
}

static int 
check_overlap(p4est_qcoord_t ox, p4est_qcoord_t oy,
#ifdef P4_TO_P8
p4est_qcoord_t oz, 
#endif
p4est_qcoord_t len, p4est_qcoord_t cx, p4est_qcoord_t cy,
#ifdef P4_TO_P8
p4est_qcoord_t cz, 
#endif
p4est_qcoord_t radix){
  p4est_qcoord_t sx, sy;
  p4est_qcoord_t temp_x, temp_y;
  p4est_qcoord_t temp;

#ifdef P4_TO_P8
  p4est_qcoord_t sz, temp_z;
  sz = oz + len/2;
#endif
  sx = ox + len/2;
  sy = oy + len/2;

  if (sx <= cx){
    temp_x = sx + len/2; 
    if (sy <= cy){
      temp_y = sy + len/2; 
#ifdef P4_TO_P8
      if (sz <= cz){
        temp_z = sz + len/2; 
      } else{
        temp_z = sz - len/2; 
      }
#endif
    } else{
      temp_y = sy - len/2; 
#ifdef P4_TO_P8
      if (sz <= cz){
        temp_z = sz + len/2; 
      } else{
        temp_z = sz - len/2; 
      }
#endif
    }
  } else{
    temp_x = sx - len/2; 
    if (sy <= cy){
      temp_y = sy + len/2; 
#ifdef P4_TO_P8
      if (sz <= cz){
        temp_z = sz + len/2; 
      } else{
        temp_z = sz - len/2; 
      }
#endif
    } else{
      temp_y = sy - len/2; 
#ifdef P4_TO_P8
      if (sz <= cz){
        temp_z = sz + len/2; 
      } else{
        temp_z = sz - len/2; 
      }
#endif
    }
  }

  temp = (temp_x-cx)*(temp_x-cx) + (temp_y-cy)*(temp_y-cy)
#ifdef P4_TO_P8
         + (temp_z-cz)*(temp_z-cz)
#endif
; 
  if (temp <= radix*radix)
    return 1; 
  else 
    return 0;
}

static int
refine_sphere_test (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  int                 qid;
  int                 tilelen;
  int                 offsi, offsj, offsk;
  int                 i, j, k;
  int                 tx, ty, tz, temp;
  int                 markedCount;
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

printf(" %d %d hieght is %d level %d \n", q->x, q->y,  current_len, q->level);
  
  if (check_overlap(q->x, q->y,
#ifdef P4_TO_P8
q->z, 
#endif
current_len, center, center,
#ifdef P4_TO_P8
center_z, 
#endif
radix1)){
    if (q->level >= refine_level){
      q->p.user_int = MARKED;
    }
    return 1;
  }
  q->p.user_int = UNMARKED;
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

static int
refine_fn1 (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  int                 qid;
  int                 tilelen;
  int                 i, j, k;
  int                 markedCount;
  p4est_qcoord_t      unit_len, current_len;
  double                 offsi, offsj, offsk;
  double                 tx, ty, tz, temp;


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
        if ( p1 <= tx && tx <= p2 && ty <= p2 ){
          markedCount += 1;
          //return 1;
        }
      }
    }
#ifdef P4_TO_P8
  }
#endif
  //printf ("x: %d y: %d markedCount: %d level %d tilelen: %d\n", q->x, q->y, markedCount, q->level, tilelen);
  //printf ("quadrant : %d p1 %d  p2  %d\n", P4EST_QUADRANT_LEN(0), p1, p2);
  if (!markedCount) {
          q->p.user_int = UNMARKED;
          return 0;
  } else if (tilelen*tilelen*tilelen == 1 && markedCount == 1) {
          q->p.user_int = MARKED;
#ifdef P4_TO_P8
          printf ("marked value tree1 (x: %f y: %f z: %f) level %d tilelen: %d\n", offsi-0.5, offsj-0.5, offsk-0.5, q->level, tilelen);
#endif
          return 0;
  } else {
          q->p.user_int = UNMARKED;
          return 1;
  }
  return 0;
}

static int
refine_fn2 (p4est_t * p4est, p4est_topidx_t which_tree,
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
        if ( p1 <= ty && ty <= p2 && tx <= p2 ){
          markedCount += 1;
          //return 1;
        }
      }
    }
#ifdef P4_TO_P8
  }
#endif
  //printf ("x: %d y: %d z: %d markedCount: %d level %d tilelen: %d\n", q->x, q->y, q->z, markedCount, q->level, tilelen);
  //printf ("quadrant len %d tilelen*unit_len: %d\n", current_len, tilelen*unit_len);
  if (!markedCount) {
          q->p.user_int = UNMARKED;
          return 0;
  } else if (tilelen*tilelen*tilelen == 1 && markedCount == 1) {
          q->p.user_int = MARKED;
#ifdef P4_TO_P8
          printf ("marked value tree2 (x: %f y: %f z:%f) level %d tilelen: %d\n", offsi-0.5, offsj-0.5, offsk-0.5, q->level, tilelen);
#endif
          return 0;
  } else {
          q->p.user_int = UNMARKED;
          return 1;
  }
  return 0;
}

static int
refine_bowl (p4est_t * p4est, p4est_topidx_t which_tree,
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
        //if ( temp >= radix1*radix1 && temp <= radix2*radix2 ){
#ifdef P4_TO_P8
        //printf(" tx %f ty %f tz %f center %d radix %d \n ", tx, ty, tz, center, radix1);
#endif
        if ( temp >= radix1*radix1 && temp <= radix2*radix2 ){
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


static int
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
        //if ( temp >= radix1*radix1 && temp <= radix2*radix2 ){
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

static int
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
             //+ (tz-center_z)*(tz-center_z);
        //if ( temp >= radix1*radix1 && temp <= radix2*radix2 ){
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

static int
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

/* *INDENT-OFF* */
static const timings_regression_t regression_oldschool[] =
{
#ifndef P4_TO_P8
  { P4EST_CONFIG_UNIT, 1, 10, 0x6e3e83c4U },
  { P4EST_CONFIG_UNIT, 1, 11, 0x334bc3deU },
  { P4EST_CONFIG_UNIT, 64, 14, 0xad908ce4U },
  { P4EST_CONFIG_UNIT, 256, 15, 0x9e7da646U },
  { P4EST_CONFIG_STAR, 1, 6, 0x14107b57U },
  { P4EST_CONFIG_STAR, 4, 6, 0x14107b57U },
  { P4EST_CONFIG_STAR, 52, 13, 0xc86c74d9U },
  { P4EST_CONFIG_STAR, 64, 13, 0xc86c74d9U },
#else
  { P4EST_CONFIG_UNIT, 1, 5, 0xe1ffa67bU },
  { P4EST_CONFIG_UNIT, 1, 6, 0x2cad814dU },
  { P4EST_CONFIG_UNIT, 3, 8, 0xeb252238U },
  { P4EST_CONFIG_PERIODIC, 1, 5, 0x99874fedU },
  { P4EST_CONFIG_PERIODIC, 2, 5, 0x575af6d5U },
  { P4EST_CONFIG_PERIODIC, 7, 6, 0xbc35524aU },
  { P4EST_CONFIG_ROTWRAP, 2, 6, 0x372f7402U },
  { P4EST_CONFIG_ROTWRAP, 7, 6, 0xa2f1ee48U },
  { P4EST_CONFIG_TWOCUBES, 5, 6, 0xa8b1f54eU },
  { P4EST_CONFIG_TWOCUBES, 8, 5, 0x98d3579dU },
  { P4EST_CONFIG_ROTCUBES, 1, 5, 0x404e4aa8U },
  { P4EST_CONFIG_ROTCUBES, 7, 6, 0x4c381706U },
  { P4EST_CONFIG_SHELL, 1, 4, 0x8c56f159U },
  { P4EST_CONFIG_SHELL, 3, 5, 0xafbc4f8cU },
  { P4EST_CONFIG_SHELL, 5, 6, 0xf6d9efb8U },
#endif
  { P4EST_CONFIG_NULL, 0, 0, 0 }
};

static const timings_regression_t regression_latest[] =
{
#ifndef P4_TO_P8
  { P4EST_CONFIG_UNIT, 1, 10, 0x6e3e83c4U },
  { P4EST_CONFIG_UNIT, 1, 11, 0x334bc3deU },
  { P4EST_CONFIG_UNIT, 64, 14, 0xad908ce4U },
  { P4EST_CONFIG_UNIT, 256, 15, 0x9e7da646U },
  { P4EST_CONFIG_STAR, 1, 6, 0x14107b57U },
  { P4EST_CONFIG_STAR, 4, 6, 0x14107b57U },
  { P4EST_CONFIG_STAR, 52, 13, 0xc86c74d9U },
  { P4EST_CONFIG_STAR, 64, 13, 0xc86c74d9U },
#else
  { P4EST_CONFIG_UNIT, 1, 5, 0xc2012e84U },
  { P4EST_CONFIG_UNIT, 1, 6, 0x2cad814dU },
  { P4EST_CONFIG_UNIT, 3, 8, 0xeb252238U },
  { P4EST_CONFIG_PERIODIC, 1, 5, 0x2776c9b7U },
  { P4EST_CONFIG_PERIODIC, 2, 5, 0x2776c9b7U },
  { P4EST_CONFIG_PERIODIC, 7, 6, 0x4f281079U },
  { P4EST_CONFIG_ROTWRAP, 2, 6, 0x372f7402U },
  { P4EST_CONFIG_ROTWRAP, 7, 6, 0x372f7402U },
  { P4EST_CONFIG_TWOCUBES, 5, 6, 0xa8b1f54eU },
  { P4EST_CONFIG_TWOCUBES, 8, 5, 0x0aad11d0U },
  { P4EST_CONFIG_ROTCUBES, 1, 5, 0x404e4aa8U },
  { P4EST_CONFIG_ROTCUBES, 7, 6, 0x4c381706U },
  { P4EST_CONFIG_SHELL, 1, 4, 0x8c56f159U },
  { P4EST_CONFIG_SHELL, 3, 5, 0xafbc4f8cU },
  { P4EST_CONFIG_SHELL, 5, 6, 0xf6d9efb8U },
#endif
  { P4EST_CONFIG_NULL, 0, 0, 0 }
};
/* *INDENT-ON* */

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
  const timings_regression_t *r, *regression;
  timings_config_t    config;
  sc_statinfo_t       stats[TIMINGS_NUM_STATS];
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
  printf (" @@@@@@@@@@ test coordinates eighth %d %d %d \n", eighth, center, center_z);


  /* get command line argument: maximum refinement level */
  //level_shift = 2;
  //refine_level = 5;

  /* print general setup information */
  P4EST_GLOBAL_STATISTICSF
    ("Processors %d configuration %s level %d shift %d\n", mpi->mpisize,
     config_name, refine_level, level_shift);

  /* start overall timing */
  mpiret = sc_MPI_Barrier (mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  sc_flops_start (&fi);

  /* create connectivity and forest structures */
  regression = NULL;
#ifndef P4_TO_P8
    
      connectivity = p4est_connectivity_new_unitsquare ();
#else
    
      printf(" jjjjjjjjj unit cube from connectivity!!! \n");
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
  p4est_refine (p4est2, 1, refine_sphere3, NULL);
  //p4est_refine (p4est2, 1, refine_bowl, NULL);
  //p4est_refine (p4est1, 1, refine_fn1, NULL);
  //p4est_refine (p4est2, 1, refine_fn2, NULL);
  sc_flops_shot (&fi, &snapshot);
  //p4est_intersection(p4est1, p4est2, p4est_out);
  p4est_union(p4est1, p4est2, p4est_out);
  //p4est_diff(p4est2, p4est1, p4est_out, p4est_diff_aafn, coarsen_fn);

  p4est_remove(p4est_out);
  //p4est_remove(p4est1);
  //p4est_remove(p4est2);

  //sc_stats_set1 (&stats[TIMINGS_REFINE], snapshot.iwtime, "Refine");


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
  p4est_balance (p4est1, P4EST_CONNECT_FULL, NULL);
  sc_flops_shot (&fi, &snapshot);

  //sc_stats_set1 (&stats[TIMINGS_BALANCE], snapshot.iwtime, "Balance");
  //printf ("ddddd try balance %f \n", snapshot.iwtime );



#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est1, NULL, "timings_balanced");
#endif
  count_balanced = p4est1->global_num_quadrants;
  crc = p4est_checksum (p4est1);

  /* time rebalance - is a noop on the tree */
  sc_flops_snap (&fi, &snapshot);
  p4est_balance (p4est1, P4EST_CONNECT_FULL, NULL);
  sc_flops_shot (&fi, &snapshot);
  //sc_stats_set1 (&stats[TIMINGS_REBALANCE], snapshot.iwtime, "Rebalance");

  
  P4EST_ASSERT (count_balanced == p4est1->global_num_quadrants);
  P4EST_ASSERT (crc == p4est_checksum (p4est1));

  /* time a uniform partition */
  sc_flops_snap (&fi, &snapshot);
  p4est_partition (p4est1, 0, NULL);
  sc_flops_shot (&fi, &snapshot);
  //sc_stats_set1 (&stats[TIMINGS_PARTITION], snapshot.iwtime, "Partition");
#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est1, NULL, "timings_partitioned");
#endif
  P4EST_ASSERT (crc == p4est_checksum (p4est1));

  /* time building the ghost layer */
  sc_flops_snap (&fi, &snapshot);
  ghost = p4est_ghost_new (p4est1, P4EST_CONNECT_FULL);
  sc_flops_shot (&fi, &snapshot);
  //sc_stats_set1 (&stats[TIMINGS_GHOSTS], snapshot.iwtime, "Ghost layer");
  gcrc = p4est_ghost_checksum (p4est1, ghost);

  /* time the node numbering */
  if (!skip_nodes) {
    sc_flops_snap (&fi, &snapshot);
    nodes = p4est_nodes_new (p4est1, ghost);
    sc_flops_shot (&fi, &snapshot);
    //sc_stats_set1 (&stats[TIMINGS_NODES], snapshot.iwtime, "Nodes");
  }
  else {
    //sc_stats_set1 (&stats[TIMINGS_NODES], 0., "Nodes");
  }
  /* set this anyway so the output format is dimension independent */
  //sc_stats_set1 (&stats[TIMINGS_TRILINEAR_OBSOLETE], 0., "Unused");

  if (!skip_nodes) {
    p4est_nodes_destroy (nodes);
  }

  p4est_ghost_destroy (ghost);

  /* time a partition with a shift of all elements by one processor */
  for (i = 0, next_quadrant = 0; i < p4est1->mpisize; ++i) {
    prev_quadrant = next_quadrant;
    next_quadrant = (p4est1->global_num_quadrants * (i + 1)) / p4est1->mpisize;
    quadrant_counts[i] = (p4est_locidx_t) (next_quadrant - prev_quadrant);
  }
  if (p4est1->mpisize > 1) {
    quadrant_counts[0] += quadrant_counts[p4est1->mpisize - 1];  /* same type */
    quadrant_counts[p4est1->mpisize - 1] = 0;
  }

  sc_flops_snap (&fi, &snapshot);
  global_shipped = p4est_partition_given (p4est1, quadrant_counts);
  sc_flops_shot (&fi, &snapshot);
  //sc_stats_set1 (&stats[TIMINGS_REPARTITION], snapshot.iwtime, "Repartition");

  P4EST_GLOBAL_PRODUCTIONF
    ("Done " P4EST_STRING "_partition_given shipped %lld quadrants %.3g%%\n",
     (long long) global_shipped,
     global_shipped * 100. / p4est1->global_num_quadrants);
  P4EST_ASSERT (crc == p4est_checksum (p4est1));

  /* verify forest checksum */
  if (regression != NULL && mpi->mpirank == 0) {
    for (r = regression; r->config != P4EST_CONFIG_NULL; ++r) {
      if (r->config != config || r->mpisize != mpi->mpisize
          || r->level != refine_level)
        continue;
      SC_CHECK_ABORT (crc == r->checksum, "Checksum mismatch");
      P4EST_GLOBAL_INFO ("Checksum regression OK\n");
      break;
    }
  }


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


  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
