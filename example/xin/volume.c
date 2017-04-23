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

#include <p4est_to_p8est.h>

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

#include <sys/time.h>

#include "api.h"

int*** test_array1;
int*** test_array2;
int cube_len = -1;
int sk_level = 4;

int print_sign = 0;

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

void test_break(){}

void sort_p4est(p4est_t *p4est){
	sc_array_t *quads;
	p4est_tree_t *tree;
        tree = p4est_tree_array_index(p4est->trees, 0);
	quads = &(tree->quadrants);
        sc_array_sort(quads, p4est_quadrant_compare);
}

void check_data_fields(p4est_t *p4est){
        sc_array_t *quads;
	p4est_tree_t *tree;
        sc_array_t *array;
        p4est_quadrant_t *q;
        int i;
	char *s;
        int myc = 0;

        tree = p4est_tree_array_index(p4est->trees, 0);
	quads = &(tree->quadrants);
	array = quads;
        s = array->array;
        for (i=0; i<array->elem_count; i++){
          q = (p4est_quadrant_t *) s;
          printf("%d ", q->p.user_int); 
          s += array->elem_size;
          myc++;
        } 
        printf("\n quadrant counts: %d, data: \n", myc);
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
        //printf ("**** in main set operation x: %d y: %d z: %d level %d tilelen: %d \n", quadrant->x/unit_len, quadrant->y/unit_len, quadrant->z/unit_len, quadrant->level, tilelen);
        // 1st input tree
        tree1 = p4est_tree_array_index(p4est_in1->trees, 0);
        quads1 = &(tree1->quadrants);
        //if (print_sign == 1) 
        //  printf(" tree 1 element count: %zu, ele size: %zu\n", quads1->elem_count, quads1->elem_size);
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
        if (print_sign == 1)
          sc_array_sort(quads2, p4est_quadrant_compare);
        q2_index = sc_array_bsearch(quads2, quadrant, p4est_quadrant_compare);
        if (print_sign == 1){ 
          //printf(" tree 2 element count: %zu, ele size: %zu\n", quads2->elem_count, quads2->elem_size);
          //if (q2_index != -1)  printf("======>>>>>>>>>>>>>> found one in the second tree!!!! %d\n ", q2_index);
        }
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
                        //if (print_sign == 1)
                        //printf(" tree1 is marked, but tree2 needs to be refined!!! enter sub tree !!!!!!!\n ");
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
                                //if (print_sign == 1)
                                //printf(" tree2 is unmarked, tree1 had been refined!!! \n ");
                                quadrant->p.user_int = MARKED;
                                return 0;
                        } else { // entering "subtree" case; refine
                                //if (print_sign == 1)
                                //printf(" tree2 is unmarked, tree1 needs refined, sub sub\n ");
                                return -1;
                        }
                } else { // shouldn't enter this clause
                        printf("Unexpected case triggered in p4est_diff_refine_fn\n");
                        return 0;
                }
        }
}

int valid_id(int x, int y, int z){
  int res = 0;
  int full_len = 1 << refine_level;
  if ( x >= 0 && x < full_len 
    && y >= 0 && y < full_len
    && z >= 0 && z < full_len )
    res = 1;
  return res;
}

int
refine_offset (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  int                 qid;
  int                 tilelen;
  int                 i, j, k;
  int                 markedCount;
  double                 offsi, offsj, offsk;
  double                 tx, ty, tz, temp;
  // the voxel that resides in the sphere of tx,ty,taz(offset)
  double                 px, py, pz;  
  int ix, iy, iz;
  int full_len;
  p4est_qcoord_t      unit_len, current_len;

  if ((int) q->level > refine_level) {
    return 0;
  }
  if ((int) q->level < min_level) {
    return 1;
  }
  full_len = 1 << refine_level;
  tilelen = 1 << (refine_level - q->level);       /* Pixel size of quadrant */
  unit_len = P4EST_QUADRANT_LEN (refine_level);
  current_len = P4EST_QUADRANT_LEN (q->level); 
  offsi = (double)(q->x) / unit_len + 0.5;       /* Pixel x offset */
  offsj = (double)(q->y) / unit_len + 0.5;       /* Pixel y offset */
  offsk = (double)(q->z) / unit_len + 0.5;       /* Pixel z offset */
  
  P4EST_ASSERT (offsi >= 0 && offsj >= 0 && offsk >= 0);
  markedCount = 0; 
  for (k = 0; k < tilelen; ++k) {
    for (j = 0; j < tilelen; ++j) {
      for (i = 0; i < tilelen; ++i) {
        tx = offsi + i;
        ty = offsj + j;
        tz = offsk + k; 
        for (px = tx - offset; px <= tx + offset; px++)
        for (py = ty - offset; py <= ty + offset; py++)
        for (pz = tz - offset; pz <= tz + offset; pz++){
          ix = (int)(px-0.5);
          iy = (int)(py-0.5);
          iz = (int)(pz-0.5);
          if ( ix >= 0 && ix < full_len 
               && iy >= 0 && iy < full_len
               && iz >= 0 && iz < full_len ){
            if (test_array1[ix][iy][iz] == 1){
              temp = (tx-px)*(tx-px) + (ty-py)*(ty-py) + (tz-pz)*(tz-pz);
              if (temp <= offset * offset) markedCount += 1;
            }
          }
          /*
          */
        }
        //ix = (int)(tx-0.5);
        //iy = (int)(ty-0.5);
        //iz = (int)(tz-0.5);
        //if (test_array1[ix][iy][iz] == 1) markedCount += 1;
      }
    }
  }
  if (!markedCount) {
          q->p.user_int = UNMARKED;
          return 0;
  } else if (tilelen*tilelen*tilelen == 1 && markedCount == 1) {
          q->p.user_int = MARKED;
          return 0;
  } else {
          return 1;
  }
  return 0;
}

int
refine_full (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  
  if ((int) q->level >= sk_level){
    q->p.user_int = MARKED;
    return 0;
  }
  if ((int) q->level < min_level) {
    return 1;
  }
  return 1;
}

int
refine_test (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  int                 qid;
  int                 tilelen;
  int                 i, j, k;
  int                 markedCount, hardCount;
  double                 offsi, offsj, offsk;
  double                 tx, ty, tz, temp;
  int ix, iy, iz;
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
  offsk = (double)(q->z) / unit_len + 0.5;       /* Pixel z offset */
  
  P4EST_ASSERT (offsi >= 0 && offsj >= 0 && offsk >= 0);
  markedCount = 0; 
  for (k = 0; k < tilelen; ++k) {
    for (j = 0; j < tilelen; ++j) {
      for (i = 0; i < tilelen; ++i) {
        tx = offsi + i;
        ty = offsj + j;
        tz = offsk + k; 
        ix = (int)(offsi-0.5) + i;
        iy = (int)(offsj-0.5) + j;
        iz = (int)(offsk-0.5) + k;
        if ( (ix-center)*(ix-center) + (iy-center)*(iy-center) + (iz-center_z)*(iz-center_z) <= radix1*radix1 ){
          markedCount += 1;
        }
      }
    }
  }
  
  if (!markedCount) {
          q->p.user_int = UNMARKED;
          return 0;
  } else if (tilelen*tilelen*tilelen == 1 && markedCount == 1) {
          q->p.user_int = MARKED;
          return 0;
  } else {
          return 1;
  }
  return 0;
}

int
refine_sphere1 (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  int                 qid;
  int                 tilelen;
  int                 i, j, k;
  int                 markedCount, hardCount;
  double                 offsi, offsj, offsk;
  double                 tx, ty, tz, temp;
  int ix, iy, iz;
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
  offsk = (double)(q->z) / unit_len + 0.5;       /* Pixel z offset */
  
  P4EST_ASSERT (offsi >= 0 && offsj >= 0 && offsk >= 0);
  markedCount = 0; 
  hardCount = 0;
  for (k = 0; k < tilelen; ++k) {
    for (j = 0; j < tilelen; ++j) {
      for (i = 0; i < tilelen; ++i) {
        tx = offsi + i;
        ty = offsj + j;
        tz = offsk + k; 
        ix = (int)(offsi-0.5) + i;
        iy = (int)(offsj-0.5) + j;
        iz = (int)(offsk-0.5) + k;
        if (test_array1[ix][iy][iz] == 1||test_array1[ix][iy][iz] == 3) markedCount += 1;
        if (test_array1[ix][iy][iz] == 2) hardCount += 1;
      }
    }
  }
  if (hardCount == tilelen*tilelen*tilelen){
    q->p.user_int = 2;
    return 0;
  }
  if (!markedCount) {
          q->p.user_int = UNMARKED;
          return 0;
  } else if (tilelen*tilelen*tilelen == 1 && markedCount == 1) {
          q->p.user_int = MARKED;
          return 0;
  } else {
          return 1;
  }
  return 0;
}

int check_data(int*** data_array, int si, int sj, int sk, int level){
  int times;   
  int ti, tj, tk;
  int i, j, k; 
  int t1, t2, t3;
  if ( level > refine_level ){
    printf("****** need to handle a higher refine level!!!!*****\n");
    return 0;
  }
  times = 1 << (refine_level - level);
  ti = si*times;
  tj = sj*times;
  tk = sk*times;
  for (i=0; i<times; i++) 
  for (j=0; j<times; j++) 
  for (k=0; k<times; k++) {
    t1 = ti + i;
    t2 = tj + j;
    t3 = tk + k;
    if (data_array[t1][t2][t3] == 1) return 1; 
  }
  return 0;
}

int
refine_border (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  int                 qid;
  int                 tilelen;
  int                 i, j, k;
  int                 markedCount;
  double                 offsi, offsj, offsk;
  double                 tx, ty, tz, temp;
  int ix, iy, iz;
  p4est_qcoord_t      unit_len, current_len;
  int cur_level;
  // the border does not have to be a high refine_level
  cur_level = sk_level;

  if ((int) q->level > cur_level) {
    return 0;
  }
  if ((int) q->level < min_level) {
    return 1;
  }
  tilelen = 1 << (cur_level - q->level);       /* Pixel size of quadrant */
  unit_len = P4EST_QUADRANT_LEN (cur_level);

  offsi = (double)(q->x) / unit_len + 0.5;       /* Pixel x offset */
  offsj = (double)(q->y) / unit_len + 0.5;       /* Pixel y offset */
  offsk = (double)(q->z) / unit_len + 0.5;       /* Pixel z offset */
  
  P4EST_ASSERT (offsi >= 0 && offsj >= 0 && offsk >= 0);
  markedCount = 0; 
  for (k = 0; k < tilelen; ++k) {
    for (j = 0; j < tilelen; ++j) {
      for (i = 0; i < tilelen; ++i) {
        tx = offsi + i;
        ty = offsj + j;
        tz = offsk + k; 
        ix = (int)(offsi-0.5) + i;
        iy = (int)(offsj-0.5) + j;
        iz = (int)(offsk-0.5) + k;
        //if (test_array2[ix][iy][iz] == 1) markedCount += 1;
        if (check_data(test_array2, ix, iy, iz, cur_level) == 1) markedCount += 1;
      }
    }
  }
  if (!markedCount) {
          q->p.user_int = UNMARKED;
          return 0;
  } else if (tilelen*tilelen*tilelen == 1 && markedCount == 1) {
          q->p.user_int = MARKED;
          return 0;
  } else {
          return 1;
  }
  return 0;
}

int*** initialize_array(int len){
  int*** data_array;
  int i, j, k;
  data_array = (int***) malloc(len * sizeof(int **));
  for (i=0; i<len; i++){
    data_array[i] = (int **) malloc ( len* sizeof(int*) );
    for (j=0; j<len; j++)
      data_array[i][j] = (int *) malloc ( len* sizeof(int) );
  }
  for ( i=0; i<len; i++ )
  for ( j=0; j<len; j++ )
  for ( k=0; k<len; k++ )
    data_array[i][j][k] = 0;
  return data_array;
}

int*** read_array(const char *file_name){
  int i, j, k;
  int len;
  FILE * fp;
  char *line = NULL;
  ssize_t read;
  size_t slen = 0;
  int data;
  int*** test_array;

  fp = fopen(file_name, "r");
  //fp = fopen("/home/xin/Dropbox/3d-printing-paper/vs-projects/Project1/Project1/stl-files/head-volume", "r");
  //fp = fopen("/home/xin/Dropbox/3d-printing-paper/vs-projects/Project1/Project1/stl-files/candle-volume", "r");
  //fp = fopen("/home/xin/Dropbox/3d-printing-paper/vs-projects/Project1/Project1/stl-files/teapot-volume", "r");
  if (fp == NULL){
    printf(" can not open file ! \n");
    exit(1);
  }
  read = getline(&line, &slen, fp);
  len = (int)read - 2;
  test_array = (int***) malloc(len * sizeof(int **));
  if (cube_len == -1) {
    cube_len = len;
    printf( " >>>>> reading bit arrays: %d !!!\n ", cube_len );
  }

  for (i=0; i<len; i++){
    test_array[i] = (int **) malloc ( len* sizeof(int*) );
    for (j=0; j<len; j++)
      test_array[i][j] = (int *) malloc ( len* sizeof(int) );
  }
  //printf ("length %d line: %s \n", len, line);
  fseek(fp, 0, SEEK_SET);

  for ( i=0; i<len; i++ )
  for ( j=0; j<len; j++ ) {
    if ((read = getline(&line, &slen, fp)) != -1){
      //printf ("length %zu line: %s \n", read, line);
      for ( k=0; k<len; k++ ){
        test_array[i][j][k] = (int)line[k] - 48;
        // 
        //if (test_array[i][j][k] == 0) test_array[i][j][k] = 2;
      }
    } 
  }
  fclose(fp);
  printf ("finish loading the volume file into array!\n");
  return test_array;
}

//has segmentation faults
void fill_neighbors(int*** data_array, int i, int j, int k, int len){
  int s1, s2, s3;
  int ip, jp, kp;
  int temp_array[29];
  int count = 0;

  for (s1=0; s1<3; s1++)
  for (s2=0; s2<3; s2++)
  for (s3=0; s3<3; s3++){
    ip = i-1+s1;
    jp = j-1+s2;
    kp = k-1+s3; 
    if (ip>=0 && ip<len && jp>=0 && jp<len && kp>=0 && kp<len){
      printf(" xxxxxx %d, %d, %d\n ", ip, jp, kp);
      if (data_array[ip][jp][kp] == 2) {
        data_array[ip][jp][kp] = 0; 
        //fill_neighbors(data_array, ip, jp, kp, len);
        temp_array[count] = ip;
        temp_array[count+1] = jp;
        temp_array[count+2] = kp;
        count += 3;
      } 
      //data_array[ip][jp][kp] = 
    } 
  }
  for (s1=0; s1<count/3; s1++){
    fill_neighbors(data_array, temp_array[3*s1], temp_array[3*s1+1], temp_array[3*s1+2], len);
  }
}

//void fill_out_cubes(int*** data_array, int len){
//  fill_neighbors(data_array, 0, 0, 0, len);
//}

void fill_one_line(int* data_array, int len){
  int lid, rid;
  int next_solid;
  lid = 0;
  rid = len-1;
  next_solid = 0;
  while (lid < len && rid >= 0 && lid < rid){
    if (next_solid == 1){
      while(lid < len && data_array[lid] == 0 ){ data_array[lid] = 2; lid++; }
      while(rid >= 0  && data_array[rid] == 0 ){ data_array[rid] = 2; rid--; }
      next_solid = 0;
    } else{
      while( lid <len && data_array[lid] == 0 ){ lid++; }
      while( rid >= 0 && data_array[rid] == 0 ){ rid--; }
      next_solid = 1;
    }
    if (lid >= len || rid < 0 || lid >= rid)break;
    if (data_array[lid] == data_array[rid] && data_array[lid] == 1){
      lid++; rid--;
    } else{
      printf("error here: left %d: %d, right %d: %d \n", 
	lid, data_array[lid], rid, data_array[rid] );
      return;
    }
  }
}

// assume the object is continuous
void fill_out_cube(int*** data_array, int len){
  int i, j, k;
  int *data_temp;
  // the start is always an element of 0
  // this indicates that the inner part has been signed
  int data;
  for (i=0; i<len; i++){
    for (j=0; j<len; j++){
      data_temp = data_array[i][j];
      fill_one_line(data_temp, len); 
    }
  }
}

void grow_element(int*** data_array, int i, int j, int k, int dis, int len){
  int min_i, min_j, min_k;
  int max_i, max_j, max_k;
  double s1, s2, s3, d1, d2, d3, r;
  int si, sj, sk;
  if (dis >= len){
    printf("the growing distance is much larger than len!!!\n");
    return;
  }
  if ( i<0 || j<0 || k<0 || i>=len || j>=len || k>=len ){
    printf("index is wrong!!! %d, %d, %d \n", i, j, k);
    return;
  }
  if ( data_array[i][j][k] != 1 ){
    printf("the element that needs to grow does not contain data !!!\n");
    return;
  }
  d1 = i+0.5;
  d2 = j+0.5;
  d3 = k+0.5;
  r = dis;
  min_i = i - dis;
  min_j = j - dis;
  min_k = k - dis;
  max_i = i + dis;
  max_j = j + dis;
  max_k = k + dis;
  if (min_i < 0) min_i = 0;
  if (min_j < 0) min_j = 0;
  if (min_k < 0) min_k = 0;
  if (max_i > len-1) max_i = len-1; 
  if (max_j > len-1) max_j = len-1; 
  if (max_k > len-1) max_k = len-1; 
  for (si=min_i; si<=max_i; si++)
  for (sj=min_j; sj<=max_j; sj++)
  for (sk=min_k; sk<=max_k; sk++){
    s1 = (double)si+0.5; 
    s2 = (double)sj+0.5; 
    s3 = (double)sk+0.5; 
    if ( (s1-d1)*(s1-d1) + (s2-d2)*(s2-d2) + (s3-d3)*(s3-d3) <= r*r && 
          data_array[si][sj][sk] == 0){
      data_array[si][sj][sk] = 3;
    } 
  }
}

// current version, the distance grows on double sides
void volume_offset(int*** data_array, int dis, int len){
  int i, j, k;
  for (i=0; i<len; i++)
  for (j=0; j<len; j++)
  for (k=0; k<len; k++){
    if (data_array[i][j][k] == 1){
      grow_element(data_array, i, j, k, dis, len); 
    }
  } 
}

void path_offset(int*** data_array, int p, int dis, int len){
  int i, j, k;
  int start_sd = p;
  for (i=0; i<len; i++)
  for (j=0; j<len; j++)
  for (k=start_sd; k<start_sd+1; k++){
    if (data_array[i][j][k] == 1){
      grow_element(data_array, i, j, k, dis, len); 
    }
  } 
}

void set_outer(int*** data_array, int dis, int len){
  int i, j, k;
  for (i=0; i<len; i++)
  for (j=0; j<len; j++)
  for (k=0; k<len; k++){
    if ( i<=dis || j<=dis || k<=dis 
      || i>=len-dis || j>=len-dis || k>=len-dis ){
      data_array[i][j][k] = 1;
    } 
  }
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

  p4est_t            *p4est3;
  p4est_t            *p4est_cut;
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
  const char *fname1 = "/home/xin/Dropbox/3d-printing-paper/vs-projects/Project1/Project1/stl-files/head-volume";
  //const char *fname1 = "/home/xin/Dropbox/3d-printing-paper/vs-projects/Project1/Project1/stl-files/teapot-volume";
  //const char *fname1 = "/home/xin/Dropbox/3d-printing-paper/vs-projects/Project1/Project1/stl-files/candle-volume";
  test_array1 = read_array(fname1);  
  fill_out_cube(test_array1, cube_len);
  //test_array2 = read_array(fname2);  

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

  p4est1 = p4est_new_ext (mpi->mpicomm, connectivity, 0, min_level, 1, 0, NULL, NULL);
  p4est2 = p4est_copy(p4est1, 1); 
  p4est_out = p4est_copy(p4est1, 1); 
  p4est3 = p4est_copy(p4est1, 1); 
  p4est_cut = p4est_copy(p4est1, 1); 

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

  printf("start now ... offset\n");

  p4est_refine (p4est1, 1, refine_sphere1, NULL);
  test_array2 = initialize_array(cube_len);
  set_outer(test_array2, 32, cube_len); 
  p4est_refine (p4est2, 1, refine_border, NULL);
  p4est_diff(p4est2, p4est1, p4est_out, p4est_diff_aafn);

gettimeofday(&t1, NULL);
  sort_p4est(p4est_out);
  //volume_offset(test_array1, 5, cube_len); 
  //for ( i=5; i < cube_len; ){
  //path_offset(test_array1, i, 5, cube_len); 
  //i += 10;
  //} 
  //path_offset(test_array1, 18, 10, cube_len); 
  //path_offset(test_array1, 19, 10, cube_len); 
gettimeofday(&t2, NULL);
  elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
  printf(" >>>>>>>  operation p4est offset elapsedTime %f  ms.\n", elapsedTime);

//printf( " start cutting !!!!!!!!!!!!!!!!!!!\n" );
  p4est_refine (p4est3, 1, refine_full, NULL);
  p4est_diff(p4est3, p4est_out, p4est_cut, p4est_diff_aafn);

  /*
  p4est_refine (p4est1, 1, refine_test, NULL);
  p4est_refine (p4est2, 1, refine_full, NULL);
  p4est_diff(p4est2, p4est1, p4est_out, p4est_diff_aafn);
  //check_data_fields(p4est_out);
  p4est_refine (p4est3, 1, refine_full, NULL);
  print_sign = 0;
  printf( " start cutting !!!!!!!!!!!!!!!!!!!\n" );
  test_break();
  p4est_diff(p4est3, p4est_out, p4est_cut, p4est_diff_aafn);
  */

  //check_data_fields(p4est3);
  
  //p4est_diff(p4est3, p4est_out, p4est_out1, p4est_diff_aafn);
  //p4est_refine (p4est1, 1, refine_offset, NULL);
  //p4est_intersection(p4est1, p4est2, p4est_out);
  //p4est_union(p4est1, p4est2, p4est_out);
  //p4est_diff(p4est2, p4est1, p4est_out, p4est_diff_aafn);

  printf( " ====> removing all empty quadrants!\n" );
  p4est_remove(p4est1);
  //p4est_remove(p4est2);
  printf( " ====> removing all empty for border!\n" );
  p4est_remove(p4est_out);
  p4est_remove(p4est3);
  p4est_remove(p4est_cut);

  //p4est_remove(p4est_out1);
  //sc_stats_set1 (&stats[TIMINGS_REFINE], snapshot.iwtime, "Refine");
  printf("finish operations here, start writing file !!!\n");

#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est1, NULL, "tree_head");
  //p4est_vtk_write_file (p4est1, NULL, "tree1_teapot");
  p4est_vtk_write_file (p4est2, NULL, "tree_second");
  p4est_vtk_write_file (p4est3, NULL, "tree_full");
  p4est_vtk_write_file (p4est_out, NULL, "tree_border");
  p4est_vtk_write_file (p4est_cut, NULL, "tree_cut");
#endif
  count_refined = p4est1->global_num_quadrants;

  /* destroy the p4est and its connectivity structure */
  P4EST_FREE (quadrant_counts);
  P4EST_FREE (p4est1->inspect);
  p4est_destroy (p4est1);
  p4est_destroy (p4est2);
  p4est_destroy (p4est_out);

  p4est_destroy (p4est3);
  p4est_destroy (p4est_cut);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
