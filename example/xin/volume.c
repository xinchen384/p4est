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
#include <math.h>

typedef struct
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

typedef struct
{
  double         x;
  double         y;
  double         z;
}
Point3d;


int*** test_array1;
int*** test_array2;
int cube_len = -1;
int sk_level = 4;

int print_sign = 0;

// tool configurations
double tool_len[5];
double tool_r[5];
int tool_pieces = 5;
double tool_length = 0;
double val = 180.0/PI;
int map_row = 64;
int map_col = 128;
//int map_row = 4;
//int map_col = 8;
int** map;
Point3d **rotate_vector;


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
        printf("\n quadrant counts: %d, the element coutn in the list: %d\n", myc, (int)(array->elem_count));
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


int check_tool_border(int id, double dist){
  double min, max, minr, maxr;
  int pid;
  double radius; 
  pid = tool_pieces - 1; 
  if ( id > pid || id < 0 ){
    printf(" *** the tool id is out of range !!! error: %d\n", id );
    return 0;
  }
  min = 0;
  while (pid > id){
    //printf(" the steps: %f, ", tool_len[pid]);
    min += tool_len[pid]; 
    pid--; 
  }
  max = min + tool_len[id]; 
  radius = tool_r[id];
  minr = min*min + radius*radius;
  maxr = max*max + radius*radius;
  if (dist*dist >= maxr){
    if (pid <= 0){
      printf(" *** calculation is wrong, does not belong this piece id %d, %f %f > %f !!!\n", id, max, radius, dist);
      return 0;
    }else{
      // on the upper level
      return -1;
    }
  }
  if (dist*dist < minr){
    // go to the lower level 
    return 1;
  }
  if (dist*dist < maxr && dist*dist >= minr){
    return 0;
  }
  printf(" *** unacceptable condition in checking tool border !!!\n");
  return 0;
}

double get_tool_intersect(double dist){
  // assume the dist starts from the central point of the tool
  double s = dist;
  double sum, angle ;
  int pid = tool_pieces - 1;
  double pieceR = -1.0;
  int i, ret, k;

  if (dist >= tool_length){
    //can not reach
    return 0;
  }

  while (s > 0 && pid >= 0){
    s = s - tool_len[pid];
    pid--;
  }
  
  //if (dist > 224) printf("left %f pid %d\n", s, pid);
  //printf("end checking the length, ");
  if (s < 0)pid++;

  if (pid < 0){
    printf(" the tool is not long enough %f to reach the dist %f !!!\n", tool_length, dist);
    //indicate that  it is not accessible at any cases
    return 360.0;
  } else{
    //printf("start checking the arc, id %d, dist %f,  ", pid, dist);
    ret = check_tool_border(pid, dist);
    while (pid<= (tool_pieces-1) && ret == 1) {
      ret = check_tool_border(pid, dist);
      pid++;
    } 
    
    //printf("end checking the arc, id %d, dist %f,  ", pid, dist);
    if (pid > tool_pieces - 1){
      printf(" there must be something wrong with the tool border calculation!!! \n ");
      return -1.0;
    }
    if (ret == 0){
      pieceR = tool_r[pid]; 
      angle = asin(pieceR/dist) * val; 
      if (isnan(angle)){ 
        printf(" short radius %f at id: %d with distance %f \n", pieceR, pid, dist); 
      }
      return angle;
    }else if (ret == -1){
      sum = 0;
      k = tool_pieces - 1;
      while (k >= pid){
        sum += tool_len[k]; 
        k--; 
      }    
      angle = acos(sum/dist) * val; 
      if (isnan(angle)){ 
        printf(" short radius %f at id: %d with distance %f \n", pieceR, pid, dist); 
      }
      return angle;
    }
  }
  printf(" calculating the tool border angle cause error !!!\n");
  return -1.0; 
}

void initialize_rotate_vector(){
	double sita, sigma;
	int i, j;
	rotate_vector = (Point3d **)malloc(map_row * sizeof(Point3d *));
        for (i=0; i<map_row; i++)
	  rotate_vector[i] = (Point3d*) malloc(map_col * sizeof(Point3d));
	for (i=0; i<map_row; i++)
	for (j=0; j<map_col; j++){
	    sita = j * 180.0/map_col / val;
	    sigma = i * 360.0/map_row / val; 
	    rotate_vector[i][j].x = sin(sita)*cos(sigma);
	    rotate_vector[i][j].y = sin(sita)*sin(sigma);
	    rotate_vector[i][j].z = cos(sita);
	}
	//for (i=0; i<map_row; i++){
	//  for (j=0; j<map_col; j++)
	//    printf("(%f, %f, %f), ", rotate_vector[i][j].x, rotate_vector[i][j].y, rotate_vector[i][j].z); 
	//  printf("\n");
	//}
}

void mark_angle_range(double sita1, double sita2, double sigma1, double sigma2, 
		double offsi, double offsj, double offsk, double dist, double ica){
	double sita, sigma;
	double vec_x, vec_y, vec_z;
	double vk, temp;
	double sita_unit = 180.0/(double)map_col;
	double sigma_unit = 360.0/(double)map_row;
	int rid, cid;
	int cid1, cid2;
	int rid1, rid2;
	int rt, ct, r, c;

	cid = (int)(sita1/sita_unit);
	if (sita1 <= ((double)cid)*sita_unit) cid1 = cid;
	else cid1 = cid + 1;
	cid2 = (int)(sita2/sita_unit);

	rid = (int)(sigma1/sigma_unit);
	if (sigma1 <= ((double)rid)*sigma_unit) rid1 = rid;
	else rid1 = rid + 1;
	rid2 = (int)(sigma2/sigma_unit);

	vk = cos(ica/val);
	for ( c=cid1; c<=cid2; c++ )
	for ( r=rid1; r<=rid2; r++ ){ 
          rt = r;
	  ct = c;
	  //if (rt < 0) rt = rt+map_row;
	  //if (rt >= map_row) rt = rt-map_row;
	  //note: sigma may be out of range
	  //if (rt < 0) rt = rt+map_row;
	  //if (rt >= map_row) rt = rt-map_row;
	  //sita = ct * 180.0/map_col / val;
	  //sigma = rt * 360.0/map_row / val; 
	  //vec_x = sin(sita)*cos(sigma);
	  //vec_y = sin(sita)*sin(sigma);
	  //vec_z = cos(sita);
	  vec_x = rotate_vector[rt][ct].x;
	  vec_y = rotate_vector[rt][ct].y;
	  vec_z = rotate_vector[rt][ct].z;
	  temp = (offsi*vec_x + offsj*vec_y + offsk*vec_z)/dist;
	  //angle = acos(temp) * val;
	  if ( map[rt][ct]==0 && temp>=vk ){
	    map[rt][ct] = 1;
	    //printf(" ======>>>>> check %f angle %f, ==== %f ica is %f\n", temp, angle, vk, ica);
	  }
	} 
	//printf("end \n");
}

void mark_map(p4est_t *p4est, int x, int y, int z){
        sc_array_t *quads;
	p4est_tree_t *tree;
        sc_array_t *array;
        p4est_quadrant_t *q;
	p4est_qcoord_t      unit_len;
  	struct timeval t1, t2;
  	double elapsedTime;
        int i, j, k;
        int elem_size;
        int myc;
	int tilelen;
	int ica_id;

  	double offsi, offsj, offsk, radius;
	double rad, dist;
        double ica1, ica2, ica;

	double sita, sigma;
	int cid, cid1, cid2;
	int rid, rid1, rid2;
	int c, r, ct, rt;
	double x1, x2, y1, y2;
	double sita_unit = 180.0/(double)map_col;
	double sigma_unit = 360.0/(double)map_row;

	double vec_x, vec_y, vec_z;
        double temp, angle, ica_angle;
	double start_angle, end_angle, angle1;
	double *ICA_list;
	double radS, v_sig, vk;
	char *s;
	
        unit_len = P4EST_QUADRANT_LEN (refine_level);
        tree = p4est_tree_array_index(p4est->trees, 0);
	quads = &(tree->quadrants);
	array = quads;

gettimeofday(&t1, NULL);
	myc = 0;
        s = array->array;
	for (i=0; i<array->elem_count; i++){
          q = (p4est_quadrant_t *) s;
	  if (q->p.user_int == 1){
            myc++;
	  }
          //printf("%d ", q->p.user_int); 
          s += array->elem_size;
        }
        printf("\n total quadrant counts: %zu, data fields count: %d\n", array->elem_count, myc);
gettimeofday(&t2, NULL);
	elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
	elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
	printf(" >>>>>>>  one iteration takes %f  ms.\n", elapsedTime);

	// initialize the accessibility map
	// initialize the map as accessible (0)
	map = (int **)malloc(map_row * sizeof(int *));
        for (i=0; i<map_row; i++)
	  map[i] = (int*) malloc(map_col * sizeof(int));
	for (i=0; i<map_row; i++)
	for (j=0; j<map_col; j++){
	  map[i][j] = 0;
	}
	
	// initialize ICA list
        elem_size = (int) array->elem_count;
        ICA_list = malloc( myc*sizeof(double) );
	for (i=0; i<myc; i++)
	  ICA_list[i] = 0.0;
	// generate ICA for each element
	s = array->array;
	ica_id = 0;
        for (i=0; i<array->elem_count; i++){
          q = (p4est_quadrant_t *) s;
	  if (q->p.user_int == 1){
   	  offsi = (double)(q->x) / unit_len + 0.5 - x;       /* Pixel x offset */
	  offsj = (double)(q->y) / unit_len + 0.5 - y;       /* Pixel y offset */
	  offsk = (double)(q->z) / unit_len + 0.5 - z;       /* Pixel z offset */
          dist = sqrt(offsi*offsi + offsj*offsj + offsk*offsk);
          //printf(" start calculating ICA 1 with dist: %f \n ", dist);
          ica1 = get_tool_intersect(dist); 
  	  tilelen = 1 << (refine_level - q->level);       
          rad = (double)tilelen/2;
          ica2 = asin(rad/dist) * val;
	  ica = ica1 + ica2;
          if (isnan(ica)){
     	    printf(" detecting nan number for ica: 1: %f, 2: %f \n", ica1, ica2);
	  }
          ICA_list[ica_id] = ica; 
	  ica_id++;
	  if (ica >= 90) printf("note: large ica %f!!!!!\n", ica);
	  
	  for ( r=0; r<map_row; r++ )
	  for ( c=0; c<map_col; c++ ){ 
	    //sita = c * 180.0/map_col / val;
	    //sigma = r * 360.0/map_row / val; 
	    //vec_x = sin(sita)*cos(sigma);
	    //vec_y = sin(sita)*sin(sigma);
	    //vec_z = cos(sita);
	    vec_x = rotate_vector[r][c].x;
	    vec_y = rotate_vector[r][c].y;
	    vec_z = rotate_vector[r][c].z;
	    temp = (offsi*vec_x + offsj*vec_y + offsk*vec_z)/dist;
	    angle = acos(temp) * val;
	    if ( angle <= ica ) map[r][c] = 1;
	  } 

	  /*
	  sita = acos( offsk/dist ) * val;
	  sigma = atan( offsj/offsi ) * val;

	  if (sita - ica < 0){
	    mark_angle_range(0, sita+ica, 0, 360-sigma_unit, offsi, offsj, offsk, dist, ica);
	  } 
	  else if (sita + ica > 180){
	    mark_angle_range(sita-ica, 180-sita_unit, 0, 360-sigma_unit, offsi, offsj, offsk, dist, ica);
	  } else {
	    v_sig = asin( sin(ica/val)/sin(sita/val) ) * val;
	    if (sigma-v_sig < 0 || sigma+v_sig >= 360)
	      mark_angle_range(sita-ica, sita+ica, 0, 360-sigma_unit, offsi, offsj, offsk, dist, ica);
	    else
	      mark_angle_range(sita-ica, sita+ica, sigma-v_sig, sigma+v_sig, offsi, offsj, offsk, dist, ica);
	  }
	  */

	  /*
	  if (start_angle < 0 || end_angle > 180){
	  for ( r=0; r<map_row; r++ )
	  for ( c=0; c<map_col; c++ ){ 
	    //sita = c * 180.0/map_col / val;
	    //sigma = r * 360.0/map_row / val; 
	    //vec_x = sin(sita)*cos(sigma);
	    //vec_y = sin(sita)*sin(sigma);
	    //vec_z = cos(sita);
	    vec_x = rotate_vector[r][c].x;
	    vec_y = rotate_vector[r][c].y;
	    vec_z = rotate_vector[r][c].z;
	    temp = (offsi*vec_x + offsj*vec_y + offsk*vec_z)/dist;
	    angle = acos(temp) * val;
	    if ( angle <= ica ) map[r][c] = 1;
	  }
	  } else {
	    v_sig = asin( sin(ica/val)/sin(sita/val) ) * val;
	    //printf(" two angles cal: %f, ica %f \n ", v_sig, ica );
	    //v_sig = asin( (ica/360.0)/sin(sita/val) ) * val;
	    if (sigma-v_sig < 0 || sigma+v_sig >= 360)
	      mark_angle_range(start_angle, end_angle, 0, 360-sigma_unit, offsi, offsj, offsk, dist, ica);
	    else
	      mark_angle_range(start_angle, end_angle, sigma-v_sig, sigma+v_sig, offsi, offsj, offsk, dist, ica);
	  }
	  */
	  
          //if (ica_id % 10000 == 0) printf("progress %d \n", ica_id); 
	  }
          s += array->elem_size;
        }
        printf(" finish generating ICA list ! %d %d\n ", ica_id, myc);
	for (i=0; i<myc; i++){
          //if (isnan(ICA_list[i]))
          if (ICA_list[i] >= 360.0)
	  printf("360 element: %f ", ICA_list[i]);
	}
	//printf("\nend ica list two counts: %d %d\n", ica_id, myc);
        printf(" finished calculating the access map \n");
}

// note that the coordinates should be normalized
void mark_accessibility_map(p4est_t *p4est, int x, int y, int z){
        sc_array_t *quads;
	p4est_tree_t *tree;
        sc_array_t *array;
        p4est_quadrant_t *q;
	p4est_qcoord_t      unit_len;
  	struct timeval t1, t2;
  	double elapsedTime;
        int i, j, k;
        int elem_size;
        int myc;
	int tilelen;
	int ica_id;

  	double offsi, offsj, offsk, radius;
	double rad, dist;
        double ica1, ica2;

	double sita, sigma;
	double vec_x, vec_y, vec_z;
        double temp, angle, ica_angle;
	double *ICA_list;
	char *s;

        unit_len = P4EST_QUADRANT_LEN (refine_level);
        tree = p4est_tree_array_index(p4est->trees, 0);
	quads = &(tree->quadrants);
	array = quads;

gettimeofday(&t1, NULL);
	myc = 0;
        s = array->array;
	for (i=0; i<array->elem_count; i++){
          q = (p4est_quadrant_t *) s;
	  if (q->p.user_int == 1){
            myc++;
	  }
          //printf("%d ", q->p.user_int); 
          s += array->elem_size;
        }
        printf("\n total quadrant counts: %zu, data fields count: %d\n", array->elem_count, myc);
gettimeofday(&t2, NULL);
	elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
	elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
	printf(" >>>>>>>  one iteration takes %f  ms.\n", elapsedTime);

        elem_size = (int) array->elem_count;
        ICA_list = malloc( myc*sizeof(double) );
	for (i=0; i<myc; i++)
	  ICA_list[i] = 0.0;

	// generate ICA for each element
	s = array->array;
	ica_id = 0;
        for (i=0; i<array->elem_count; i++){
          q = (p4est_quadrant_t *) s;
	  if (q->p.user_int == 1){
   	  offsi = (double)(q->x) / unit_len + 0.5;       /* Pixel x offset */
	  offsj = (double)(q->y) / unit_len + 0.5;       /* Pixel y offset */
	  offsk = (double)(q->z) / unit_len + 0.5;       /* Pixel z offset */
          dist = (offsi-x)*(offsi-x) + (offsj-y)*(offsj-y) + (offsk-z)*(offsk-z);
          dist = sqrt(dist);
          //printf(" start calculating ICA 1 with dist: %f \n ", dist);
          ica1 = get_tool_intersect(dist); 
  	  tilelen = 1 << (refine_level - q->level);       
          rad = (double)tilelen/2;
          ica2 = asin(rad/dist) * val;
          //if (isnan(ica)){
     	  //  printf(" detecting nan number: 1: %f, 2: %f \n", ica1, ica2);
	  //}
          ICA_list[ica_id] =ica1 + ica2; 
	  ica_id++;
	  }
          //printf("%d ", q->p.user_int); 
          s += array->elem_size;
        }
        printf(" finish generating ICA list ! %d %d\n ", ica_id, myc);
	for (i=0; i<myc; i++){
          //if (isnan(ICA_list[i]))
          if (ICA_list[i] >= 360.0)
	  printf("360 element: %f ", ICA_list[i]);
	}
	//printf("\nend ica list two counts: %d %d\n", ica_id, myc);

	map = (int **)malloc(map_row * sizeof(int *));
        for (i=0; i<map_row; i++)
	  map[i] = (int*) malloc(map_col * sizeof(int));
	// initialize the map as accessible (0)
	for (i=0; i<map_row; i++)
	for (j=0; j<map_col; j++){
	  map[i][j] = 0;
	}

        for ( i=0; i<map_row; i++ ){
        for ( j=0; j<map_col; j++ ){
	  sita = j * 180.0/map_col / val;
	  sigma = i * 360.0/map_row / val; 
	  vec_x = sin(sita)*cos(sigma);
	  vec_y = sin(sita)*sin(sigma);
	  vec_z = cos(sita);

          s = array->array;
	  ica_id = 0;
          for (k=0; k<array->elem_count; k++){
            q = (p4est_quadrant_t *) s;
	    if (q->p.user_int == 1){
   	    offsi = (double)(q->x) / unit_len + 0.5;       /* Pixel x offset */
	    offsj = (double)(q->y) / unit_len + 0.5;       /* Pixel y offset */
	    offsk = (double)(q->z) / unit_len + 0.5;       /* Pixel z offset */
            dist = (offsi-x)*(offsi-x) + (offsj-y)*(offsj-y) + (offsk-z)*(offsk-z);
            dist = sqrt(dist);
	    temp = ((offsi-x)*vec_x + (offsj-y)*vec_y + (offsk-z)*vec_z)/dist;
	    angle = acos(temp) * val;
            ica_angle = ICA_list[ica_id];
            ica_id++;
	    //if (ICA_list[ica_id] < 360.0){
            //  printf ("the angle is %f, ica angle is %f\n", angle, ICA_list[ica_id]);
 	    //}
            if (angle <= ica_angle){
              //printf ("check,");
	      map[i][j] = 1;
	      break;
            }  
	    }
            //printf("%d ", q->p.user_int); 
            s += array->elem_size;
          }
	}
          printf(" progress i: %d, j: %d \n", i, j );
	}
        printf(" finished calculating the access map \n");
}

void print_map(){
  int i, j;
  for (i=0; i<map_row; i++){
    for (j=0; j<map_col; j++){
      printf("%d", map[i][j]);
      //if (map[i][j] == 1) printf("a");
    }
    printf("\n");
  }
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

  double x, ret;
tool_len[0] = 22; // start 0
tool_len[1] = 78;  // start 22
tool_len[2] = 15.94; // start at 100
tool_len[3] = 5.08;  // start at 115.94
tool_len[4] = 2.28;  // start at 121.02

tool_r[0] = 31.5;
tool_r[1] = 25;
tool_r[2] = 1.59;
tool_r[3] = 1.59;
tool_r[4] = 0.79375;

  for (i=0; i<tool_pieces; i++){
    tool_len[i] *= 2;
    tool_r[i] *= 2;
    tool_length += tool_len[i];
  }
  //*********
  struct timeval t1, t2, t3, t4, t5;
  double elapsedTime;
  //*********
  const char *fname1 = "/home/xin/Dropbox/3d-printing-paper/vs-projects/Project1/Project1/stl-files/head-volume";
  //const char *fname1 = "/home/xin/Dropbox/3d-printing-paper/vs-projects/Project1/Project1/stl-files/teapot-volume";
  //const char *fname1 = "/home/xin/Dropbox/3d-printing-paper/vs-projects/Project1/Project1/stl-files/candle-volume";
  test_array1 = read_array(fname1);  
  //fill_out_cube(test_array1, cube_len);
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

  x = 0.5;
  ret = asin(x) * val;
  printf (" test asin in C, sin(%f) = %f \n", ret, x );
  //return 1;
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

  initialize_rotate_vector();
  p4est_refine (p4est1, 1, refine_sphere1, NULL);
  // ======>  used to test cut layer by layer 
  //test_array2 = initialize_array(cube_len);
  //set_outer(test_array2, 32, cube_len); 
  //p4est_refine (p4est2, 1, refine_border, NULL);
  //p4est_diff(p4est2, p4est1, p4est_out, p4est_diff_aafn);
  // this step becomes mandatory as the diff operation does not guarantee the sorting
  //sort_p4est(p4est_out);
  //printf( " start cutting !!!!!!!!!!!!!!!!!!!\n" );
  //p4est_refine (p4est3, 1, refine_full, NULL);
  //p4est_diff(p4est3, p4est_out, p4est_cut, p4est_diff_aafn);

  //  =========> used to test grow 
  //volume_offset(test_array1, 5, cube_len); 
  //for ( i=5; i < cube_len; ){
  //path_offset(test_array1, i, 5, cube_len); 
  //i += 10;
  //} 
  //path_offset(test_array1, 18, 10, cube_len); 
  //path_offset(test_array1, 19, 10, cube_len); 

gettimeofday(&t1, NULL);
  //mark_accessibility_map(p4est1, 128, 128, 255);
  mark_map(p4est1, 128, 128, 255);
gettimeofday(&t2, NULL);
  elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
  printf(" ===========  >>>>>>>>>>>  operation p4est accessMap elapsedTime %f  ms.\n", elapsedTime);

  print_map();
  //check_data_fields(p4est_out);
  //check_data_fields(p4est1);
 
  printf( " ====> removing all empty quadrants!\n" );
  p4est_remove(p4est1);
  //p4est_remove(p4est2);
  //printf( " ====> removing all empty for border!\n" );
  //p4est_remove(p4est_out);
  //p4est_remove(p4est3);
  //p4est_remove(p4est_cut);

  //p4est_remove(p4est_out1);
  //sc_stats_set1 (&stats[TIMINGS_REFINE], snapshot.iwtime, "Refine");
  printf("finish operations here, start writing file !!!\n");

#ifdef P4EST_TIMINGS_VTK
  p4est_vtk_write_file (p4est1, NULL, "tree_head");
  //p4est_vtk_write_file (p4est1, NULL, "tree1_teapot");
  //p4est_vtk_write_file (p4est2, NULL, "tree_second");
  //p4est_vtk_write_file (p4est3, NULL, "tree_full");
  //p4est_vtk_write_file (p4est_out, NULL, "tree_border");
  //p4est_vtk_write_file (p4est_cut, NULL, "tree_cut");
#endif
  count_refined = p4est1->global_num_quadrants;

  /* destroy the p4est and its connectivity structure */
  P4EST_FREE (quadrant_counts);
  P4EST_FREE (p4est1->inspect);
  p4est_destroy (p4est1);
  p4est_destroy (p4est2);
  p4est_destroy (p4est3);
  p4est_destroy (p4est_out);
  p4est_destroy (p4est_cut);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
