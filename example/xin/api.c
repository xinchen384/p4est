#include "api.h"


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

int
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

int
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

int
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

int
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






