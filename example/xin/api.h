#ifndef API_H
#define API_H
#endif

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


#define PI 3.14159265
#define P4EST_TIMINGS_VTK 

#define P4EST_REFINE_LEVEL 8 
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
static const p4est_qcoord_t offset = P4EST_REFINE_LENGTH/64;

/*
int
coarsen_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * q[]);

int in_the_mid(p4est_qcoord_t v1, p4est_qcoord_t v2, p4est_qcoord_t x);


static int 
check_overlap(p4est_qcoord_t ox, p4est_qcoord_t oy,
#ifdef P4_TO_P8
p4est_qcoord_t oz, 
#endif
p4est_qcoord_t len, p4est_qcoord_t cx, p4est_qcoord_t cy,
#ifdef P4_TO_P8
p4est_qcoord_t cz, 
#endif
p4est_qcoord_t radix);


int
refine_sphere_test (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q);

int p4est_diff_aafn (p4est_t *p4est_in1, p4est_t *p4est_in2, p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant);

int
refine_fn1 (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q);

int
refine_fn2 (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q);

int
refine_bowl (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q);

int
refine_sphere1 (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q);

int
refine_sphere2 (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q);

int
refine_sphere3 (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q);

*/



