if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += \
  example/set_operations/p4est_setop1 \
  example/set_operations/p4est_setop2
  
example_set_operations_p4est_setop1_SOURCES = example/set_operations/p4est_setop1.c
example_set_operations_p4est_setop2_SOURCES = example/set_operations/p4est_setop2.c

LINT_CSOURCES += $(example_set_operations_p4est_setop1_SOURCES)
LINT_CSOURCES += $(example_set_operations_p4est_setop2_SOURCES)
endif

EXTRA_DIST += \
  example/steps/hw32.h \
  example/set_operations/rect_horizontal1.h \
  example/set_operations/rect_vertical1.h \
  example/set_operations/rect_union.h
