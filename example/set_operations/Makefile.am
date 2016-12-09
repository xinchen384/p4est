if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += \
  example/set_operations/p4est_setop2 \
  example/set_operations/p4est_setop3
  
example_set_operations_p4est_setop2_SOURCES = example/set_operations/p4est_setop2.c
example_set_operations_p4est_setop3_SOURCES = example/set_operations/p4est_setop3.c


LINT_CSOURCES += $(example_set_operations_p4est_setop2_SOURCES)
LINT_CSOURCES += $(example_set_operations_p4est_setop3_SOURCES)
endif

EXTRA_DIST += \
  example/set_operations/rect_horizontal1.h \
  example/set_operations/rect_vertical1.h \
  example/set_operations/cat1.h \
  example/set_operations/cat2.h
