###############################################################
# This Makefile is written by Shihong Liao to create the static
# library named libqhull_r.a, which is used in compiling the
# main programs in the parent folder (../). It is different
# from the original Makefile included in /libqhull_r by C.B. 
# Barber.
#
# E-mail: liaoshong@gmail.com
# Date: 2017/08/27
###############################################################

CC = gcc
OBJS = geom2_r.o geom_r.o global_r.o io_r.o libqhull_r.o \
       mem_r.o merge_r.o poly2_r.o poly_r.o qset_r.o random_r.o \
       rboxlib_r.o stat_r.o usermem_r.o userprintf_r.o \
       userprintf_rbox_r.o user_r.o
LIBQHULL_NAME = libqhull_r.a

$(LIBQHULL_NAME) : $(OBJS)
	ar -r $(LIBQHULL_NAME) $(OBJS)

$(OBJS) : geom_r.h io_r.h libqhull_r.h mem_r.h merge_r.h poly_r.h \
          qhull_ra.h qset_r.h random_r.h stat_r.h user_r.h

clean :
	rm -f *.o $(LIBQHULL_NAME)
