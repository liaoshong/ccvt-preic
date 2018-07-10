CC = gcc
CFLAGS = -Wall -fopenmp -std=c99
OBJS = ccdt.o initialize.o main.o minimize.o output.o vars.o voronoi.o
EXEC = ccvt-preic
LIBQHULL_NAME = libqhull_r.a
LIBQHULL_PATH = ./libqhull_r

$(EXEC) : $(LIBQHULL_NAME) $(OBJS)
	$(CC) $(CFLAGS) -L$(LIBQHULL_PATH) $(OBJS) -lqhull_r -lm -o $(EXEC)

$(LIBQHULL_NAME) :
	cd $(LIBQHULL_PATH) && $(MAKE)

$(OBJS) : ccdt.h configure.h initialize.h minimize.h output.h vars.h voronoi.h

clean :
	rm -f $(LIBQHULL_PATH)/*.o $(LIBQHULL_PATH)/$(LIBQHULL_NAME)
	rm -f $(EXEC) $(OBJS)
