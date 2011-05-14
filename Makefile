PROG =	laplace

SRCS =	array_alloc.c grid.c test.c timer.c

OBJS =	array_alloc.o grid.o test.o timer.o

LIBS =	

CC = mpicc
CFLAGS = -O -Wall -W -pedantic -std=c89 -g
FC = gfortran
FFLAGS = -O -Wall -W -Wextra -pedantic -std=f2003 -fbounds-check -g -fbacktrace -fdump-core 
F90 = gfortran
F90FLAGS = $(FFLAGS)
LDFLAGS = $(CFLAGS)

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

grid.o: grid.h array_alloc.h timer.h
test.o: grid.h timer.h
