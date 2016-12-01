# Makefile for building C stuff with GSL

# name of program to be compiled
# this must be the name of the source code file but WITHOUT any extension
TARGET = fwm11.5NewMatchingFunction


# flags for compiling
CCFLAGS=-lm -I/usr/local/include/gsl -L/usr/local/lib -lgsl -lgslcblas -g

# default compiler
CC=gcc

# specifying a default action
all: ${TARGET}


# here are make's instructions for compiling
# very simple here because there is only one source file

${TARGET}: ${TARGET}.c
	${CC} ${CCFLAGS} $< -o $@


# how to delete a build
# again, this is very simple
clean:
	rm -f ${TARGET}


