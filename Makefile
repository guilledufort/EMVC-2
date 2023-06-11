OS_NAME := $(shell uname -s | tr A-Z a-z)

DEBUG_FLAGS = -pg
CFLAGS = -std=c99 -O3 -march=native -fstrict-aliasing -ffast-math -fomit-frame-pointer -Wall -lm
FILES = src/EMVC_project.c src/EMVC_functions_project.c
INCLUDES=

ifeq ($(OS_NAME),linux)
	CXX = gcc
else
	CXX = gcc-9
endif

all: emvc

emvc: 
		$(CXX) $(CFLAGS) $(INCLUDES) $(FILES) -o emvc-2

clean:
		rm -f emvc-2 *.o