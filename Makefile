OS_NAME := $(shell uname -s | tr A-Z a-z)

DEBUG_FLAGS = -pg
CFLAGS = -std=c99 -O3 -march=native -fstrict-aliasing -ffast-math -fomit-frame-pointer -Wall -lm
FILES = src/EMVC_project.c src/EMVC_functions_project.c
INCLUDES=

ifeq ($(OS_NAME),linux)
	CC = gcc
else
	CC = gcc
endif

all: candidate_variants_finder

candidate_variants_finder: 
		$(CC) $(CFLAGS) $(INCLUDES) $(FILES) -o candidate_variants_finder

clean:
		rm -f candidate_variants_finder *.o