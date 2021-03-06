#
# ARGweaver
# Matt Rasmussen
# Copyright 2012-2013
#
# Makefile
#

# install prefix paths
prefix = /usr


# programs
ifdef MPI
	CXX=mpic++
        CFLAGS := $(CFLAGS) -DARGWEAVER_MPI
else
        CXX=g++
endif
PYTHON = python2

# C++ compiler options
CFLAGS := $(CFLAGS) \
    -Wall -fPIC \
    -Isrc

GTEST_URL = 'http://googletest.googlecode.com/files/gtest-1.7.0.zip'
GTEST_DIR = gtest-1.7.0

#=============================================================================
# optional CFLAGS

# profiling
ifdef PROFILE
	CFLAGS := $(CFLAGS) -pg
endif

# debugging
ifdef DEBUG
	CFLAGS := $(CFLAGS) -g -DDEBUG
else
	CFLAGS := $(CFLAGS) -O3 -funroll-loops -Wno-strict-overflow
endif


#=============================================================================
# ARGweaver files

# package
PKG_VERSION:=$(shell $(PYTHON) -c 'import argweaver; print argweaver.PROGRAM_VERSION_TEXT' 2>/dev/null)
PKG_NAME=argweaver
DIST=dist
PKG=$(DIST)/$(PKG_NAME)-$(PKG_VERSION).tar.gz
PKG_DIR=$(DIST)/$(PKG_NAME)-$(PKG_VERSION)

# program files
SCRIPTS = bin/*
PROGS = bin/arg-sample bin/arg-likelihood bin/arg-summarize bin/smc2bed 
BINARIES = $(PROGS) $(SCRIPTS)

ARGWEAVER_SRC = $(shell ls src/argweaver/*.cpp)

#    src/prior.cpp
#    src/expm/matrix_exponential.cpp
#    src/expm/r8lib.cpp


ALL_SRC = \
    $(ARGWEAVER_SRC) \
    src/arg-sample.cpp \
    src/arg-summarize.cpp \
    src/smc2bed.cpp \
    src/popsize-post.cpp \
    src/compress-sites.cpp \
    src/arg-likelihood.cpp


ARGWEAVER_OBJS = $(ARGWEAVER_SRC:.cpp=.o)
ALL_OBJS = $(ALL_SRC:.cpp=.o)

LIBS =
# `gsl-config --libs`
#-lgsl -lgslcblas -lm


# ARGweaver C-library files
LIBARGWEAVER = lib/libargweaver.a
LIBARGWEAVER_SHARED_NAME = libargweaver.so
LIBARGWEAVER_SHARED = lib/$(LIBARGWEAVER_SHARED_NAME)
LIBARGWEAVER_SHARED_INSTALL = $(prefix)/lib/$(LIBARGWEAVER_SHARED_NAME)
LIBARGWEAVER_OBJS = $(ARGWEAVER_OBJS)


# ARGweaver C++ tests
CFLAGS_TEST = -I $(GTEST_DIR)/include
LIBS_TEST = -Llib -lgtest -lgtest_main -lpthread
GTEST_SRC = gtest-1.7.0
TEST_SRC = \
	src/tests/test.cpp \
	src/tests/test_local_tree.cpp \
	src/tests/test_prob.cpp

TEST_OBJS = $(TEST_SRC:.cpp=.o)

#=============================================================================
# tskit files
TSK_OBJECTS=kastore.o tskit_tables.o tskit_core.o tskit_trees.o \
        tskit_stats.o tskit_genotypes.o tskit_convert.o

CFLAGS_tskit += -Ikastore/c -I tskit/c
CC = gcc
#-----------------------------
# tskit component

libtskit.a: ${TSK_OBJECTS}
	${AR} rcs $@ ${TSK_OBJECTS} 

kastore.o: kastore
	${CC} -c ${CFLAGS_tskit} kastore/c/kastore.c -o kastore.o

tskit_%.o: tskit
	${CC} -c ${CFLAGS_tskit} tskit/c/tskit/$*.c -o $@

kastore:
	git clone https://github.com/tskit-dev/kastore.git
	# NB!!! Make sure to checkout at a version tag!
	cd kastore && git checkout C_2.0.0

tskit: 
	git clone https://github.com/tskit-dev/tskit.git
	# NB!!! Make sure to checkout at a version tag!
	cd tskit && git checkout 0.2.3

CFLAGS += -Ikastore/c -Itskit/c
CFLAGS += -L. -ltskit

#=============================================================================

# rSPR component

# TODO: need to figure out why this part of makefile doesn't work!!!

# rSPR_OBJECTS=rspr.o spr_supertree.o fill_matrix.o
# CFLAGS_rSPR=-O2 -std=c++0x -march=native


# librspr.a: $(rSPR_OBJECTS)
# 	$(AR) rcs $@ $(rSPR_OBJECTS) 

# rspr.o: rspr
# 	$(CXX) -c $(CFLAGS_rSPR) rspr/rspr.cpp -o rspr.o

# spr_supertree.o: rspr
# 	$(CXX) -c $(CFLAGS_rSPR) rspr/spr_supertree.cpp -o spr_supretree.o

# fill_matrix.o: rspr
# 	$(CXX) -c $(CFLAGS_rSPR) -o rspr/fill_matrix.cpp -o fill_matrix.o

# rspr:
# 	git clone https://github.com/cwhidden/rspr.git

CFLAGS += -Irspr -lrspr -L./rspr

#==============================================================================
# targets

.PHONY: all pkg test ctest cq install clean cleanobj lib pylib gtest

# default targets



all: $(PROGS) $(LIBARGWEAVER) $(LIBARGWEAVER_SHARED)

bin/arg-sample: src/arg-sample.o $(LIBARGWEAVER)
	$(CXX) -o bin/arg-sample src/arg-sample.o $(LIBARGWEAVER) $(CFLAGS)

bin/smc2bed: src/smc2bed.o $(LIBARGWEAVER)
	$(CXX) -o bin/smc2bed src/smc2bed.o $(LIBARGWEAVER) $(CFLAGS)

bin/arg-summarize: src/arg-summarize.o $(LIBARGWEAVER)
	$(CXX) -o bin/arg-summarize src/arg-summarize.o $(LIBARGWEAVER) $(CFLAGS)

bin/popsize-post: src/popsize-post.o $(LIBARGWEAVER)
	$(CXX) -o bin/popsize-post src/popsize-post.o $(LIBARGWEAVER) $(CFLAGS)

bin/compress-sites: src/compress-sites.o $(LIBARGWEAVER)
	$(CXX) -o bin/compress-sites src/compress-sites.o $(LIBARGWEAVER) $(CFLAGS)

bin/arg-likelihood: src/arg-likelihood.o $(LIBARGWEAVER)
	$(CXX) -o bin/arg-likelihood src/arg-likelihood.o $(LIBARGWEAVER) $(CFLAGS)

#-----------------------------
# ARGWEAVER C-library
lib: $(LIBARGWEAVER) $(LIBARGWEAVER_SHARED)

$(LIBARGWEAVER): $(LIBARGWEAVER_OBJS)
	mkdir -p lib
	$(AR) -r $(LIBARGWEAVER) $(LIBARGWEAVER_OBJS)

$(LIBARGWEAVER_SHARED): $(LIBARGWEAVER_OBJS)
	mkdir -p lib
	$(CXX) -o $(LIBARGWEAVER_SHARED) -shared $(LIBARGWEAVER_OBJS) $(LIBS)


#-----------------------------
# packaging

pkg: $(PKG)

$(PKG):
	mkdir -p $(DIST)
	git archive --format=tar --prefix=$(PKG_NAME)-$(PKG_VERSION)/ HEAD | \
	gzip > $(PKG)

#-----------------------------
# testing

test: $(LIBARGWEAVER_SHARED)
	nosetests -v test

cq:
	nosetests -v test/test_codequality.py

ctest: src/tests/test
	src/tests/test

src/tests/test: $(TEST_OBJS) $(LIBARGWEAVER)
	$(CXX) -o src/tests/test $(TEST_OBJS) $(LIBS_TEST) $(LIBARGWEAVER)

$(TEST_OBJS): %.o: %.cpp
	$(CXX) -c $(CFLAGS) $(CFLAGS_TEST) -o $@ $<

# Download and install gtest unit-testing framework.
gtest:
	wget $(GTEST_URL) -O gtest.zip
	rm -rf $(GTEST_SRC)
	unzip gtest
	cd $(GTEST_SRC) && ./configure && make
	cp $(GTEST_SRC)/lib/.libs/libgtest_main.a lib/
	cp $(GTEST_SRC)/lib/.libs/libgtest.a lib/

#-----------------------------
# install

install: $(BINARIES)
	mkdir -p $(prefix)/bin
	cp $(BINARIES) $(prefix)/bin
	$(PYTHON) setup.py install --prefix=$(prefix)

pylib: $(LIBARGWEAVER_SHARED_INSTALL)
	$(PYTHON) setup.py install --prefix=$(prefix)


$(LIBARGWEAVER_SHARED_INSTALL): $(LIBARGWEAVER_SHARED)
	mkdir -p $(prefix)/lib
	cp $(LIBARGWEAVER_SHARED) $(LIBARGWEAVER_SHARED_INSTALL)

#=============================================================================
# basic rules

$(ALL_OBJS): %.o: %.cpp
	$(CXX) -c $(CFLAGS) -o $@ $<

clean:
	rm -f $(ALL_OBJS) $(LIBARGWEAVER) $(LIBARGWEAVER_SHARED) $(TEST_OBJS) $(PROGS)

clean-test:
	rm -f $(TEST_OBJS)

clean-obj:
	rm -f $(ALL_OBJS) $(TEST_OBJS)
