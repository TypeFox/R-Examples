{###############################################################################
# test=gather.R
# This file is part of the R package harvestr.
# 
# Copyright 2012 Andrew Redd
# Date: 6/2/2012
# 
# DESCRIPTION
# ===========
# unit testing for gather and other process flow functions.
# 
# LICENSE
# ========
# harvestr is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# dostats is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with 
# dostats. If not, see http://www.gnu.org/licenses/.
# 
}###############################################################################
library(testthat)
library(plyr)
context("main functions")
options(harvestr.time=FALSE)

test_that("gather is replicable", {
    seeds  <- gather(10, seed=1234)
    seeds0 <- gather(10, seed=1234)
    expect_identical(seeds, seeds0)
})
test_that("gather gives the appropriate kinds of seed.", {
    seed <- gather(1)[[1]]
    replace.seed(seed)
    expect_identical(RNGkind(), c("L'Ecuyer-CMRG", "Inversion"))
})
test_that("gather replaces seed.", {
    set.seed(123, "Mersenne-Twister", "Box-Muller")
    l <- get.seed()
    seeds <- gather(10)
    k <- get.seed()
    expect_identical(l,k)
    expect_identical(RNGkind(), c("Mersenne-Twister", "Box-Muller"))
})
test_that("farm is replicable", {
    seeds <- gather(10)
    e <- farm(seeds, rnorm(10), cache=F, time=F)
    f <- farm(seeds, rnorm(10), cache=F, time=F)
    expect_identical(e,f)
    expect_equivalent(e,f)
})
test_that('farm is indifferent to order.', {
    seeds <- gather(10)
    o <- sample(seq_along(seeds))
    e <- farm(seeds, rnorm(10))
    g <- farm(seeds[o], rnorm(10))[order(o)]
    expect_equivalent(e,g)
})
test_that("reap is reproducible", {
    seed <- gather(1)[[1]]
    x <- plant(list(1:10), list(seed))[[1]]
    a <- reap(x, sample)
    b <- withseed(seed, sample(1:10))
    expect_identical(a,b)
    expect_identical(reap(a, sample), reap(a, sample))
})
test_that("harvest", {
    seeds <- gather(10)
    e <- farm(seeds, rnorm(10))
    x <- harvest(e, sample, replace=T)
    y <- harvest(e, sample, replace=T)
    expect_equivalent(x,y)
})
test_that("Permutation", {
    seeds <- gather(10)
    e <- farm(seeds, rnorm(10))
    x <- harvest(e, sample, replace=T)
    o <- sample(seq_along(e))
    y <- harvest(e[o], sample, replace=T)[order(o)]
    expect_equivalent(x,y)
})
test_that("using with", {
  data <- farm(gather(3), data.frame(x123=runif(100), y456=rnorm(100)))
  m1 <- harvest(data, with, mean(x123))
  m2 <- lapply(data, with, mean(x123))
  expect_equivalent(m1, m2)
})
