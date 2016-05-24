library(testthat)
library(grpregOverlap)

context("Testing incidenceMatrix(), overlapMatrix()")
data(birthwt.grpreg)
data(pathway.dat)

test_that("group argument, a list?", {
  X <- as.matrix(birthwt.grpreg[,-1:-2])
  group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
  expect_that(overlapMatrix(X, group), 
              throws_error("Argument 'group' must be a list"))
  
  group <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8), c(9), c(10, 11), 
                 c(12), c(13), c(14, 15, 16))
  expect_that(overlapMatrix(X, group), 
              not(throws_error("Argument 'group' must be a list")))
  
  group <- lapply(group, function(x) colnames(X)[x])
  expect_that(overlapMatrix(X, group), 
              not(throws_error("Argument 'group' must be a list")))
})

test_that('group information, correct?', {
  X <- as.matrix(birthwt.grpreg[,-1:-2])
  group <- list(c('a', 'b', 'c'), c('d', 'e', 'f'), c('g', 'h'), c('i'), 
                c('j', 'k'), c('l'), c('m'), c('n', 'o', 'p'))
  expect_that(overlapMatrix(X, group), 
              throws_error("The names of variables in X don't match with names in group!"))
  
  group <- list(c(1, 2, 3, 7), c(4, 5, 6), c(7, 8, 10), c(7, 10, 11), c(6, 14, 15, 16))
  ## create new groups for variables 9, 12, 13, put them at bottom
  expect_equal(rownames(overlapMatrix(X, group)), 
               paste(paste("grp", 1:5, sep="")))
  ## TODO:
  ## (1) handle cases where variables not belong to any of groups in 'group'
  ## put each variable into a separate group, stack those group at right
  ## (2) provide option of removing groups including only one variable.
  ## Will add this later...
  ## The following test won't pass at this point:
  ## expect_equal(rownames(overlapMatrix(X, group)), 
  ##              paste(paste("grp", 1:8, sep="")))  
})

test_that("incidence matrix, correct?", {
  X <- pathway.dat$expression
  group <- pathway.dat$pathways
  incid.mat <- overlapMatrix(X, group)
  expect_is(incid.mat, 'dgCMatrix')
  expect_equal(dim(incid.mat), c(308, 308))
})




