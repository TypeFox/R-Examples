library(testthat)
library(phyloTop)
library(ape)

############################
# create some test objects
############################

trees <- rmtree(3,50) # three random trees on 50 tips
tree1 <- trees[[1]] # extract first tree
p <- phyloTop(trees)[1,] # extract first row of phyloTop output
pn <- phyloTop(trees, normalise=TRUE)[1,] # extract first row of phyloTop output

############################
# test that individual functions give same result as main phyloTop function
############################

test_that("avgLadder function gives same result as first column of phyloTop function", {
  expect_equal(avgLadder(tree1),p[[1]])
})

test_that("cherries function gives same result as second column of phyloTop function", {
  expect_equal(cherries(tree1),p[[2]])
})

test_that("colless.phylo function gives same result as third column of phyloTop function", {
  expect_equal(colless.phylo(tree1),p[[3]])
})

test_that("ILnumber function gives same result as fourth column of phyloTop function", {
  expect_equal(ILnumber(tree1),p[[4]])
})

test_that("maxHeight function gives same result as fifth column of phyloTop function", {
  expect_equal(maxHeight(tree1),p[[5]])
})

test_that("pitchforks function gives same result as sixth column of phyloTop function", {
  expect_equal(pitchforks(tree1),p[[6]])
})

test_that("sackin.phylo function gives same result as seventh column of phyloTop function", {
  expect_equal(sackin.phylo(tree1),p[[7]])
})

test_that("stairs1 function gives same result as eighth column of phyloTop function", {
  expect_equal(stairs(tree1)[[1]],p[[8]])
})

test_that("stairs2 function gives same result as eighth column of phyloTop function", {
  expect_equal(stairs(tree1)[[2]],p[[9]])
})


############################
# test that individual functions give same normalised result as main, normalised phyloTop function
############################

test_that("normalised avgLadder function gives same result as first column of normalised phyloTop function", {
  expect_equal(avgLadder(tree1, normalise=TRUE),pn[[1]])
})

test_that("normalised cherries function gives same result as second column of normalised phyloTop function", {
  expect_equal(cherries(tree1, normalise=TRUE),pn[[2]])
})

test_that("normalised colless.phylo function gives same result as third column of normalised phyloTop function", {
  expect_equal(colless.phylo(tree1, normalise=TRUE),pn[[3]])
})

test_that("normalised ILnumber function gives same result as fourth column of normalised phyloTop function", {
  expect_equal(ILnumber(tree1, normalise=TRUE),pn[[4]])
})

test_that("normalised maxHeight function gives same result as fifth column of normalised phyloTop function", {
  expect_equal(maxHeight(tree1, normalise=TRUE),pn[[5]])
})

test_that("normalised pitchforks function gives same result as sixth column of normalised phyloTop function", {
  expect_equal(pitchforks(tree1, normalise=TRUE),pn[[6]])
})

test_that("normalised sackin.phylo function gives same result as seventh column of normalised phyloTop function", {
  expect_equal(sackin.phylo(tree1, normalise=TRUE),pn[[7]])
})

