# LIBS
library(treeman)
library(testthat)

# DATA
tree <- randTree(10)
trees <- cTrees(tree, tree, tree)

# RUNNING
context('Testing \'TreeMen Class\'')
test_that(".checkTreeMen() works", {
  trees@treelist[[4]] <- "this is not a tree"
  res <- treeman:::.checkTreeMen(trees)
  expect_false(res)
})
test_that("print() works", {
  print(trees)
})
test_that("show() works", {
  show(trees)
})
test_that("as.character() works", {
  res <- as.character(trees)
  expect_true(is(res, "character"))
})
test_that("[ works", {
  expect_that(trees['ntips'], equals(30))
  expect_that(trees['ntrees'], equals(3))
})
test_that("[[ works", {
  expect_true(is(trees[[1]], "TreeMan"))
})
test_that(".cMenToMen() works", {
  trees <- treeman:::.cMenToMen(trees, trees)
  expect_true(is(trees, "TreeMen"))
})
test_that(".cMenToMan() works", {
  trees <- treeman:::.cMenToMan(trees, tree)
  expect_true(is(trees, "TreeMen"))
})
test_that(".cMenToAny() works", {
  trees <- treeman:::.cMenToAny(trees, tree)
  expect_true(is(trees, "TreeMen"))
})
test_that(".cTreeObjs() works", {
  trees <- treeman:::.cTreeObjs(trees, tree, tree, trees)
  expect_true(is(trees, "TreeMen"))
})
test_that("setAs() works", {
  trees <- list(tree, tree, tree)
  trees <- as(trees, "TreeMen")
  expect_true(is(trees, "TreeMen"))
})