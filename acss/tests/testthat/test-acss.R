
context("acss")

test_that("acss returns what it should return", {
  res1 <- structure(c(22.0030069458495, NA, 24.75269450012, 24.9239869906869, 
      2.37922171238823e-07, NA, 3.53750027225721e-08, 3.14146556992709e-08
      ), .Dim = c(2L, 4L), .Dimnames = list(c("01011100", "00030101"
      ), c("K.2", "K.4", "D.2", "D.4")))
  
  expect_equal(acss(c("01011100", "00030101"), alphabet = c(2, 4)), res1)
  
  # taken from the AcssGuide.pdf
  res2 <- acss("010120123", alphabet = 5)[1,]
  expect_equivalent(round(res2[1], 5), c(30.03474))
  expect_equivalent(res2[2], 9.091617e-10)
})

test_that("acss structural aspects", {
  expect_is(acss("01011100"), "matrix")
  expect_identical(nrow(acss("01011100")), 1L)
  # test if one string works like 2 strings:
  two <- acss(c("01011100", "01011100"))
  one <- acss(c("01011100"))
  expect_true(all(rownames(two) == rownames(one)))
  expect_equal(two[1,], one[1,])
  expect_equal(two[2,], one[1,])
  vector <- acss(c("010", "010"))[,2]
  expect_true(is.vector(vector, "numeric"))
  expect_named(vector)
})

context("local_complexity")

test_that("local_complexity structural", {
  l1 <- local_complexity(c("01011010111"), span=5, alphabet = 5)
  expect_that(l1, is_a("list"))
  expect_that(length(l1), is_identical_to(1L))
  expect_named(l1)
  l2 <- local_complexity(c("01011010111" , "01011010111"), span=5, alphabet = 5)
  expect_that(l2, is_a("list"))
  expect_that(length(l2), is_identical_to(2L))
  expect_named(l2)
  expect_that(l2[[1]], equals(l2[[2]]))
  expect_that(l2[[1]], equals(l1[[1]]))
  l2b <- local_complexity(c("01011010111" ,"GHHGGHGHUE"),span=7)
  expect_false(length(l2b[[1]]) == length(l2b[[2]]))
  expect_false(isTRUE(all.equal(l2b[[1]], l2b[[2]])))
})

test_that("local_complexity: guide", {
  expect_equal(round(mean(local_complexity("010120123", alphabet = 5)[[1]]), 5),  16.81957)  
})



context("prob_random")

test_that("prob_random structural", {
  l1 <- prob_random(c("HEHHEE"))
  expect_named(l1)
  expect_true(is.vector(l1, "numeric"))
  expect_identical(length(l1), 1L)
  l2 <- prob_random(c("HEHHEE", "HEHHEE"))
  expect_named(l2)
  expect_identical(length(l2), 2L)
  expect_true(is.vector(l2, "numeric"))
  expect_identical(l2[1], l2[2])
  expect_identical(l2[1], l1)
  l3 <- prob_random(c("HEHHEE", "GHHGGHGHUE", "HSHSHHSHSS"))
  expect_named(l3)
  expect_identical(length(l3), 3L)
  expect_true(is.vector(l3, "numeric"))
  expect_that(l3[1], is_more_than(l3[2]))
  expect_that(l3[1], is_more_than(l3[3]))
  expect_that(l3[2], is_more_than(l3[3]))
})



