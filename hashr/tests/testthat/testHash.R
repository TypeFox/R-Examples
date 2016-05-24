
context("hashing")

test_that("",{
  expect_true(is.integer(hash(10)))
  expect_equal(length(hash(c("aap","noot"))),2L)
  expect_equal(length(hash(c("aap","noot"),recursive=FALSE)),1L)
  expect_equal(length(hash(list(1,c("aap","noot")))),2L)
  expect_equal(length(hash(list(1,c("aap","noot")),recursive=FALSE)),1L)
  m <- lm(Sepal.Width ~ Sepal.Length, data=iris)
  expect_true(is.integer(hash(m)))
  expect_equal(length(hash(m)),1L)
  x <- c("call any vegetable","and the chances are good","that the vegetable will respond to you")
  L <- strsplit(x," ",fixed=TRUE)
  expect_true(all(sapply(hash(L),is.integer)))
})




