context("Row and column names")

test_that("gettors retrieve names", {
  p <- mutaframe(a = 1:10)
  expect_that(names(p), equals("a"))
  expect_that(colnames(p), equals("a"))
  expect_that(rownames(p), is_identical_to(as.character(1:10)))
})

test_that("settors set names", {
  p <- mutaframe(a = 1:10)
  names(p) <- "b"
  expect_that(names(p), equals("b"))
  
  colnames(p) <- "c"
  expect_that(names(p), equals("c"))
  
  rownames(p) <- letters[1:10]
  expect_that(rownames(p), is_identical_to(letters[1:10]))
})

