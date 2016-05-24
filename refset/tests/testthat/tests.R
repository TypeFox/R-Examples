
library(testthat)

context("Basic refset")

test_that("refsets refer to original", {
  dfr <- data.frame(a=1:5, b=6:10)
  refset(ss, dfr, 1:3, 1:2)
  expect_equivalent(ss$b, 6:8)
  dfr$b <- 1:5
  expect_equivalent(ss$b, 1:3)
})

test_that("refsets change original", {
  dfr <- data.frame(a=1:5, b=6:10)
  refset(ss, dfr, 1:3, 1:2)
  ss$b <- 1:3
  expect_equivalent(dfr$b, c(1:3,9:10))
})

test_that("subsets of refsets work and change original", {
  dfr <- data.frame(a=1:5, b=6:10)
  refset(ss, dfr, 3:5, 1:2)
  dfr$a <- 6:10
  expect_equivalent(ss[1:2,1], 8:9)
  expect_equivalent(ss[[3,1]], 10L)
  expect_equivalent(ss[c(TRUE,TRUE,FALSE),], dfr[3:4,])
  ss$a <- 1:3
})

# Feature 3: subsets of subsets (etc...) refer to original
test_that("refsets of refsets work", {
  dfr <- data.frame(a=1:5, b=6:10)
  refset(ss, dfr, 3:5, 1:2)
  refset(sss, ss, 1:2, 1)
  expect_equivalent(sss, 3:4)
  dfr$a <- 6:10
  expect_equivalent(sss, 8:9)
  ss$a <- 1:3
  expect_equivalent(sss, 1:2)
})

test_that("refsets of refsets change original", {
  dfr <- data.frame(a=1:5, b=6:10)
  refset(ss, dfr, 3:5, 1:2)
  refset(sss, ss, 1:2, 1)
  expect_equivalent(sss, dfr$a[3:4])
  sss <- 8:9
  expect_equivalent(dfr$a[3:4], 8:9)
  expect_equivalent(ss$a[1:2], 8:9)
})

test_that("reassignment breaks link", {
  dfr <- data.frame(a=1:5, b=6:10)
  refset(ss, dfr, 1:3, 1:2)
  copyss <- ss 
  dfr$b <- 1:5
  expect_equivalent(ss$b,  1:3)
  expect_equivalent(copyss$b, 6:8)
})

test_that("name-based indexing works", {
  dfr <- data.frame(a=1:5, b=6:10)  
  refset(ss, dfr, "a")
  expect_equivalent(ss$a[1:3], 1:3)
  dfr$a <- 5:1
  expect_equivalent(ss$a[1:3], 5:3)
  ss$a <- 1:5
  expect_equivalent(dfr$a, 1:5)
})

test_that("dynamic indexing works", {
  dfr <- data.frame(a=1:5, b=6:10)  
  refset(ss, dfr, dfr$a <= 3, 1:2, drop=FALSE)
  refset(ss.static, dfr, dfr$a <= 3, 1:2, dyn.idx=FALSE)
  expect_equivalent(nrow(ss), 3L)
  dfr$a <- 2:6
  expect_equivalent(nrow(ss), 2L)
  ss$a <- 3:4
  expect_equivalent(nrow(ss),1L)
  expect_equivalent(nrow(ss.static), 3L)
  expect_equivalent(ss.static$a, c(3L, 4L, 4L))
})

test_that("data frame indexing works", {
  dfr <- data.frame(a=1:5, b=6:10)  
  refset(ss, dfr, a <= 3, 1:2)
  expect_equivalent(ss$b, 6:8)
  dfr$a <- rep(1,5)
  expect_equivalent(ss$b, 6:10)
})

test_that("assigning into an environment works", {
  e2 <- new.env(emptyenv())
  assign("dfr", data.frame(a=1:5, b=6:10), e2)
  refset(ss, dfr, 1:3, 1:2, assign.env=e2, eval.env=e2)
  expect_equivalent(e2$ss$a, 1:3)
  e2$dfr$a <- 6:10
  expect_equivalent(e2$ss$a, 6:8)
  e2$ss$a <- 1:3
  expect_equivalent(e2$dfr$a, c(1:3, 9:10))
})

test_that("empty arguments work", {
  dfr <- data.frame(a=1:5, b=6:10, c=1:5)  
  refset(ss, dfr, 1:3, )
  refset(ss2, dfr, , c(1,3))
  refset(ss2.static, dfr, , c(1,3), dyn.idx=FALSE)
  expect_equivalent(dim(ss), c(3L, 3L))
  expect_equivalent(dim(ss2), c(5L, 2L))
  
  dfr <- rbind(dfr, rep(10,3))
  expect_equivalent(dim(ss2), c(6L, 2L))
  expect_equivalent(dim(ss2.static), c(5L, 2L))
})

test_that("drop works", {
  dfr <- data.frame(a=1:5, b=6:10)  
  refset(ssd, dfr, ,"a", drop=TRUE)  
  refset(ssnd, dfr, ,"a", drop=FALSE)  
  expect_null(dim(ssd))
  expect_equivalent(dim(ssnd), c(5L, 1L)) 
  
  refset(ssd, dfr, 1, , drop=TRUE)
  refset(ssnd, dfr, 1, , drop=FALSE)
  expect_is(ssd, "list")
  expect_is(ssnd, "data.frame")
  
  refset(ssd, dfr, 1:2, , drop=TRUE)
  refset(ssnd, dfr, 1:2, , drop=FALSE)
  expect_equivalent(ssd, ssnd)
})

# Feature 7: any kind of subsettable data within a view
test_that("vectors work", {
  vec <- 1:5
  refset(ssv, vec, 1:2)
  expect_equivalent(ssv, 1:2)
  vec[2:3] <- 1
  expect_equivalent(ssv, c(1,1))
  ssv <- 4:5
  expect_equivalent(vec, c(4:5, 1, 4:5))
})

test_that("lists work", {
  lst <- list(a=1, b="foo", "bar")
  refset(ssl, lst, c('a', 'b'))
  refset(ssl2, lst, 3)
  refset(ssl3, lst, sapply(lst, is.character))
  
  lst$a <- 10
  expect_equal(ssl$a, 10)
  lst[[3]] <- "baz"
  expect_equal(ssl2[[1]], "baz")
  expect_equal(ssl3[[2]], "baz")
  lst$a <- "bop"
  expect_equivalent(length(ssl3), 3L)
  ssl3[[1]] <- "bip"
  expect_equal(lst$a, "bip")
  ssl3[[1]] <- NA
  expect_equivalent(length(ssl3), 2L)
})

test_that("matrices work", {
  mat <- matrix(1:9, nrow=3)
  refset(ssm, mat, 1:2, 1:2)
  expect_equivalent(ssm, matrix(c(1:2, 4:5), nrow=2))
  mat[1:2,1:2] <- 4
  expect_equivalent(ssm, matrix(4, ncol=2, nrow=2))
  ssm[] <- 2:5
  expect_equal(mat[1:2,1:2], matrix(2:5, nrow=2)) # not identical (int v dbl)
  
  rownames(mat) <- LETTERS[1:3]
  expect_equivalent(rownames(ssm), LETTERS[1:2])
  
  mat[] <- 1:9
  colnames(mat) <- LETTERS[24:26]
  refset(ssm2, mat, LETTERS[1:2], LETTERS[25:26])
  rownames(mat) <- LETTERS[3:1]
  expect_equivalent(ssm2, mat[3:2, 2:3])
  
  skip("can't set rownames of subset")
  rownames(ssm) <- LETTERS[25:26] 
  expect_equivalent(rownames(ssm), LETTERS[25:26])
  rownames(mat) == LETTERS[c(25, 26, 3)] # no dice!
})

test_that("views work in all environments", {
  g <- function() {
    dfr <- data.frame(a=1:5, b=6:10)
    refset(ss, dfr, 1:3, 1:2)
    dfr$b <- 1:5
    expect_equivalent(ss$b, 1:3)
    ss$a <- 6:8 
    expect_equivalent(dfr$a, c(6:8, 4:5))
  } 
  g()
})

test_that("dynamic indexing works in all environments", {
  g <- function() {
    dfr <- data.frame(a=1:5, b=6:10)
    refset(ss, dfr, dfr$a <= 3, 1:2)
    expect_equivalent(ss$b, 6:8)
    dfr$a <- 5:1
    expect_equivalent(ss$b, 8:10)
  } 
  g()
})

test_that("two-argument form works", {
  dfr <- data.frame(a=1:4, b=5:8)
  refset(ss, dfr[1,1])
  expect_equivalent(ss, 1L)
  refset(ss2, dfr[c(FALSE, TRUE), "b"])
  expect_equivalent(ss2, c(6L, 8L))
  dfr$a <- 4:1
  dfr$b <- 8:5
  expect_equivalent(ss, 4L)
  expect_equivalent(ss2, c(7L, 5L))
  ss <- 10L
  expect_equivalent(dfr[1,1], 10L)
  ss2 <- c(NA, NA)
  expect_equivalent(dfr$b, c(8L, NA, 6L, NA))

})

test_that("two-argument form with [[ works", {
  dfr <- data.frame(a=1:4, b=5:8)
  refset(ss, dfr[["a"]])
  expect_equivalent(ss, 1:4)
  dfr$a <- 4:1
  expect_equivalent(ss, 4:1)
  ss <- 1:4
  expect_equivalent(dfr$a, 1:4)
})

test_that("two-argument form with $ works", {
  dfr <- data.frame(a=1:4, b=5:8)
  refset(ss, dfr$a)
  expect_equivalent(ss, 1:4)
  dfr$a <- 4:1
  expect_equivalent(ss, 4:1)
  ss <- 1:4
  expect_equivalent(dfr$a, 1:4)
})


test_that("complex two-argument forms work", {
  dfr <- data.frame(a=1:4, b=5:8, c=9:12)
  refset(ss, dfr$b[2:3])
  expect_equivalent(ss, 6:7)
  dfr$b <- LETTERS[1:4]
  expect_equivalent(ss, LETTERS[2:3])
  ss <- LETTERS[25:26]
  expect_equivalent(dfr$b, c("A", "Y", "Z", "D"))
  
  dfr <- data.frame(a=1:4, b=5:8, c=9:12)
  refset(ss2, dfr[c(TRUE, FALSE), c("a", "c")])
  expect_equivalent(ss2, dfr[c(1,3), c("a", "c")])
  dfr$a <- 4:1
  expect_equivalent(ss2, dfr[c(1,3), c("a", "c")])
  ss2 <- 10
  expect_equivalent(dfr$a, c(10, 3, 10, 1))
})

test_that("changing nrow of original object works", {
  dfr <- data.frame(a=1:8, b=1:8)
  refset(ss, dfr[1:(nrow(dfr)/2),])
  expect_equivalent(nrow(ss), 4L)
  dfr <- dfr[1:4,]
  expect_equivalent(nrow(ss), 2L)
  dfr <- rbind(dfr, dfr, dfr)
  expect_equivalent(nrow(ss), 6L)
})

test_that("changing ncol of original object works", {
  mx <- matrix(1:16, ncol=4)
  refset(ss, mx[, c(FALSE, TRUE)])
  expect_equivalent(ncol(ss), 2L)
  mx <- cbind(mx, mx)
  expect_equivalent(ncol(ss), 4L)
})


