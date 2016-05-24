library("testthat")
context("hclust2 vs single linkage")

test_that("single_iris_distmat", {
   library("datasets")
   data("iris")

   d <- as.matrix(iris[,1:4])
   d[,] <- jitter(d) # otherwise we get a non-unique solution
   d <- dist(d)

   h1 <- hclust2(d, thresholdGini=1.0)
   h2 <- hclust(d, method='single')

   expect_equal(h1$merge, h2$merge)
   expect_equal(h1$order, h2$order)
})


test_that("single_iris_defaultdist", {
   library("datasets")
   data("iris")

   d <- as.matrix(iris[,2:3])
   d[,] <- jitter(d) # otherwise we get a non-unique solution

   h1 <- hclust2(objects=d, thresholdGini=1.0)
   h2 <- hclust(dist(d), method='single')

   expect_equal(h1$merge, h2$merge)
   expect_equal(h1$order, h2$order)
})
