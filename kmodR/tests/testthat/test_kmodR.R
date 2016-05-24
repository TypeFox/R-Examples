context("helper functions")

if (exists("iris") == TRUE) {
  dataset = iris[,1:4]
} else {
  dataset <- rbind(matrix(rnorm(200, mean = 1, sd = 0.3), ncol = 4),
                   matrix(rnorm(200, mean = 2, sd = 0.3), ncol = 4),
                   matrix(rnorm(200, mean = 3, sd = 0.3), ncol = 4))
}

k = 4
l = 3
C = dataset[sample(nrow(dataset),k),]
x = dataset[sample(nrow(dataset),1),]

dXC = dist_sqr_XC(dataset,C)
Q = Q_reorder(dXC)


# C_zero tests

test_that("C_zero is same width as dataset",
          {
            Cz = C_zero(dataset, k, method = "random", c_type = "clustroid")
            expect_equal(ncol(Cz),ncol(dataset))
          })

test_that("C_zero is length k",
          {
            Cz = C_zero(dataset, k, method = "random", c_type = "clustroid")
            expect_equal(nrow(Cz),k)
          })


# dist_sqr_xC tests

test_that("dist_sqr_xC is length 1",
          {
            dsxC = dist_sqr_xC(x, C)
            expect_equal(nrow(dsxC), 1)
          })

test_that("dist_sqr_xC is width 2",
          {
            dsxC = dist_sqr_xC(x, C)
            expect_equal(ncol(dsxC), 2)
          })


# dist_sqr_XC tests

test_that("dist_sqr_XC is same length as dataset",
          {
            dsXC = dist_sqr_XC(dataset, C)
            expect_equal(nrow(dsXC), nrow(dataset))
          })

test_that("dist_sqr_xC is width 2",
          {
            dsXC = dist_sqr_XC(dataset, C)
            expect_equal(ncol(dsXC), 2)
          })
