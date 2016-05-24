library("lsbclust")

## Check whether extended-precision long doubles are available
donttest <- .Machine$sizeof.longdouble == 0

## Check overall means, row and column means for consistency

context("Overall means")

## Create output object for overall means
# data("supermarkets")
# set.seed(432)
# ovl <- orc.lsbclust(data = supermarkets, margin = 3, delta = c(1, 0, 1, 0), 
#                    nclust = 12, nstart = 500, type = "overall")
# save(ovl, file = "./tests/testthat/ovl.rda")

## Load saved result
load("ovl.rda")

## Calculate new version
data("supermarkets")
set.seed(432)
ovltest <- orc.lsbclust(data = supermarkets, margin = 3, delta = c(1, 0, 1, 0), 
                    nclust = 12, nstart = 500, type = "overall")

## Check different aspects of results
test_that("cluster vector", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_identical(ovl$cluster, ovltest$cluster)
  })
test_that("cluster means", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(ovl$centers, ovltest$centers)
  })

## Remove objects
rm(ovl, ovltest)

context("Row means")

## Create output object for row means
# data("supermarkets")
# set.seed(4321)
# rows <- orc.lsbclust(data = supermarkets, margin = 3, delta = c(0, 1, 0, 0), 
#                    nclust = 9, nstart = 500, type = "rows")
# save(rows, file = "./tests/testthat/rows.rda")

## Load saved result
load("rows.rda")

## Calculate new version
data("supermarkets")
set.seed(4321)
rowstest <- orc.lsbclust(data = supermarkets, margin = 3, delta = c(0, 1, 0, 0), 
                        nclust = 9, nstart = 500, type = "rows")

## Check different aspects of results
test_that("cluster vector", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_identical(rows$cluster, rowstest$cluster)
  })
test_that("cluster means", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(rows$centers, rowstest$centers)
  })

## Remove objects
rm(rows, rowstest)

context("Column means")

## Create output object for column means
# data("supermarkets")
# set.seed(43251)
# cols <- orc.lsbclust(data = supermarkets, margin = 3, delta = c(1, 0, 1, 0), 
#                    nclust = 7, nstart = 500, type = "column")
# save(cols, file = "./tests/testthat/cols.rda")

## Load saved result
load("cols.rda")

## Calculate new version
data("supermarkets")
set.seed(43251)
colstest <- orc.lsbclust(data = supermarkets, margin = 3, delta = c(1, 0, 1, 0), 
                         nclust = 7, nstart = 500, type = "columns")

## Check different aspects of results
test_that("cluster vector", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_identical(cols$cluster, colstest$cluster)
  })
test_that("cluster means", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(cols$centers, colstest$centers)
  })

## Remove objects
rm(cols, colstest)