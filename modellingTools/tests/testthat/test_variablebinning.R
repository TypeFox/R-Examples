require(modellingTools, quietly = TRUE, warn.conflicts = FALSE)
require(dplyr, quietly = TRUE, warn.conflicts = FALSE)

context("Variable Binning")

### vector_bin

x <- seq(0.01,9,by = 0.01)
x <- c(x,rep(NA,100))

vbx_10 <- vector_bin(x,bins = 10)
vbx_w10 <- vector_bin(x,bins = 10,type = "width")

vbx_nax10 <- vector_bin(x,bins = 10,na_include = FALSE)
vbx_wnax10 <- vector_bin(x,bins = 10,type = "width",na_include = FALSE)

too_many_bins <- vector_bin(x,bins = 10000)
custom_bins <- vector_bin(x,c(0.1,0.3,0.5,0.7))

test_that("vector_bin returns correct number of bins", {
  expect_equal(length(levels(vbx_10)),11)
  expect_equal(length(levels(vbx_w10)),11)

  expect_equal(length(levels(vbx_nax10)),10)
  expect_equal(length(levels(vbx_wnax10)),10)
})

test_that("vector_bin returns the original vector if there are too many bins asked for", {
  expect_identical(too_many_bins,x)
})

test_that("regardless of value of na_include, there are no na's in the output", {
  expect_false(any(is.na(vbx_10)))
  expect_false(any(is.na(vbx_w10)))

  expect_false(any(is.na(vbx_nax10)))
  expect_false(any(is.na(vbx_wnax10)))
})

# test_that("custom binnning works", {
#   expect_identical(vbx_10,
#                    vector_bin(x,bins = get_vector_cutpoints(vbx_10)))
#   expect_identical(vbx_w10,
#                    vector_bin(x,bins = get_vector_cutpoints(vbx_w10)))
#
#   expect_identical(vbx_nax10,
#                    vector_bin(x,bins = get_vector_cutpoints(vbx_nax10),na_include = FALSE))
#   expect_identical(vbx_wnax10,
#                    vector_bin(x,bins = get_vector_cutpoints(vbx_wnax10),na_include = FALSE))
# })

test_that("custom binning provides the correct factor levels", {
  expect_identical(levels(vbx_10),
                   levels(vector_bin(x,bins = get_vector_cutpoints(vbx_10))))
  expect_identical(levels(vbx_w10),
                   levels(vector_bin(x,bins = get_vector_cutpoints(vbx_w10))))

  expect_identical(levels(vbx_nax10),
                   levels(vector_bin(x,bins = get_vector_cutpoints(vbx_nax10),na_include = FALSE)))
  expect_identical(levels(vbx_wnax10),
                   levels(vector_bin(x,bins = get_vector_cutpoints(vbx_wnax10),na_include = FALSE)))
})

### get_vector_cutpoints

custom_bins <- c(1,2,3)
gvx <- cut(x,custom_bins)

test_that("get_vector_cutpoints returns the correct type", {
  expect_true(is.numeric(get_vector_cutpoints(gvx)))
})

test_that("get_vector_cutpoints returns the correct values", {
  expect_equal(get_vector_cutpoints(gvx),c(1,2,3))
  expect_equal(get_vector_cutpoints(c(-1,1)),c(-1,1))
  expect_equal(get_vector_cutpoints(c("123")),c(123))
})

test_that("scientific notation works", {
  expect_equal(get_vector_cutpoints(c(1e06,1e-06,-1e06,-1e-06)),
               c(-1e06,-1e-06,1e-06,1e06))
})

### binned_data_cutpoints

d <- data_frame(v1 = cut(x,c(1,2,3)),
                v2 = cut(-x,c(-1,-2,-3)),
                v3 = x,
                v4 = factor(x))

bdc <- binned_data_cutpoints(d)

test_that("binned_data_cutpoints returns a named list with apprpriate names", {
  expect_true(is.list(bdc))
  expect_named(bdc,
               expected = c("v1","v2","v4"),
               ignore.order = TRUE,
               ignore.case = TRUE)
})

test_that("binned_data_cutpoints returns the correct cutpoints", {
  expect_equal(bdc$v1,c(1,2,3))
  expect_equal(bdc$v2,c(-3,-2,-1))
  expect_equal(bdc$v4,x[!is.na(x)])
})

### simple_bin

tr <- data_frame(v1 = seq(0,9,by = 1),
                 v2 = letters[1:10],
                 v3 = factor(v2),
                 id = v1)
tt <- tr

bin1 <- simple_bin(tr,
                   tt,
                   exclude_vars = "id",
                   bins = 3)
trbin1 <- simple_bin(tr,exclude_vars = "id",
                     bins = 3)
ttbin1 <- simple_bin(tt,bins = binned_data_cutpoints(trbin1),exclude_vars = "id")

bin2 <- simple_bin(tr,
                   tt,
                   exclude_vars = "id",
                   bins = 3,
                   type = "width")
trbin2 <- simple_bin(tr,exclude_vars = "id",type = "width",
                     bins = 3)
ttbin2 <- simple_bin(tt,bins = binned_data_cutpoints(trbin2),exclude_vars = "id")

bin3 <- simple_bin(tr,
                   tt,
                   exclude_vars = "id",
                   na_include = FALSE,
                   bins = 3)
trbin3 <- simple_bin(tr,exclude_vars = "id",na_include = FALSE,
                     bins = 3)
ttbin3 <- simple_bin(tt,bins = binned_data_cutpoints(trbin3),exclude_vars = "id")

test_that("simple_bin returns factors for all numeric and factor columns", {
  expect_true(is.factor(column_vector(trbin1,"v1")))
  expect_true(is.factor(column_vector(trbin1,"v3")))
  expect_true(is.factor(column_vector(ttbin1,"v1")))
  expect_true(is.factor(column_vector(ttbin1,"v3")))

  expect_true(is.factor(column_vector(trbin2,"v1")))
  expect_true(is.factor(column_vector(trbin2,"v3")))
  expect_true(is.factor(column_vector(ttbin2,"v1")))
  expect_true(is.factor(column_vector(ttbin2,"v3")))

  expect_true(is.factor(column_vector(trbin3,"v1")))
  expect_true(is.factor(column_vector(trbin3,"v3")))
  expect_true(is.factor(column_vector(ttbin3,"v1")))
  expect_true(is.factor(column_vector(ttbin3,"v3")))

})

test_that("simple_bin still returns non-binnable columns", {
  expect_true(all(c("v2","id") %in% colnames(trbin1)))
  expect_true(all(c("v2","id") %in% colnames(ttbin1)))

  expect_true(all(c("v2","id") %in% colnames(trbin2)))
  expect_true(all(c("v2","id") %in% colnames(ttbin2)))

  expect_true(all(c("v2","id") %in% colnames(trbin3)))
  expect_true(all(c("v2","id") %in% colnames(ttbin3)))
})

test_that("binning on the test set works as expected", {
  expect_identical(binned_data_cutpoints(trbin1),
                   binned_data_cutpoints(ttbin1))

  expect_identical(binned_data_cutpoints(trbin2),
                   binned_data_cutpoints(ttbin2))

  expect_identical(binned_data_cutpoints(trbin3),
                   binned_data_cutpoints(ttbin3))
})

test_that("test sets are same from manual and automatic generation", {
  expect_identical(bin1$test,ttbin1)
  expect_identical(bin2$test,ttbin2)
  expect_identical(bin3$test,ttbin3)

})


