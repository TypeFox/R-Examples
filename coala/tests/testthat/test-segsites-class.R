context("SegSites Class")

test_that("normal initialization of segsites works", {
  snps <- matrix(c(1, 1, 0,  0, 1, 0,  0, 1, 1), 3, 3)
  pos <- c(.1, .3, .5)

  segsites <- create_segsites(snps, pos)
  expect_true(is_segsites(segsites))
  expect_equal(dim(segsites), c(3, 3))
  expect_equal(nrow(segsites), 3)
  expect_equal(ncol(segsites), 3)
  expect_equal(get_positions(segsites), pos)
  expect_equal(get_trio_locus(segsites), c(0, 0, 0))

  expect_error(create_segsites(snps, 1:4 / 10))
  expect_error(create_segsites(snps, pos, 1:4 / 10))
})


test_that("fixed positions are removed on initialization", {
  snps <- matrix(c(1, 1, 0,  1, 1, 1,  0, 0, 0), 3, 3)
  segsites <- create_segsites(snps, 1:3 / 10)
  expect_equal(dim(segsites), c(3, 1))
  expect_equal(get_positions(segsites), 0.1)
  expect_equivalent(as.matrix(segsites), matrix(c(1, 1, 0), 3, 1))
})


test_that("subsetting of segsites works", {
  segsites <- create_segsites((matrix(c(1, 1, 0, 1, 1,
                                        0, 0, 0, 1, 0,
                                        0, 1, 1, 0, 1), 3, 5, byrow = TRUE)),
                              c(.1, .2, .5, .7, .75))

  expect_equal(dim(segsites), c(3, 5))
  ss_subset <- segsites[ , 1:2]
  expect_equal(dim(ss_subset), c(3, 2))
  expect_equivalent(as.matrix(ss_subset), as.matrix(segsites)[ , 1:2])
  expect_equal(get_positions(ss_subset), get_positions(segsites)[1:2])

  expect_equal(dim(segsites[ , 1]), c(3, 1))
  expect_equal(dim(segsites[1:2 , ]), c(2, 3))
  expect_equal(dim(segsites[-2 , 1:3]), c(2, 2))
})


test_that("segsites can be converted to ms-output", {
  segsites <- create_test_segsites()
  output <- conv_to_ms_output(segsites)
  expect_equal(output, c("segsites: 5",
                         "positions: 0.10 0.20 0.50 0.70 0.75",
                         "11011",
                         "00010",
                         "01101"))
})


test_that("converting segsites to matrix works", {
  x <- create_test_segsites()
  mat <- as.matrix(x)
  expect_true(is.matrix(mat))
  expect_false(is_segsites(mat))
})
