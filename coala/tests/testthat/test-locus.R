context("Locus Class")

test_that("Creating loci works", {
  locus <- locus_class$new(1013)
  expect_true(is.locus(locus))
  expect_equal(locus$get_length(), 1013)
  expect_equal(locus$get_number(), 1)

  locus <- locus_class$new(1014, 10)
  expect_true(is.locus(locus))
  expect_equal(locus$get_length(), 1014)
  expect_equal(locus$get_number(), 10)

  expect_error(locus_class$new("abc", 10))
  expect_error(locus_class$new(-5, 10))
  expect_error(locus_class$new(10, "10"))
  expect_error(locus_class$new(10, -3))
})


test_that("Adding a locus to a model works", {
  model2 <- coal_model(10) + locus_averaged(10, 1010)
  expect_equal(get_locus_number(model2), 10)
  expect_equal(get_locus_length(model2, 1), 1010)

  model2 <- model2 +
    locus_single(1017) +
    locus_single(1018) +
    locus_single(1019)

  expect_equal(get_locus_number(model2), 13)
  expect_equal(get_locus_length(model2, 1), 1010)
  expect_equal(get_locus_length(model2, 13), 1019)
})


test_that("Adding a locus trio works", {
  m <- coal_model(10) +
    locus_trio(locus_length = c(10, 30, 50),
               distance = c(20, 40)) +
    locus_trio(locus_length = c(middle = 30, left = 10, right = 50),
               distance = c(20, 40)) +
    locus_trio(locus_length = c(10, 30, 50),
               distance = c(middle_right = 40, left_middle = 20))

  expect_equivalent(get_locus_length_matrix(m),
                    cbind(matrix(1:5 * 10, 3, 5, byrow = TRUE), 1))
})


test_that("locus positions are converted correctly", {
  model <- coal_model(5:6) +
    locus_trio(locus_length = c(10, 30, 50), distance = c(20, 40)) +
    locus_trio(locus_length = c(50, 30, 10), distance = c(40, 20)) +
    locus_averaged(2, 100)

  expect_equal(conv_middle_to_trio_pos(.5, model, 1), 45 / 150)
  expect_equal(conv_middle_to_trio_pos(.5, model, 2), 105 / 150)
  expect_equal(conv_middle_to_trio_pos(.5, model, 3), .5)
  expect_equal(conv_middle_to_trio_pos(.5, model, 4), .5)

  expect_equal(conv_middle_to_trio_pos(15, model, 1, relative_in = FALSE),
               45 / 150)
  expect_equal(conv_middle_to_trio_pos(15, model, 2, relative_in = FALSE),
               105 / 150)
  expect_equal(conv_middle_to_trio_pos(15, model, 3, relative_in = FALSE), .15)

  expect_equal(conv_middle_to_trio_pos(.5, model, 1, relative_out = FALSE), 45)
  expect_equal(conv_middle_to_trio_pos(.5, model, 2, relative_out = FALSE), 105)
  expect_equal(conv_middle_to_trio_pos(.5, model, 3, relative_out = FALSE), 50)

  expect_equal(conv_middle_to_trio_pos(10, model, 1,
                                       relative_out = FALSE,
                                       relative_in = FALSE), 40)

  ss <- create_segsites(matrix(0, 6, 6),
                        c(0.1, 0.5, 0.2, 0.6, 0.5, 1),
                        rep(c(-1, 0, 1), each = 2),
                        check = FALSE)
  expect_equal(get_snp_positions(list(ss, ss, ss, ss), model),
               list(c(1, 5, 36, 48, 125, 150) / 150,
                    c(5, 25, 96, 108, 145, 150) / 150,
                    c(0.1, 0.5, 0.2, 0.6, 0.5, 1),
                    c(0.1, 0.5, 0.2, 0.6, 0.5, 1)))
  expect_equal(get_snp_positions(list(ss, ss, ss, ss),
                                 model, relative = FALSE),
               list(c(1, 5, 36, 48, 125, 150),
                    c(5, 25, 96, 108, 145, 150),
                    c(10, 50, 20, 60, 50, 100),
                    c(10, 50, 20, 60, 50, 100)))
})
