context("Sumstat Tajima's D")

test_that("calulation works", {
  # example from http://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-genomics-fall-2005/study-materials/tajimad1.pdf #nolint
  ss <- create_segsites(matrix(c(1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0,
                                 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                                 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0,
                                 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                                 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
                                 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0), #nolint
                               10, 16, byrow = TRUE), 1:16 / 16)

  tajd <- sumstat_tajimas_d(population = 1)
  model <- coal_model(10, 1)
  expect_true(abs(tajd$calculate(list(ss), NULL, NULL, model) + 1.446172)
              < 1e-6)

  tajd <- sumstat_tajimas_d(population = "all")
  expect_true(abs(tajd$calculate(list(ss), NULL, NULL, model) + 1.446172)
              < 1e-6)

  expect_error(tajd$calculate(list(ss), NULL, coal_model(1, 1)))
  expect_equal(tajd$calculate(list(create_empty_segsites(10)),
                              NULL, NULL, model),
               NaN)
})
