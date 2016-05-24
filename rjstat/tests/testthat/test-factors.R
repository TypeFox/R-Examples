context("Factors")

dataset <- readLines("dataset.json")

test_that("factors are factors", {
    expect_is(fromJSONstat(dataset, naming = "label",
                           use_factors = TRUE)[[1]][[1]],
              "factor")
    expect_is(fromJSONstat(dataset, naming = "id",
                           use_factors = TRUE)[[1]][[1]],
              "factor")

    expect_is(fromJSONstat(dataset, naming = "label",
                           use_factors = TRUE)[[1]][[2]],
              "factor")
    expect_is(fromJSONstat(dataset, naming = "id",
                           use_factors = TRUE)[[1]][[2]],
              "factor")

    expect_is(fromJSONstat(dataset, naming = "label",
                           use_factors = TRUE)[[1]][[3]],
              "factor")
    expect_is(fromJSONstat(dataset, naming = "id",
                           use_factors = TRUE)[[1]][[3]],
              "factor")
})

test_that("factor levels are correct", {
    expect_equal(nlevels(fromJSONstat(dataset, naming = "label",
                                      use_factors = TRUE)[[1]][[1]]),
                 2)
    expect_equal(levels(fromJSONstat(dataset, naming = "label",
                                     use_factors = TRUE)[[1]][[1]]),
                 c("Category 11", "Category 12"))
    expect_equal(nlevels(fromJSONstat(dataset, naming = "id",
                                      use_factors = TRUE)[[1]][[1]]),
                 2)
    expect_equal(levels(fromJSONstat(dataset, naming = "id",
                                     use_factors = TRUE)[[1]][[1]]),
                 c("testcategory11", "testcategory12"))

    expect_equal(nlevels(fromJSONstat(dataset, naming = "label",
                                      use_factors = TRUE)[[1]][[2]]),
                 2)
    expect_equal(levels(fromJSONstat(dataset, naming = "label",
                                     use_factors = TRUE)[[1]][[2]]),
                 c("Category 21", "Category 22"))
    expect_equal(nlevels(fromJSONstat(dataset, naming = "id",
                                      use_factors = TRUE)[[1]][[2]]),
                 2)
    expect_equal(levels(fromJSONstat(dataset, naming = "id",
                                     use_factors = TRUE)[[1]][[2]]),
                 c("testcategory21", "testcategory22"))

    expect_equal(nlevels(fromJSONstat(dataset, naming = "label",
                                      use_factors = TRUE)[[1]][[3]]),
                 1)
    expect_equal(levels(fromJSONstat(dataset, naming = "label",
                                     use_factors = TRUE)[[1]][[3]]),
                 "Category 3")
    expect_equal(nlevels(fromJSONstat(dataset, naming = "id",
                                      use_factors = TRUE)[[1]][[3]]),
                 1)
    expect_equal(levels(fromJSONstat(dataset, naming = "id",
                                     use_factors = TRUE)[[1]][[3]]),
                 "testcategory3")
})

test_that("factor integer codes are correct", {
    expect_equivalent(unclass(fromJSONstat(dataset,
                                           use_factors = TRUE)[[1]][[1]]),
                      c(1, 1, 2, 2))
    expect_equivalent(unclass(fromJSONstat(dataset,
                                           use_factors = TRUE)[[1]][[2]]),
                      c(1, 2, 1, 2))
    expect_equivalent(unclass(fromJSONstat(dataset,
                                           use_factors = TRUE)[[1]][[3]]),
                      c(1, 1, 1, 1))
})
