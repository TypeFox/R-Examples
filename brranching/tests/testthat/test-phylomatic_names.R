context("phylomatic_names")

mynames <- c("Poa annua", "Salix goodingii", "Helianthus annuus")

test_that("phylomatic_names - rsubmit works", {
  skip_on_cran()

  nms <- suppressMessages(phylomatic_names(mynames, format = 'rsubmit'))

  expect_is(nms, "character")
  expect_true(grepl("\\%", nms[1]))
  expect_equal(strsplit(nms[1], "\\%")[[1]][1], "poaceae")
})

test_that("phylomatic_names - isubmit (default) works", {
  skip_on_cran()

  nms <- suppressMessages(phylomatic_names(taxa = mynames, format = 'isubmit'))

  expect_is(nms, "character")
  expect_false(grepl("\\%", nms[1]))
  expect_true(grepl("\\/", nms[1]))
  expect_equal(strsplit(nms[1], "\\/")[[1]][1], "poaceae")
})

test_that("phylomatic_names - db=apg works", {
  skip_on_cran()

  nms <- suppressMessages(phylomatic_names(mynames, db = "apg"))

  expect_is(nms, "character")
  expect_false(grepl("\\%", nms[1]))
  expect_true(grepl("\\/", nms[1]))
  expect_equal(strsplit(nms[1], "\\/")[[1]][1], "poaceae")
})

test_that("phylomatic fails as expected", {
  skip_on_cran()

  # fails when no taxa given
  expect_error(phylomatic_names(), "argument \"taxa\" is missing")

  # fails when format isn't of allowed values
  expect_error(phylomatic_names(mynames, format = "asdfadf"),
               "'arg' should be one of")

  # fails when db isn't of allowed values
  expect_error(phylomatic_names(mynames, db = "things"),
               "'arg' should be one of")
})
