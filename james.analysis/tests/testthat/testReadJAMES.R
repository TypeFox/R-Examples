context("Read JSON files")

test_that("returned data has class 'james'", {
  # get path to raw JSON file included in package distribution
  json.file <- system.file("extdata", "james.json", package = "james.analysis") 
  # read results from file
  james <- readJAMES(json.file)
  # verify
  expect_is(james, "james")
})

test_that("readJAMES complains on invalid arguments", {
  expect_error(readJAMES(NULL), "non-empty string")
  expect_error(readJAMES(NA), "non-empty string")
  expect_error(readJAMES(""), "non-empty string")
  expect_error(readJAMES(c("a", "b")), "non-empty string")
})

test_that("readJAMES complains if the input file is not found", {
  expect_error(readJAMES("i-do-not-exist"), "does not exist")
})

test_that("readJAMES checks JSON structure", {
  # write temporary JSON file with bogus structure
  file <- tempfile()
  content <- '{"employees":[{"firstName":"John","lastName":"Doe"},{"firstName":"Anna","lastName":"Smith"},{"firstName":"Peter","lastName":"Jones"}]}'
  writeLines(content, file)
  # verify
  expect_error(readJAMES(file), "unexpected JSON structure")
  # overwrite with empty array
  writeLines("[]", file)
  # verify
  expect_error(readJAMES(file), "unexpected JSON structure")
  # remove temporary file
  file.remove(file)
})
