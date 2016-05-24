context("Input checks")

test_that("Prompts errors correctly", {
  expect_error(check_type(doc = "http://http://cran.at.r-project.org/", which = factor(4)))
})

test_that("check_type produces class output", {

  x <- check_type(doc = "http://christianrubba.com/htmltab/ex/wiki_indian_election2014.html", which = "//table[5]")
  expect_that(x, is_a("XMLInternalDocument"))

  expect_error(check_type(doc = "http://christianrubba.com/htmltab/ex/wiki_indian_election2014.html", which = 1))

  parsed <- XML::htmlParse("http://christianrubba.com/htmltab/ex/wiki_indian_election2014.html")
  z <- check_type(doc = parsed, which = 3)
  expect_that(z, is_a("XMLInternalDocument"))
})
