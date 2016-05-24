# tests for bold_seq fxn in bold
context("bold_seq")

test_that("bold_seq returns the correct dimensions/classes", {
  skip_on_cran()
  
  a <- bold_seq(taxon='Coelioxys')
  b <- bold_seq(bin='BOLD:AAA5125')
  c <- bold_seq(taxon='Coelioxys', response=TRUE)
  
  expect_equal(c$status_code, 200)
  expect_equal(c$headers$`content-type`, "application/x-download")

  expect_is(a, "list")
  expect_is(b, "list")

  expect_is(a[[1]], "list")
  expect_is(a[[1]]$id, "character")
  expect_is(a[[1]]$sequence, "character")

  expect_is(c, "response")
  expect_is(c$headers, "insensitive")
})

test_that("bold_seq returns correct error when parameters empty or not given", {
  skip_on_cran()
  
  expect_error(bold_seq(taxon=''), "must provide a non-empty value")
  expect_error(bold_seq(), "must provide a non-empty value")
})
