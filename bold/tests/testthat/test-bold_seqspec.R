# tests for bold_seqspec fxn in bold
context("bold_seqspec")

test_that("bold_seqspec returns the correct dimensions or values", {
  skip_on_cran()

  a <- bold_seqspec(taxon='Osmia')
  b <- bold_seqspec(taxon='Osmia', response=TRUE)
  c <- bold_seqspec(taxon='Osmia', sepfasta=TRUE)

  expect_equal(b$status_code, 200)
  expect_equal(b$headers$`content-type`, "application/x-download")

  expect_is(a, "data.frame")
  expect_is(b, "response")
  expect_is(c, "list")
  expect_is(c$data, "data.frame")
  expect_is(c$fasta, "list")
  expect_is(c$fasta[[1]], "character")

  expect_is(a$recordID, "integer")
  expect_is(a$directions, "character")

  expect_is(b$headers, "insensitive")
})

test_that("bold_seq returns correct error when parameters empty or not given", {
  skip_on_cran()

  expect_error(bold_seqspec(taxon=''), "must provide a non-empty value")
  expect_error(bold_seqspec(), "must provide a non-empty value")
})
