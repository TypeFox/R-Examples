
# Test esummary() ---------------------------------------------------------

context("Testing 'esummary()'")

if (getOption('reutils.test.remote')) {
  x <- esearch(term = "Chlamydia psittaci", db = "nuccore", retmax = 2)
  a <- esummary(x)
  
  test_that("esummary() returns an 'esummary' object", {
    expect_is(a, "esummary")
  })
  
  test_that("'retmode' returns 'xml'", {
    expect_that(retmode(a), equals("xml"))
  })
  
  test_that("'content()' returns an XMLInternalDocument or data.fame", {
    expect_that(content(a), is_a("XMLInternalDocument"))
    expect_that(content(a, 'json'), throws_error("Cannot return data of retmode.+"))
    expect_that(content(a, 'parsed'), is_a("data.frame"))
  })
}
