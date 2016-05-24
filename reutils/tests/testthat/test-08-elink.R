
# Test elink() -----------------------------------------------------------

context("Testing 'elink()'")

if (getOption('reutils.test.remote')) {
  x <- elink(c("24475906", "34577062"), dbFrom="nuccore")
  
  test_that("elink() returns an 'elink' object", {
    expect_is(x, "elink")
  })
  
  test_that("'content()' returns an entrez_linkset or an XMLInternalDocument", {
    expect_that(content(x, "text"), is_a("character"))
    expect_that(content(x, "xml"), is_a("XMLInternalDocument"))
    expect_that(content(x, 'parsed'), is_a("entrez_linkset"))
  })
  
  test_that("'uid', 'database', and 'linkset' return the appropriate results", {
    expect_equal(uid(x), c("24475906", "34577062"))
    expect_equal(database(x), "nuccore")
    expect_that(linkset(x), is_a("list"))
  })
  
}
