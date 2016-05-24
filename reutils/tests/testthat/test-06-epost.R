
# Test epost() -----------------------------------------------------------

context("Testing 'epost()'")

if (getOption('reutils.test.remote')) {
  uids <- esearch("Chlamydia", "bioproject", retmax = 20)
  
  test_that("epost() works with 'esearch' objects", {
    p1 <- epost(uids)
    expect_is(p1, 'epost')
    expect_match(webenv(p1), "NCID_.+")
    expect_equal(querykey(p1), 1)
  })
  
  test_that("epost() works with 'entrez_uid' objects", {
    p2 <- epost(content(uids, 'parsed'))
    expect_is(p2, 'epost')
    expect_match(webenv(p2), "NCID_.+")
    expect_equal(querykey(p2), 1)
  })
  
  test_that("epost() works with character/numeric vectors", {
    p3 <- epost(c("194680922", "50978626", "28558982", "9507199", "6678417"), "protein")
    expect_match(webenv(p3), "NCID_.+")
    expect_equal(querykey(p3), 1)
  })
}
