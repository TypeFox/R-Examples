# tests for alm_signposts fxn in alm
context("alm_signposts")

dat1 <- alm_signposts(doi="10.1371/journal.pone.0029797")
dois <- c('10.1371/journal.pone.0001543','10.1371/journal.pone.0040117',
 '10.1371/journal.pone.0029797','10.1371/journal.pone.0039395')
dat2 <- alm_signposts(doi=dois)

test_that("alm_signposts returns the correct class", {
  expect_that(dat1, is_a("data.frame"))
  expect_that(dat2, is_a("data.frame"))
  expect_that(dat1$viewed, is_a("integer"))
  expect_that(dat2$cited, is_a("integer"))
  expect_that(dat2$id, is_a("character"))
})

test_that("alm_signposts returns the correct dimensions", {
  expect_that(nrow(dat1), equals(1))
  expect_that(nrow(dat2), equals(4))
})
