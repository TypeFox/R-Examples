# tests for almevents fxn in alm
context("almevents")

out <- alm_events(doi="10.1371/journal.pone.0029797")
out <- out[!out %in% c("sorry, no events content yet","parser not written yet")] # remove those with no data

test_that("almevents returns the correct class", {
	expect_that(out, is_a("list"))
	expect_that(out[["pmc"]], is_a("list"))
})

test_that("almevents returns correct things when two dois passed in", {
  dois <- c('10.1371/journal.pone.0001543','10.1371/journal.pone.0040117')
  out2 <- alm_events(doi=dois)
  expect_is(out2, "list")
  expect_equal(length(out2), 2)
  expect_equal(length(out2), 2)
})

test_that("almevents returns correctly when one specific source requested", {
  out3 <- alm_events(source_id="reddit")
  expect_equal(length(out3), 50)
  expect_match(names(out3[[1]]), "reddit")
})
