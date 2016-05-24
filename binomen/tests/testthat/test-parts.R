context("parts")

out <- make_taxon(genus = "Poa", epithet = "annua", authority = "L.",
   family = 'Poaceae', clazz = 'Poales', kingdom = 'Plantae',
   variety = 'annua')

test_that("name w/ taxon objects", {
  aa <- name(out)

  expect_is(aa, "character")
  expect_is(aa[1], "character")
  expect_equal(aa[1], "Plantae")
  expect_equal(length(aa), 6)
})

test_that("name w/ taxon objects - piping", {
  aa <- out %>% name()

  expect_is(aa, "character")
  expect_is(aa[1], "character")
  expect_equal(aa[1], "Plantae")
  expect_equal(length(aa), 6)
})

test_that("uri w/ taxon objects", {
  aa <- uri(out)

  expect_is(aa, "character")
  expect_is(aa[1], "character")
  expect_equal(aa[1], "none")
  expect_equal(length(aa), 6)
})

test_that("uri w/ taxon objects - piping", {
  aa <- out %>% uri()

  expect_is(aa, "character")
  expect_is(aa[1], "character")
  expect_equal(aa[1], "none")
  expect_equal(length(aa), 6)
})

test_that("rank w/ taxon objects", {
  aa <- rank(out)

  expect_is(aa, "character")
  expect_is(aa[1], "character")
  expect_equal(aa[1], "kingdom")
  expect_equal(length(aa), 6)
})

test_that("rank w/ taxon objects - piping", {
  aa <- out %>% rank()

  expect_is(aa, "character")
  expect_is(aa[1], "character")
  expect_equal(aa[1], "kingdom")
  expect_equal(length(aa), 6)
})

test_that("taxonid w/ taxon objects", {
  aa <- taxonid(out)

  expect_is(aa, "character")
  expect_is(aa[1], "character")
  expect_equal(aa[1], "none")
  expect_equal(length(aa), 6)
})

test_that("taxonid w/ taxon objects - piping", {
  aa <- out %>% taxonid()

  expect_is(aa, "character")
  expect_is(aa[1], "character")
  expect_equal(aa[1], "none")
  expect_equal(length(aa), 6)
})

test_that("name fails well", {
  expect_error(name(), "no applicable method")
  expect_error(name("ad"), "no applicable method")
})

test_that("uri fails well", {
  expect_error(uri(), "no applicable method")
  expect_error(uri("ad"), "no applicable method")
})

test_that("rank fails well", {
  expect_error(rank(), "no applicable method")
  expect_error(rank("ad"), "no applicable method")
})

test_that("taxonid fails well", {
  expect_error(taxonid(), "no applicable method")
  expect_error(taxonid("ad"), "no applicable method")
})
