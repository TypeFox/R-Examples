context("strain")

out <- make_taxon(genus="Poa", epithet="annua", authority="L.",
   family='Poaceae', clazz='Poales', kingdom='Plantae', variety='annua')

test_that("less than works", {
  aa <- out %>% strain(. < family)

  expect_is(aa, "taxon")
  expect_is(aa$binomial, "binomial")
  expect_is(aa$binomial$authority, "character")
  expect_is(aa$grouping, "grouping")
  expect_equal(names(aa$grouping), c('genus', 'species', 'variety'))
  expect_identical(aa$binomial, out$binomial)

  expect_equal(length(out$grouping), 6)
  expect_equal(length(aa$grouping), 3)
})

test_that("greater than works", {
  aa <- out %>% strain(. > family)

  expect_is(aa, "taxon")
  expect_is(aa$binomial, "binomial")
  expect_is(aa$binomial$authority, "character")
  expect_is(aa$grouping, "grouping")
  expect_equal(names(aa$grouping), c('kingdom', 'family'))
  expect_identical(aa$binomial, out$binomial)

  expect_equal(length(out$grouping), 6)
  expect_equal(length(aa$grouping), 2)
})
