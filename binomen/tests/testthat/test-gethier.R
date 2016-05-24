context("gethier")

test_that("gethier works with grouping class input - eg1", {
  out <- grouping(kingdom = taxonref("kingdom", "Animalia"),
                species = taxonref("species", "Homo sapiens"))
  aa <- gethier(out)

  expect_is(out, "grouping")
  expect_is(aa, "data.frame")
  expect_is(aa$rank, "character")
  expect_equal(NROW(aa), 2)
})

test_that("gethier works with grouping class input - eg2", {
  out <- grouping(order = taxonref("order", "hymenoptera"),
                  genus = taxonref("genus", "Bembidion"),
                  species = taxonref("species", "Bembidion punctatus"))
  aa <- gethier(out)

  expect_is(out, "grouping")
  expect_is(aa, "data.frame")
  expect_is(aa$rank, "character")
  expect_equal(NROW(aa), 3)
})

test_that("gethier works with taxon class input - eg1", {
  out <- taxon(binomial(genus = "Poa", epithet = "annua", authority = "L."),
              grouping(family = taxonref("family", 'Poaceae'),
                       genus = taxonref("genus", 'Poa')))
  aa <- gethier(out)

  expect_is(out, "taxon")
  expect_is(out$binomial, "binomial")
  expect_is(out$grouping, "grouping")
  expect_is(aa, "data.frame")
  expect_is(aa$rank, "character")
  expect_equal(NROW(aa), 2)
})

test_that("gethier check for class works", {
  expect_error(gethier(), "no applicable method")
  expect_error(gethier.grouping(), "argument \"x\" is missing, with no default")
  expect_error(gethier.taxon(), "argument \"x\" is missing, with no default")
})
