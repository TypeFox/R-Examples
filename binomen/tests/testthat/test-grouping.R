context("grouping")

test_that("grouping basic functionality works", {
  aa <- grouping(kingdom=taxonref("kingdom", "Animalia"),
                       species=taxonref("species", "Homo sapiens"))

  expect_is(aa, "grouping")
  expect_is(aa$kingdom, "taxonref")
  expect_is(aa$species, "taxonref")

  expect_named(aa, c("kingdom", "species"))
  expect_named(aa$kingdom, c("rank", "name", "id", "uri"))

  expect_equal(aa$kingdom$rank, "kingdom")
  expect_equal(aa$species$name, "Homo sapiens")
  expect_equal(length(aa), 2)
  expect_equal(length(aa$kingdom), 4)
})

test_that("grouping fails well", {
  expect_equal(length(grouping()), 0)
  expect_error(grouping(stuff = 5), "unused argument")
  expect_error(grouping(kingdom = "stuff"),
               "One or more inputs was not of class taxonref")
  expect_error(grouping(division = "stuff"),
               "One or more inputs was not of class taxonref")
})
