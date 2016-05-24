context("dump")

test_that("dump_links works", {
  skip_on_cran()
  
  aa <- dump_links()
  expect_is(aa, "list")
  expect_is(aa$mammals, "list")
  expect_is(aa$reptiles, "list")
  expect_is(aa$amphibians, "list")
  expect_is(aa$fishes, "list")
  expect_is(aa$birds, "list")
  expect_named(aa, c('mammals','reptiles','amphibians','fishes','birds'))
  expect_true(grepl("https", aa$mammals$view))
  expect_true(grepl("urn:uuid", aa$fishes$data))
})
