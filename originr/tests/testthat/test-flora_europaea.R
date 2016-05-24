context("gisd functions")

test_that("gisd works", {
  skip_on_cran()

  sp <- c("Carpobrotus edulis", "Rosmarinus officinalis")
  aa <- flora_europaea(sp[1], verbose = FALSE)
  bb <- flora_europaea(sp[2], verbose = FALSE)

  expect_is(aa, "list")
  expect_named(aa, c('native', 'exotic', 'status_doubtful', 'occurrence_doubtful', 'extinct'))
  expect_equal(aa$native, NA_character_)
  expect_is(aa$exotic, "character")

  expect_is(bb, "list")
  expect_named(bb, c('native', 'exotic', 'status_doubtful', 'occurrence_doubtful', 'extinct'))
  expect_equal(bb$occurrence_doubtful, NA_character_)
  expect_is(bb$exotic, "character")
})

test_that("fails well - species not found when searching GBIF", {
  skip_on_cran()

  sp <- "asdfadsf"
  aa <- flora_europaea(sp, verbose = FALSE)
  expect_null(aa)
})
