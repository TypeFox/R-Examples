context("dc_search")


test_that("dc_search basic functionality works", {
  skip_on_cran()

  # default uses q=*:*
  aa <- dc_search(verbose = FALSE)
  # basic search
  bb <- dc_search(q = "laser", rows = 2, verbose = FALSE)
  # specify fields to get back
  bb <- dc_search(q = "laser", fl=c('doi','publicationYear'), rows=2, verbose = FALSE)
  # search a specific field
  cc <- dc_search(q = "subject:geology", fl=c('doi','subject'), rows=2, verbose = FALSE)

  expect_is(aa, "tbl_df")
  expect_is(aa$minted, "character")
  expect_is(aa$language, "character")
  expect_less_than(20, NCOL(aa))

  expect_is(bb, "tbl_df")
  expect_named(bb, c('doi', 'publicationYear'))
  expect_equal(NROW(bb), 2)

  expect_is(cc, "tbl_df")
  expect_named(cc, c('doi', 'subject'))
  expect_equal(NROW(cc), 2)
})

test_that("dc_search works w/ csv output", {
  skip_on_cran()

  aa <- dc_search(q = 'wind', fl=c('doi','title'), wt='csv', verbose = FALSE)
  bb <- dc_search(q = 'wind', fl=c('doi','title'), verbose = FALSE)

  expect_is(aa, "tbl_df")
  expect_named(aa, c('doi', 'title'))

  # csv and json output differ, some encoding problem likely
  expect_false(identical(aa, bb))
})

test_that("dc_search fails nicely", {
  skip_on_cran()
  library('httr')

  expect_error(dc_search(callopts = timeout(0.01)), "Timeout was reached")
})
