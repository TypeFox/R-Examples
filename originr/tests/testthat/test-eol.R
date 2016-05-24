context("eol functions")

test_that("eol works", {
  skip_on_cran()

  aa <- eol(name = 'Brassica oleracea', dataset = 'gisd', verbose = FALSE)
  bb <- eol(name = 'Ciona intestinalis', dataset = 'mineps', verbose = FALSE)
  cc <- eol(name = c('Lymantria dispar','Cygnus olor'), dataset = 'i3n', verbose = FALSE)

  expect_is(aa, "data.frame")
  expect_equal(aa$searched_name, 'Brassica oleracea')
  expect_is(aa$searched_name, 'character')
  expect_equal(aa$db, "gisd")

  expect_is(bb, "data.frame")
  expect_equal(bb$searched_name, 'Ciona intestinalis')
  expect_is(bb$searched_name, 'character')
  expect_equal(bb$db, "mineps")

  expect_is(cc, "data.frame")
  expect_equal(cc$searched_name, c('Lymantria dispar','Cygnus olor'))
  expect_is(cc$searched_name, 'character')
  expect_equal(cc$db[1], "i3n")
})

test_that("fails well", {
  skip_on_cran()

  expect_error(eol(), "please provide a taxonomic name")
})
