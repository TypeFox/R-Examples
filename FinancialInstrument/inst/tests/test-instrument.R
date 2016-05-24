context("instruments")

test_that("stock creates stock", {
  # stock is created and assigned in .instrument environment and the primary_id
  # of the stock is returned.
  expect_is(stock("AAA", currency("USD")), "character")
  expect_is(getInstrument("AAA", type="stock"), "stock")
  expect_true(is.null(names(stock(c("AAA", "BBB"), "USD"))))
})

test_that("stock not assigned", {
  # overwrite=FALSE is ignored because assign_i=FALSE
  s <- c("BBB", "AAA")
  ilist <- stock(s, "USD", assign_i=FALSE, overwrite=FALSE)
  expect_is(ilist, "list")
  expect_true(all(vapply(ilist, inherits, FUN.VALUE=TRUE, "stock")))
  expect_identical(names(ilist), s)
})

test_that("stock overwrite throws errors", {
  expect_error(stock("AAA", "USD", overwrite=FALSE))
  rm_stocks("BBB")
  expect_error(stock(c("BBB", "AAA"), "USD", overwrite=FALSE))
  # Make sure it didn't define BBB
  expect_true(!getInstrument("BBB", type="stock", silent=TRUE))
})

test_that("loadInstruments from list", {
  rm_instruments(keep.currencies=FALSE)
  stock(c("A", "B"), currency("USD")) # put some stuff in the .instrument env
  L <- stock(c("DD", "EE"), "USD", assign_i=FALSE)
  expect_true(!is.instrument.name("DD"))
  loadInstruments(L)
  expect_true(is.instrument.name("DD"))
  expect_true(is.instrument.name("A"))
  reloadInstruments(L)
  expect_true(!is.instrument.name("A"))
  expect_true(is.instrument.name("DD"))
})