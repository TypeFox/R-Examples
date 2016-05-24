library(wkb)
context("Conversion from string hex representation to raw vector")

# create a raw vector
refwkbvector <- as.raw(c(0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                         0xf0, 0x3f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08, 0x40))

# create a list of two raw vectors
refwkblist <- list(
  as.raw(c(0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
           0xf0, 0x3f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08, 0x40)),
  as.raw(c(0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
           0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40))
)

test_that("character string hexadecimal representation converts to raw vector", {
  # create a character string containing a hexadecimal representation
  hex <- "0101000000000000000000f03f0000000000000840"

  # convert to raw vector
  wkb <- hex2raw(hex)

  # test
  expect_identical(wkb, refwkbvector)
})

test_that("character vector hexadecimal representation converts to raw vector", {
  # create a character vector containing a hexadecimal representation
  hex <- c("01", "01", "00", "00", "00", "00", "00", "00", "00", "00", "00",
           "f0", "3f", "00", "00", "00", "00", "00", "00", "08", "40")

  # convert to raw vector
  wkb <- hex2raw(hex)

  # test
  expect_identical(wkb, refwkbvector)
})

test_that("vector of character string hexadecimal representations converts to list of raw vectors", {
  # create vector of two character strings each containing a hex representation
  hex <- c("0101000000000000000000f03f0000000000000840",
           "010100000000000000000000400000000000000040")

  # convert to list of two raw vectors
  wkb <- hex2raw(hex)

  # test
  expect_identical(wkb, refwkblist)
})
