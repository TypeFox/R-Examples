library(wkb)
library(sp)
context("Conversion to and from WKB Point representations")

# create an object of class SpatialPoints
x = c(1, 2)
y = c(3, 2)
refobj <- SpatialPoints(data.frame(x, y))

# create a AsIs list of little-endian WKB geometry representations of type Point
refwkb <- I(list(
  as.raw(c(0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
           0xf0, 0x3f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08, 0x40)),
  as.raw(c(0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
           0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40))
))

# create an AsIs list of big-endian WKB geometry representations of type Point
refwkbbe <- I(list(
  as.raw(c(0x00, 0x00, 0x00, 0x00, 0x01, 0x3f, 0xf0, 0x00, 0x00, 0x00, 0x00,
           0x00, 0x00, 0x40, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00)),
  as.raw(c(0x00, 0x00, 0x00, 0x00, 0x01, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00,
           0x00, 0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00))
))

test_that("little-endian WKB Point representation converts to SpatialPoints object", {
  # convert little-endian WKB Point representation to SpatialPoints object
  obj <- readWKB(refwkb)

  # test
  expect_equal(obj, refobj)
})

test_that("SpatialPoints object converts to little-endian WKB Point representation", {
  # convert SpatialPoints object to little-endian WKB Point representation
  wkb <- writeWKB(refobj)

  # test
  expect_equal(wkb, refwkb)
})

test_that("big-endian WKB Point representation converts to SpatialPoints object", {
  # convert big-endian WKB Point representation to SpatialPoints object
  obj <- readWKB(refwkbbe)

  # test
  expect_equal(obj, refobj)
})

test_that("SpatialPoints object converts to big-endian WKB Point representation", {
  # convert SpatialPoints object to big-endian WKB Point representation
  wkbbe <- writeWKB(refobj, endian = "big")

  # test
  expect_equal(wkbbe, refwkbbe)
})
