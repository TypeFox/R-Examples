# Verify that image magics can successfully be identified 
# 
# Author: jefferis
###############################################################################
context("Verify that image magics can be identified")

test_that("image_type identifies bmps", {
      expect_that(image_type("testdata/real.bmp"),
          equals("bmp"))
      expect_that(image_type("testdata/bmp-pretending-to-be.jpg"),
          equals('bmp'))
    })

test_that("image_type returns NA_character_ for pnms", {
      expect_that(image_type("testdata/real.pnm"),
          equals(NA_character_))
    })

test_that("image_type identifies pngs", {
      expect_that(image_type("testdata/real.png"),
          equals('png'))
      expect_that(image_type("testdata/png-pretending-to-be.bmp"),
          equals('png'))
    })

test_that("image_type identifies jpgs", {
      expect_that(image_type("testdata/real.jpg"),
          equals('jpg'))
    })
