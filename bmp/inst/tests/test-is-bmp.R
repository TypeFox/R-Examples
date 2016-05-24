context("Check that is.bmp behaves")

require(pixmap)

test_that("is.bmp returns TRUE for bmps", {
      expect_that(is.bmp("../images/5HT1bMARCM-F000001_seg001_lsm.bmp"),
          equals(TRUE))
      expect_that(is.bmp("../images/w5h3-8bit.bmp"),
          equals(TRUE))
    })

test_that("is.bmp returns FALSE for pnms", {
      expect_that(is.bmp("../images/5HT1bMARCM-F000001_seg001_lsm.pnm"),
          equals(FALSE))
      expect_that(is.bmp("../images/w5h3-8bit.pgm"),
          equals(FALSE))
    })
