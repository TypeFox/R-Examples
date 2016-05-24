context("Reading simple greyscale images")

require(pixmap)

test_that("8 bit bmp image matches pgm equivalent loaded with pixmap", {
      pgmfile="../images/w5h3-8bit.pgm"
      bmpfile=sub("pgm$","bmp",pgmfile)
      expect_that(read.bmp(bmpfile)/(2^8-1),
          is_equivalent_to(read.pnm(pgmfile)@grey))
    })

test_that("Trying to load a non BMP file throws an error", {
      pgmfile="../images/w5h3-8bit.pgm"
      expect_that(read.bmp(pgmfile),
          throws_error())
    })

# don't have an example of this
#test_that("16 bit bmp image matches pgm equivalent loaded with pixmap", {
#      pgmfile="../images/w5h3-16bit.pgm"
#      bmpfile=sub("pgm$","bmp",pgmfile)
#      expect_that(read.bmp(bmpfile)/(2^16-1),
#          is_equivalent_to(read.pnm(pgmfile)@grey))
#    })
