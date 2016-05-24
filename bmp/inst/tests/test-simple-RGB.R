context("Reading simple RGB images")

require(pixmap)

test_that("8 bit RGB BMP image matches pnm equivalent loaded with pixmap", {
      pnmfile="../images/w5h4-RGB.pnm"
      bmpfile=sub("pnm$","bmp",pnmfile)
      # test red green and blue channel separately
      b = read.bmp(bmpfile)/(2^8-1)
      p = read.pnm(pnmfile)
      expect_that(b[,,1],equals(p@red))
      expect_that(b[,,2],equals(p@green))
      expect_that(b[,,3],equals(p@blue))
    })

test_that("32 bit ARGB BMP image matches pnm equivalent loaded with pixmap", {
      pnmfile="../images/explosion_32bit.pnm"
      bmpfile=sub("pnm$","bmp",pnmfile)
      # test red green and blue channel separately
      b = read.bmp(bmpfile)/(2^8-1)
      p = read.pnm(pnmfile)
      expect_that(b[,,2],equals(p@red))
      expect_that(b[,,3],equals(p@green))
      expect_that(b[,,4],equals(p@blue))
    })

test_that("RGB BMP image converted to pixmap is identical to pnm loaded by read.pnm", {
      pnmfile="../images/w5h4-RGB.pnm"
      bmpfile=sub("pnm$","bmp",pnmfile)
      # test red green and blue channel separately
      b = read.bmp(bmpfile)
      pb=pixmapRGB(b)
      p = read.pnm(pnmfile)
      expect_that(pb,equals(p))
    })
