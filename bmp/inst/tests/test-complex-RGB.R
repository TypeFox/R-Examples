context("Reading complex RGB images")

require(pixmap)

test_that("complex RGB BMP image as pixmap is identical to pnm loaded by read.pnm", {
      pnmfile="../images/5HT1bMARCM-F000001_seg001_lsm.pnm"
      bmpfile=sub("pnm$","bmp",pnmfile)
      # test red green and blue channel separately
      b = read.bmp(bmpfile)
      pb=pixmapRGB(b)
      p = read.pnm(pnmfile)
      expect_that(pb,equals(p))
    })
