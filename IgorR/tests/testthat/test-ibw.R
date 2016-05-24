# Tests for reading of Igor Binary Wave (.ibw) files
# 
# Author: jefferis
###############################################################################

context("Verify reading of Igor ibw files")

test_that("Read Igor v5 file", {
      expect_that (w5<-read.ibw(system.file('igor/version5.ibw', package = 'IgorR')),
          is_a("numeric"),'read igor v5 wave without error')
      expect_that(attr(w5,'BinHeader')$version,
          equals(5))
      expect_that(tsp.igorwave(w5),
          equals(structure(c(0, 4, 1), .Names = c("start", "end", "frequency"))))
      
    })

test_that("Read Igor v2 file", {
      expect_that (w2<-read.ibw(system.file('igor/version2.ibw', package = 'IgorR')),
          is_a("numeric"),'read igor v2 wave without error')
      expect_that(attr(w2,'BinHeader')$version,
          equals(2))
      expect_that(tsp.igorwave(w2),
          equals(structure(c(0, 4, 1), .Names = c("start", "end", "frequency"))))
    })

test_that("Processing of dates", {
  expect_equivalent(.convertIgorDate(3444214517),as.POSIXct("2013-02-20 14:15:17"))
  # from nm20120811c1_016.pxp 
  # current implementation says 13:37 on my laptop, lmb cluster and mac in PDT
  # but mtime of file and embedded time string is 12:37 
  # and this is also what I get on a sunos machine in PDT
  # DISABLED due to OS dependent results that I cannot fix
  # expect_equivalent(.convertIgorDate(3427533421),as.POSIXct("2012-08-11 13:37:01"))
})

test_that("Processing of funny null terminated strings", {
  expect_is(bad_string_wave<-read.ibw('igor/bad_null_string.ibw'), 'numeric')
})
