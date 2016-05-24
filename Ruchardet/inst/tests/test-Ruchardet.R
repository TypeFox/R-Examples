context("test Rucharet")


test_that("test detectEncoding", {
          expect_true(file.exists(file.path(system.file("tests", package="Ruchardet"), "big5.txt")))
          expect_true(file.exists(file.path(system.file("tests", package="Ruchardet"), "gb18030.txt")))
          expect_true(file.exists(file.path(system.file("tests", package="Ruchardet"), "shift_jis.txt")))
          expect_true(file.exists(file.path(system.file("tests", package="Ruchardet"), "utf8.txt")))
          expect_true(file.exists(file.path(system.file("tests", package="Ruchardet"), "euc_kr.txt")))

          expect_equal(detectFileEncoding(file.path(system.file("tests", package="Ruchardet"), "big5.txt"),n=-1),      "Big5")
          expect_equal(detectFileEncoding(file.path(system.file("tests", package="Ruchardet"), "gb18030.txt"),n=-1),   "gb18030")
          expect_equal(detectFileEncoding(file.path(system.file("tests", package="Ruchardet"), "shift_jis.txt"),n=-1), "Shift_JIS")
          expect_equal(detectFileEncoding(file.path(system.file("tests", package="Ruchardet"), "utf8.txt"),n=-1),      "UTF-8")
          expect_equal(detectFileEncoding(file.path(system.file("tests", package="Ruchardet"), "euc_kr.txt"),n=-1),      "EUC-KR")
          
          fo <- file(file.path(system.file("tests", package="Ruchardet"), "shift_jis.txt"), 'r')
          expect_equal(detectFileEncoding(fo), "Shift_JIS")
          close(fo)

})



