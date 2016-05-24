#' # Test utf-8
#' This test is mainly interesting for windows
require(testthat)
require(PythonInR)
invisible(capture.output(pyConnect()))

rAscii <- "abcABC"
pySet("ascii", rAscii)
expect_that(pyGet("ascii"), equals(rAscii))

# is necessary for more info see 005_execute_file.R
rUtf8 <- iconv("äöüßÄöü", from="UTF-8")
    
pySet("utf8", rUtf8)
expect_that(pyGet("utf8"), equals(rUtf8))

