# Testing code for logging functions

# `file.size()` is not defined on Windows R 3.1.1, which is CRAN's current `oldrel`
# (r-oldrel-windows-ix86+x86_64). Define a workaround for use below.
if(!exists("file.size"))
  file.size <- function(...) file.info(...)$size

context("logging")

test_that("functions handle bad input", {
  expect_error(openlog())
  expect_error(openlog("test", loglevel = TRUE))
  expect_error(openlog("test", loglevel = TRUE))
  expect_error(openlog("test", append = 1))
  expect_error(openlog("test", sink = 1))

  expect_error(printlog("test", level = TRUE))
  expect_error(printlog("test", ts = 1))
  expect_error(printlog("test", cr = 1))

  # printlog and closelog should generate warnings if called with no open log
  if(exists(".loginfo", envir = .GlobalEnv)) rm(".loginfo", envir = .GlobalEnv)
  expect_warning(printlog("hi"))
  expect_warning(closelog())
})

test_that("openlog handles special cases", {
  LOGFILE <- openlog("test", sink = FALSE)
  closelog()

  # Appending
  oldsize <- file.size(LOGFILE)
  openlog("test", sink = FALSE, append = TRUE)
  closelog()
  expect_more_than(file.size(LOGFILE), oldsize)
  file.remove(LOGFILE)

  # Connections - not already open
  test <- gzfile("test.txt.gz")
  LOGFILE <- openlog(test, sink = FALSE)
  expect_true(isOpen(test))
  closelog()
  expect_error(isOpen(test))  # now closed and unavailable
  expect_true(file.exists(LOGFILE))
  file.remove(LOGFILE)

  # Connections - already open
  test <- gzfile("test.txt.gz")
  open(test, "w")
  LOGFILE <- openlog(test, sink = FALSE)
  expect_true(isOpen(test))
  close(test)  # this will screw up closelog()!
  expect_true(file.exists(LOGFILE))
  file.remove(LOGFILE)

  # Detect sink mismatch correctly
  #  capture.output({
  #     LOGFILE <- openlog("test", sink = TRUE)
  #     sink()
  #     printlog("hi")
  #     expect_warning(closelog())
  #  })
  #  file.remove(LOGFILE)
})

test_that("Basic logging works correctly", {
  # opens correctly?
  LOGFILE <- openlog("test", sink = FALSE)
  expect_is(LOGFILE, "character")
  expect_true(file.exists(LOGFILE))
  expect_equal(length(readLines(LOGFILE)), 1)

  # Blank log messages
  expect_true(printlog())
  expect_true(printlog(ts = FALSE, cr = FALSE))

  # log messages written?
  oldsize <- file.size(LOGFILE)
  oldlines <- length(readLines(LOGFILE))
  expect_true(printlog("Line 1"))
  expect_true(printlog(1, "plus", 1, "equals", 1 + 1))
  newsize <- file.size(LOGFILE)
  expect_more_than(newsize, oldsize)
  expect_equal(length(readLines(LOGFILE)) - oldlines, 2)

  # non-simple (character, numeric) objects written?
  expect_true(printlog(cars))
  expect_more_than(file.size(LOGFILE), newsize)

  closelog()
  file.remove(LOGFILE)
})

test_that("logging sinks correctly", {
  utils::capture.output({
    sn <- sink.number()
    LOGFILE <- openlog("test", sink = TRUE)
    expect_equal(sink.number(), sn + 1)
    print("line 2")
    closelog(sessionInfo = FALSE)
    expect_equal(sink.number(), sn)
    expect_equal(length(readLines(LOGFILE)), 3)
  })

  file.remove(LOGFILE)
})

test_that("closelog works correctly", {
  # sessionInfo added?
  LOGFILE <- openlog("test", sink = FALSE)
  oldsize <- file.size(LOGFILE)
  expect_is(closelog(), "numeric")
  expect_more_than(file.size(LOGFILE), oldsize)

  # flag information returned?
  LOGFILE <- openlog("test", sink = FALSE)
  expect_equal(closelog(), 0)
  LOGFILE <- openlog("test", sink = FALSE)
  printlog(flag = TRUE)
  expect_equal(closelog(), 1)

  # suppressing sessionInfo data
  LOGFILE <- openlog("test", sink = FALSE)
  closelog()
  oldsize <- file.size(LOGFILE)
  LOGFILE <- openlog("test", sink = FALSE)
  closelog(sessionInfo = FALSE)
  expect_less_than(file.size(LOGFILE), oldsize)

  file.remove(LOGFILE)
})

test_that("Priority levels work correctly", {

  # file.size() is unreliable on Windows until file is closed
  # (or maybe it's file.info(); see top of file)
  skip_on_os("windows")

  LOGFILE <- openlog("test", loglevel = 0, sink = FALSE)

  size0 <- file.size(LOGFILE)
  printlog("Line 1")
  size1 <- file.size(LOGFILE)
  expect_more_than(size1, size0)
  printlog("Line 2", level = -1)
  size2 <- file.size(LOGFILE)
  expect_equal(size2, size1)
  printlog("Line 1", level = 1)
  expect_more_than(file.size(LOGFILE), size2)

  closelog()
  file.remove(LOGFILE)
})

test_that("Directory change OK", {

  # file.size() is unreliable on Windows until file is closed
  # (or maybe it's file.info(); see top of file)
  skip_on_os("windows")

  LOGFILE <- openlog("test", sink = FALSE)

  printlog("Line 1")
  size0 <- file.size(LOGFILE)
  olddir <- getwd()
  setwd(tempdir())
  printlog("Line 2")
  size1 <- file.size(LOGFILE)
  expect_more_than(size1, size0)
  setwd(olddir)
  printlog("Line 3")
  expect_more_than(file.size(LOGFILE), size1)

  closelog()
  file.remove(LOGFILE)
})

test_that("all log files closed on error", {
  skip("don't know how to test this, because testthat traps errors!")

  detach("package:luzlogr", unload=TRUE)
  options(luzlogr.close_on_error = TRUE)
  library(luzlogr)

  nl <- nlogs()
  openlog("test", sink = FALSE)
  expect_error(stop())
  expect_equal(nlogs(), nl)
})
