# Testing code for internal functions

context("internals")

test_that("functions handle bad input", {
  expect_error(newlog(logfile = 1))
  expect_error(newlog(loglevel = '1'))
  expect_error(newlog(sink = 1))
  expect_error(newlog(description = 1))
  expect_error(newlog(closeit = 1))

  expect_error(setlogdata(datum = 1))
})


test_that("newlog and nlogs work correctly", {
  nl <- nlogs()
  expect_equal(nl, 0)

  logfile <- file("logfile")
  loglevel <- 0
  sink <- TRUE
  description <- "description"
  closeit <- TRUE

  # Make sure new entry put on stack
  newlog(logfile = logfile, loglevel = loglevel, sink = sink,
         description = description, closeit = closeit)
  expect_equal(nlogs(), nl + 1)

  # Make sure information stored correctly
  loginfo <- getloginfo()
  expect_is(loginfo, "list")

  expect_equal(logfile, loginfo$logfile)
  expect_equal(loglevel, loginfo$loglevel)
  expect_equal(sink, loginfo$sink)
  expect_equal(description, loginfo$description)
  expect_equal(closeit, loginfo$closeit)

  removelog()
  expect_warning(getloginfo())
  close(logfile)
})

test_that("removelog works correctly", {
  nl <- nlogs()
  expect_equal(nl, 0)

  expect_warning(removelog())

  logfile <- file("logfile")
  loglevel <- 0
  sink <- TRUE
  description <- "description"
  closeit <- TRUE

  newlog(logfile = logfile, loglevel = loglevel, sink = sink,
         description = description, closeit = closeit)

  # Make sure information removed correctly
  loginfo <- removelog()
  expect_is(loginfo, "list")

  expect_equal(logfile, loginfo$logfile)
  expect_equal(loglevel, loginfo$loglevel)
  expect_equal(sink, loginfo$sink)
  expect_equal(description, loginfo$description)
  expect_equal(closeit, loginfo$closeit)

  expect_equal(nlogs(), 0)
  close(logfile)
})

test_that("setlogdata works correctly", {
  logfile <- file("logfile")
  loglevel <- 0
  sink <- TRUE
  description <- "description"
  closeit <- TRUE

  newlog(logfile = logfile, loglevel = loglevel, sink = sink,
         description = description, closeit = closeit)

  # Only allowed to modify one aspect of log info
  expect_error(setlogdata("logfile", "newlogfile"))
  expect_error(setlogdata("loglevel", 1))
  expect_error(setlogdata("sink", FALSE))
  expect_error(setlogdata("description", "newdescription"))
  expect_null(setlogdata("flags", 1))

  # Make sure information changed correctly
  loginfo <-   removelog()
  expect_equal(loginfo$flags, 1)
  close(logfile)
})
