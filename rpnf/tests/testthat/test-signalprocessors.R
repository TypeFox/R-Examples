context(desc="Test xo.signalprocessor() function")

test_that(desc="Test, if xo.signalprocessor() throws errors/warnings on wrong arguments",
{
  data <- data.frame()
  expect_error(object=xo.signalprocessor(data))
})

test_that(desc="Test, if xo.signalprocessor() detects DOUBLE_TOP correctly", {
  data <- read.csv("runit-testcase-signalprocessor-doubletop.csv",colClasses=c("Date","integer","character","integer","character"))
  expect_equal(object=xo.signalprocessor(data)$signal.bs,expected=data$result.signal.bs)
})
test_that(desc="Test, if xo.signalprocessor() detects BULLISH_SIGNAL correctly", {
  data <- read.csv("runit-testcase-signalprocessor-bullishsignal.csv",colClasses=c("Date","integer","character","integer","character"))
  expect_equal(object=xo.signalprocessor(data)$signal.bs,expected=data$result.signal.bs)
})
test_that(desc="Test, if xo.signalprocessor() detects TRIPLE_BULLISH_SIGNAL correctly", {
  data <- read.csv("runit-testcase-signalprocessor-triplebullishsignal.csv",colClasses=c("Date","integer","character","integer","character"))
  expect_equal(object=xo.signalprocessor(data)$signal.bs,expected=data$result.signal.bs)
})
test_that(desc="Test, if xo.signalprocessor() detects BULLISH_CATAPULT correctly", {
  data <- read.csv("runit-testcase-signalprocessor-bullishcatapult.csv",colClasses=c("Date","integer","character","integer","character"))
  expect_equal(object=xo.signalprocessor(data)$signal.bs,expected=data$result.signal.bs)
})
test_that(desc="Test, if xo.signalprocessor() detects BULLISH_TRIANGLE correctly", {
  data <- read.csv("runit-testcase-signalprocessor-bullishtriangle.csv",colClasses=c("Date","integer","character","integer","character"))
  expect_equal(object=xo.signalprocessor(data)$signal.bs,expected=data$result.signal.bs)
})
test_that(desc="Test, if xo.signalprocessor() detects BEARISH_SIGNAL_REVERSED correctly", {
  data <- read.csv("runit-testcase-signalprocessor-bearishsignalreversed.csv",colClasses=c("Date","integer","character","integer","character"))
  expect_equal(object=xo.signalprocessor(data)$signal.bs,expected=data$result.signal.bs)
})

# TODO low pole is a tricky one!!!
# test.signalprocessor.LOW_POLE <- function() {
#   data <- read.csv("runit-testcase-signalprocessor-lowpole.csv",colClasses=c("Date","integer","character","integer","character"))
#   checkEquals(data$result.signal.bs,xo.signalprocessor(data)$signal.bs)
# }
# test.signalprocessor.BEAR_TRAP <- function() {
#   data <- read.csv("runit-testcase-signalprocessor-beartrap.csv",colClasses=c("Date","integer","character","integer","character"))
#   checkEquals(data$result.signal.bs,xo.signalprocessor(data)$signal.bs)
# }

#
# bearish signals
test_that(desc="Test, if xo.signalprocessor() detects DOUBLE_BOTTOM correctly", {
  data <- read.csv("runit-testcase-signalprocessor-doublebottom.csv",colClasses=c("Date","integer","character","integer","character"))
  expect_equal(object=xo.signalprocessor(data)$signal.bs,expected=data$result.signal.bs)
})
test_that(desc="Test, if xo.signalprocessor() detects BEARISH_SIGNAL correctly", {
  data <- read.csv("runit-testcase-signalprocessor-bearishsignal.csv",colClasses=c("Date","integer","character","integer","character"))
  expect_equal(object=xo.signalprocessor(data)$signal.bs,expected=data$result.signal.bs)
})
test_that(desc="Test, if xo.signalprocessor() detects TRIPLE_BOTTOM correctly", {
  data <- read.csv("runit-testcase-signalprocessor-triplebottom.csv",colClasses=c("Date","integer","character","integer","character"))
  expect_equal(object=xo.signalprocessor(data)$signal.bs,expected=data$result.signal.bs)
})
test_that(desc="Test, if xo.signalprocessor() detects TRIPLE_BEARISH_SIGNAL correctly", {
  data <- read.csv("runit-testcase-signalprocessor-triplebearishsignal.csv",colClasses=c("Date","integer","character","integer","character"))
  expect_equal(object=xo.signalprocessor(data)$signal.bs,expected=data$result.signal.bs)
})
test_that(desc="Test, if xo.signalprocessor() detects BEARISH_CATAPULT correctly", {
  data <- read.csv("runit-testcase-signalprocessor-bearishcatapult.csv",colClasses=c("Date","integer","character","integer","character"))
  expect_equal(object=xo.signalprocessor(data)$signal.bs,expected=data$result.signal.bs)
})
test_that(desc="Test, if xo.signalprocessor() detects BEARISH_TRIANGLE correctly", {
  data <- read.csv("runit-testcase-signalprocessor-bearishtriangle.csv",colClasses=c("Date","integer","character","integer","character"))
  expect_equal(object=xo.signalprocessor(data)$signal.bs,expected=data$result.signal.bs)
})
test_that(desc="Test, if xo.signalprocessor() detects BULLISH_SIGNAL_REVERSED correctly", {
  data <- read.csv("runit-testcase-signalprocessor-bullishsignalreversed.csv",colClasses=c("Date","integer","character","integer","character"))
  expect_equal(object=xo.signalprocessor(data)$signal.bs,expected=data$result.signal.bs)
})

# test.signalprocessor.HIGH_POLE <- function() {
#   data <- read.csv("runit-testcase-signalprocessor-highpole.csv",colClasses=c("Date","integer","character","integer","character"))
#   checkEquals(data$result.signal.bs,xo.signalprocessor(data)$signal.bs)
# }
# test.signalprocessor.BULL_TRAP <- function() {
#   data <- read.csv("runit-testcase-signalprocessor-bulltrap.csv",colClasses=c("Date","integer","character","integer","character"))
#   checkEquals(data$result.signal.bs,xo.signalprocessor(data)$signal.bs)
# }
