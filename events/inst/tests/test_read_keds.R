context('Check read_keds')

acts <- c("ARN","BIL","BRE","CAL","DAG","ENT","ERI","FOR","FRO","GON",
            "ITH","MOR","NGO","ORC","OSG","ROH", "SAM","UNO")
targs <- c("BIL","BRE","CAL","DAG","ENT","ERI","FOR","FRO","GON","ITH","ORC",
             "OSG","ROH","SAM","UNO")
srcs <- c("ARN","BIL","BRE","CAL","DAG","ERI","FOR","FRO","GON","ITH", 
            "MOR","NGO","ORC","OSG","SAM")
cds <- c("023","024","031","032","033","041","042","044","045","050",
           "060","064","066","071","081","085","095","102","112","113",
           "122","141","173","191","222","223","224")

test_that("read_keds works without one a day", {
  dd <- read_keds('./tabari-demo.keds', one.a.day=FALSE)
  expect_equal(nrow(dd), 77)
  expect_that(dd, is_a('eventdata'))
  expect_that(dd$date, is_a('Date'))
  expect_that(dd$source, is_a('factor'))
  expect_that(dd$target, is_a('factor'))
  expect_that(dd$code, is_a('factor'))

  expect_identical(actors(dd), acts)
  expect_identical(sources(dd), srcs)
  expect_identical(targets(dd), targs)
  expect_identical(codes(dd), cds)
  expect_identical(dd$date[1], as.Date("1995-01-01"))
  expect_identical(dd$date[length(dd$date)], as.Date("1995-01-19"))
})

test_that("read_keds works with with one a day", {
  dd <- read_keds('./tabari-demo.keds')
  expect_equal(nrow(dd), 66) ## duplicates removed
  expect_that(dd, is_a('eventdata'))
  expect_that(dd$date, is_a('Date'))
  expect_that(dd$source, is_a('factor'))
  expect_that(dd$target, is_a('factor'))
  expect_that(dd$code, is_a('factor'))

  expect_identical(actors(dd), acts)
  expect_identical(sources(dd), srcs)
  expect_identical(targets(dd), targs)
  expect_identical(codes(dd), cds)
  expect_identical(dd$date[1], as.Date("1995-01-01"))
  expect_identical(dd$date[length(dd$date)], as.Date("1995-01-19"))
})

test_that("testing one_a_day filter", {
  load('./tabari-demo.RData')

  dd <- read_keds('./tabari-demo.keds', one.a.day=FALSE)
  dd2 <- one_a_day(dd)
  expect_identical(dd2, tabari.demo)
})

test_that("read_keds works with multiple files", {
  dd <- read_keds(c('./tabari-demo.keds', './tabari-demo.keds'), one.a.day=TRUE)
  expect_equal(nrow(dd), 66)
  dd <- read_keds(c('./tabari-demo.keds', './tabari-demo.keds'), one.a.day=FALSE)
  expect_equal(nrow(dd), 154)
})

