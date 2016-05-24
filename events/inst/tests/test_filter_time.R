context('Check filter_time')

load('./tabari-demo.RData')

test_that("test end", {
  dd <- filter_time(tabari.demo, end="1995-01-06")
  expect_equal(nrow(dd), 43)
  expect_equal(dd$date[1], as.Date("1995-01-01"))
  expect_equal(dd$date[43], as.Date("1995-01-06"))
})

test_that("test start", {
  dd <- filter_time(tabari.demo, start="1995-01-04")
  expect_equal(nrow(dd), 40)
  expect_equal(dd$date[1], as.Date("1995-01-04"))
  expect_equal(dd$date[40], as.Date("1995-01-19"))
})

test_that("test start and end", {
  dd <- filter_time(tabari.demo, end="1995-01-06", start="1995-01-04")
  expect_equal(nrow(dd), 17)
  expect_equal(dd$date[1], as.Date("1995-01-04"))
  expect_equal(dd$date[17], as.Date("1995-01-06"))
})

