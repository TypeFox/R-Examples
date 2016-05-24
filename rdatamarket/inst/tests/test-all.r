context("DataMarket dataset info")

test_that("Dataset 17tm looks as expected", {
  info <- dminfo('17tm')
  expect_equal(info$id, '17tm')
  expect_equal(info$ds, '17tm')
  expect_equal(info$title, 'Oil: Production tonnes')
  expect_is(info$dimensions, 'list')
  dim1 <- info$dimensions[[1]]
  expect_equal(dim1$id, 'kqc')
  expect_equal(dim1$title, 'Country')
  expect_equal(dim1$type, 'simple')
  expect_is(dim1$values, 'list')
  expect_equal(dim1$values[[1]][["id"]], 'a')
  expect_equal(dim1$values[[1]][["title"]], 'Algeria')
  expect_equal(dim1$values[[2]][["id"]], '17')
  expect_equal(dim1$values[[2]][["title"]], 'Angola')
  expect_equal(dim1$values[[3]][["id"]], 'd')
  expect_equal(dim1$values[[3]][["title"]], 'Argentina')
})

test_that("Info is kept around if passed to interpret_ds", {
  mock <- list(id='17tm', ds='17tm', meta='dummy', dimensions='dummy',
    status='dummy', title='dummy')
  expect_identical(interpret_ds(mock), list(
    base='http://datamarket.com',
    qs=list(ds='17tm'),
    infos=list(`17tm`=mock)
    ))

  mock2 <- list(id='foo', ds='foo', meta='dummy', dimensions='dummy',
    status='dummy', title='dummy')
  expect_identical(interpret_ds(list(`17tm`=mock, foo=mock2)), list(
    base='http://datamarket.com',
    qs=list(ds='17tm/foo'),
    infos=list(`17tm`=mock, foo=mock2)
    ))

  # and also if unnamed
  expect_identical(interpret_ds(list(mock, mock2)), list(
    base='http://datamarket.com',
    qs=list(ds='17tm/foo'),
    infos=list(`17tm`=mock, foo=mock2)
    ))
})

test_that("Info is reused if given to dminfo", {
  mock <- list(id='17tm', ds='17tm', meta='dummy', dimensions='dummy',
    status='dummy', title='dummy')
  expect_identical(dminfo(mock), list(`17tm`=mock))
  mock <- list(`17tm`=mock)
  expect_identical(dminfo(mock), mock)
})

context("DataMarket timeseries")

test_that("Timeseries from dataset 17tm works", {
  series <- dmseries('17tm!kqc=a.17.d')
  expect_is(series, 'zoo')
  expect_identical(names(series), c('Algeria', 'Angola', 'Argentina'))
  expect_equal(as.numeric(series[1]), c(26.481, 0.655, 13.7647586207))
  expect_equal(as.numeric(series[2]), c(33.872, 0.631, 14.6439655172))
  times <- index(series)
  expect_identical(times[1], 1965L)
  expect_identical(times[2], 1966L)
})

test_that("Timeseries from dataset 17tm with old DS format works", {
  series <- dmseries('17tm|kqc=a.17.d')
  expect_is(series, 'zoo')
  expect_identical(names(series), c('Algeria', 'Angola', 'Argentina'))
  expect_equal(as.numeric(series[1]), c(26.481, 0.655, 13.7647586207))
  expect_equal(as.numeric(series[2]), c(33.872, 0.631, 14.6439655172))
  times <- index(series)
  expect_identical(times[1], 1965L)
  expect_identical(times[2], 1966L)
})

test_that("Timeseries with parameter dimension filtering works", {
  series <- dmseries('17tm', Country='Algeria')
  expect_is(series, 'zoo')
  expect_equal(as.numeric(series[1]), c(26.481))
  expect_equal(as.numeric(series[2]), c(33.872))
  times <- index(series)
  expect_identical(times[1], 1965L)
  expect_identical(times[2], 1966L)

  # or, in short:
  series_from_param <- dmseries('17tm', Country='Algeria')
  series_from_dsstring <- dmseries('17tm!kqc=a')
  expect_identical(series_from_param, series_from_dsstring)

})

test_that("Timeseries with multi-valued parameter dimension filtering works", {
  series_from_param <- dmseries('17tm', Country=c('Algeria', 'Angola',
    'Argentina'))
  series_from_dsstring <- dmseries('17tm!kqc=a.17.d')
  expect_identical(series_from_param, series_from_dsstring)
})

test_that("Timeseries with quarter granularity works", {
  series <- dmseries("1k8d")
  expect_is(series, "zoo")
  expect_identical(names(series), c("Income.Payments..Quarterly.Data."))
  expect_identical(as.numeric(series[1]), c(-0.331))
  expect_identical(as.numeric(series[2]), c(-0.314))
  times <- index(series)
  expect_identical(class(times), "yearqtr")
  expect_identical(times[1], as.yearqtr("1960 Q1"))
  expect_identical(times[2], as.yearqtr("1960 Q2"))
})

test_that("Timeseries with month granularity works", {
  series <- dmseries("1k62")
  expect_is(series, "zoo")
  expect_identical(names(series), "MZM.Own.Rate")
  expect_identical(as.numeric(series[1]), 2.711)
  expect_identical(as.numeric(series[2]), 2.705)
  times <- index(series)
  expect_identical(class(times), "yearmon")
  expect_identical(times[1], as.yearmon("1974-01"))
  expect_identical(times[2], as.yearmon("1974-02"))
})

test_that("Timeseries with week granularity works", {
  series <- dmseries("1lpt")
  expect_is(series, "zoo")
  expect_identical(names(series), "U.S..Economic.Statistics..M2.Money.Stock")
  expect_identical(as.numeric(series[1]), 1601.8)
  expect_identical(as.numeric(series[2]), 1595.2)
  times <- index(series)
  expect_identical(class(times), "factor")
  expect_equal(as.character(times[1]), "1980-W01")
  expect_equal(as.character(times[2]), "1980-W45")
})

test_that("Timeseries with date granularity works", {
  series <- dmseries("1k2e")
  expect_is(series, "zoo")
  expect_identical(names(series), "U.S....Euro.Foreign.Exchange.Rate")
  expect_identical(as.numeric(series[1]), 1.1812)
  expect_identical(as.numeric(series[2]), 1.1760)
  times <- index(series)
  expect_identical(class(times), "Date")
  expect_identical(times[1], as.Date("1999-01-04"))
  expect_identical(times[2], as.Date("1999-01-05"))
})

context("DataMarket long-form data ('list')")

test_that("Long-form data from dataset 17tm works", {
  lis <- dmlist('ds=17tm!kqc=a.17.d&maxdate=2012')
  lis <- lis[order(lis$Country, lis$Year),]
  expect_is(lis, 'data.frame')
  expect_identical(names(lis), c('Country', 'Year', 'Value'))
  expect_identical(as.character(lis$Country), c(
    replicate(48, 'Algeria'),
    replicate(48, 'Angola'),
    replicate(48, 'Argentina')
  ))
  expect_identical(lis$Year, c(replicate(3, 1965:2012)))
  expect_equal(lis$Value[1:4], c(26.481, 33.872, 39.076, 42.904))
})

test_that("Long-form data from dataset 17tm with old DS format works", {
  lis <- dmlist('ds=17tm|kqc=a.17.d&maxdate=2012')
  lis <- lis[order(lis$Country, lis$Year),]
  expect_is(lis, 'data.frame')
  expect_identical(names(lis), c('Country', 'Year', 'Value'))
  expect_identical(as.character(lis$Country), c(
    replicate(48, 'Algeria'),
    replicate(48, 'Angola'),
    replicate(48, 'Argentina')
  ))
  expect_identical(lis$Year, c(replicate(3, 1965:2012)))
  expect_equal(lis$Value[1:4], c(26.481, 33.872, 39.076, 42.904))
})

test_that("Long-form data with parameter dimension filtering works", {
  list_from_param <- dmlist('17tm', Country='Algeria')
  list_from_dsstring <- dmlist('17tm!kqc=a')
  expect_identical(list_from_param, list_from_dsstring)
})

test_that("Long-form with multi-valued parameter dimension filtering works", {
  list_from_param <- dmlist('17tm', Country=c('Algeria', 'Angola', 'Argentina'))
  list_from_dsstring <- dmlist('17tm!kqc=a.17.d')
  expect_identical(list_from_param, list_from_dsstring)
})

context("Formatting functions")

test_that("Formatting of dmdataset object works", {
  ds <- dminfo("17tm")
  expect_equal(format(ds), 'Title: "Oil: Production tonnes"
Provider: "BP"
Dimensions:
  "Country" (60 values):
    "Algeria"
    "Angola"
    "Argentina"
    "Australia"
    "Azerbaijan"
    [...]')
})

test_that("Formatting of dmdimension object works", {
  ds <- dminfo("17tm")
  expect_equal(format(ds$dimensions[[1]]), '"Country" (60 values):
    "Algeria"
    "Angola"
    "Argentina"
    "Australia"
    "Azerbaijan"
    [...]')
})

test_that("Formatting of hierarchical dmdimension object works", {
  ds <- dminfo("1hy5")
  expect_equal(format(ds$dimensions[[1]]), '"Country" (5 values):
    "Total World"
    -> "OECD"
    -> "NonOECD"
    -> "European Union"
    -> "Former Soviet Union"')
})

test_that("Formatting of dmdimvalues works", {
  ds <- dminfo("17tm")
  expect_equal(format(ds$dimensions[[1]]$values[[1]]), '"Algeria"')
  expect_equal(format(ds$dimensions[[1]]$values[2:4]), ' 17  "Angola"
  d  "Argentina"
  z  "Australia"')
})

test_that("Formatting of dmhierarchicaldimvalue object works", {
  ds <- dminfo("1hy5")
  expect_equal(format(ds$dimensions[[1]]$values[[1]]), '"Total World"')
  expect_equal(format(ds$dimensions[[1]]$values[[2]]), '-> "OECD"')
  expect_equal(format(ds$dimensions[[1]]$values[[3]]), '-> "NonOECD"')
})

context("Util functions")

test_that("dimfilter forms DS strings from named args correctly", {
  mockinfos <- list(`17tm`=list(dimensions=list(kqc=list(
    id='kqc', title='Country', values=list(
      a=list(id='a', title='Algeria'),
      `17`=list(id='17', title='Angola'),
      d=list(id='d', title='Argentina')
    )
  ))))
  result <- dimfilter('17tm', mockinfos, Country='Algeria')
  expect_equal(result, '17tm!kqc=a')
  result <- dimfilter('17tm', mockinfos, Country=c('Angola', 'Argentina'))
  expect_equal(result, '17tm!kqc=17.d')
})

