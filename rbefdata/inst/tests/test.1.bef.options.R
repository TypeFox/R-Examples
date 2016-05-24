context("Check set and get bef.options()")

test_that("the right default options are set", {
  expect_that(bef.options("url"), matches("http://china.befdata.biow.uni-leipzig.de"))
  expect_that(bef.options("tematres_url"), matches("http://tematres.befdata.biow.uni-leipzig.de"))
  expect_that(bef.options("tematres_service_url"), matches("http://tematres.befdata.biow.uni-leipzig.de/vocab/services.php"))
  expect_that(bef.options("download_dir"), matches("downloads"))
})

test_that("the options are settable", {
  bef.options(url = "www.test.a.com")
  bef.options(tematres_url = "http://befdatatesting.biow.uni-leipzig.de:7070/vocab/index.php")
  bef.options(tematres_service_url = "http://befdatatesting.biow.uni-leipzig.de:7070/vocab/services.php")
  bef.options(download_dir = "testdir")

  expect_that(bef.options("url"), matches("www.test.a.com"))
  expect_that(bef.options("tematres_url"), matches("http://befdatatesting.biow.uni-leipzig.de:7070/vocab/index.php"))
  expect_that(bef.options("tematres_service_url"), matches("http://befdatatesting.biow.uni-leipzig.de:7070/vocab/services.php"))
  expect_that(bef.options("download_dir"), matches("testdir"))
})
