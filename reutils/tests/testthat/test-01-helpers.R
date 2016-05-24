context("Test helper functions")

xml <- XML::xmlParse("data/test.xml")

test_that("'xvalue' and 'xname' work as advertised", {
  expect_equal(xvalue(xml, "//Id"), c("23927984", "23903989"))
  expect_equal(xname(xml, "//IdList/*"), c("Id", "Id"))
  expect_equal(typeof(xvalue(xml, "//Id", as="numeric")), "double")
  expect_equal(typeof(xvalue(xml, "//Id", as="integer")), "integer")
})

test_that("Setting defaults with 'value' work as advertised", {
  expect_equal(xvalue(xml, "//bla", default=NA), NA_character_)
  expect_equal(xvalue(xml, "//bla", default=""), "")
  expect_equal(xvalue(xml, "//bla", default=NULL), NULL)
  expect_equal(xvalue(xml, "//bla", default=0), "0")
  expect_equal(xvalue(xml, "//bla", as="numeric", default=0), 0)
})

test_that("ncbi_retrieval_type works with 'bioproject'", {
  expect_equal(ncbi_retrieval_type('bioproject', rettype = 'docsum'), list(rettype = "docsum", retmode = "xml"))
})
