context("read nexml")

f <- system.file("examples", "trees.xml", package = "RNeXML")
url <- "https://raw.githubusercontent.com/ropensci/RNeXML/master/inst/examples/trees.xml"

test_that("we can read nexml from a file path", {
  nex <- nexml_read(f)
  
  expect_is(nex, "nexml")
  expect_equal(nex@trees@names, "trees")
})

test_that("we can read nexml from a url", {
  nex <- read.nexml(url)
  
  expect_is(nex, "nexml")
  expect_equal(nex@trees@names, "trees")
})

test_that("we can read nexml from a character string of xml", {
  str <- paste0(readLines(f), collapse = "")
  nex <- nexml_read(str)
  
  expect_is(nex, "nexml")
  expect_equal(nex@trees@names, "trees")
})

test_that("we can read nexml from a XMLInternalDocument object", {
  library("httr")
  library("XML")
  x <- xmlParse(content(GET(url)))
  nex <- nexml_read(x)
  
  expect_is(nex, "nexml")
  expect_equal(nex@trees@names, "trees")
})

test_that("we can read nexml from a XMLInternalNode object", {
  library("httr")
  library("XML")
  x <- xmlParse(content(GET(url)))
  nex <- nexml_read(xmlRoot(x))
  
  expect_is(nex, "nexml")
  expect_equal(nex@trees@names, "trees")
})

test_that("alias for nexml_read works", {
  nex <- read.nexml(f)
  
  expect_is(nex, "nexml")
  expect_equal(nex@trees@names, "trees")
})
