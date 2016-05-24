context("concatenate method")


test_that("we can concatenate two files with unique ids", {
 f1 <- system.file("examples", "trees.xml", package="RNeXML")
 f2 <- system.file("examples", "comp_analysis.xml", package="RNeXML")
 nex1 <- read.nexml(f1)
 nex2 <- read.nexml(f2)
 expect_is(c(nex1, nex2), "nexml")

})


test_that("we get an error if the files to be concatenated have non-unique ids", {
 f1 <- system.file("examples", "trees.xml", package="RNeXML")
 nex1 <- read.nexml(f1)
 nex2 <- read.nexml(f1)
 expect_error(c(nex1, nex2),"ids are not unique")
})



test_that("we can concatenate meta elements", {
  out <-  c(meta(content="example", property="dc:title"),
            meta(content="Carl", property="dc:creator"))
  expect_is(out, "ListOfmeta")
  sapply(out, expect_is, "meta")

})


test_that("we can conatenate meta elements with empty ListOfmeta elements", {

  ## Doesn't trigger our method if x is not class `meta`
  out <- c(new("ListOfmeta"), 
           meta(content="example", property="dc:title"),
           meta(content="Carl", property="dc:creator"))

  expect_is(out, "ListOfmeta")
  sapply(out, expect_is, "meta")
  
## in any order 
  out <-  c(meta(content="example", property="dc:title"),
            meta(content="Carl", property="dc:creator"),
            new("ListOfmeta"))

  expect_is(out, "ListOfmeta")
  sapply(out, expect_is, "meta")

})


test_that("we can conatenate meta elements with ListOfmeta elements", {
out <- c(meta(content="example", property="dc:title"),
         meta(content="Carl", property="dc:creator"))

  out <- c(out, meta("skos:note", "an editorial note"))

  expect_is(out, "ListOfmeta")
  sapply(out, expect_is, "meta")

  
  out <- c(meta("skos:note", "another editorial note"), out)

  expect_is(out, "ListOfmeta")
  sapply(out, expect_is, "meta")

})


test_that("we can concatenate two ListOfmeta elements", {

  metalist <- c(meta(content="example", property="dc:title"),
                meta(content="Carl", property="dc:creator"))
  out <- c(metalist, metalist) 
  expect_is(out, "ListOfmeta")
  expect_is(out[[1]], "meta")
  expect_equal(length(out), 4)
})

test_that("we can concatenate a ListOfmeta and a meta", {

  metalist <- c(meta(content="example", property="dc:title"),
                meta(content="Carl", property="dc:creator"))

  out <- c(metalist, meta(content="a", property="b")) 
  expect_is(out, "ListOfmeta")
  expect_is(out[[1]], "meta")
  expect_equal(length(out), 3)

})


test_that("we can read in a file with existing meta and append without overwriting", 
          {
            f <- system.file("examples/biophylo.xml", package="RNeXML")
            nex <- nexml_read(f)
            g <- tempfile()
            nexml_write(nex, g)
            nex2 <- nexml_read(g)
            expect_more_than(length(get_metadata(nex2)), length(get_metadata(nex)))
          })

