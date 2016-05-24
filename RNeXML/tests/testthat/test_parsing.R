context("parsing")

# More lower-level parsing tests in inheritance

test_that("We can parse a NeXML file to an S4 RNeXML::tree object", {
  f <- system.file("examples", "trees.xml", package="RNeXML")
  doc <- xmlParse(f)
  root <- xmlRoot(doc)
  nexml <- as(root, "nexml")  ## parse the XML into S4
  expect_is(nexml,"nexml")
})



test_that("We preserve existing namespace", {

  f <- system.file("examples/biophylo.xml", package="RNeXML")
  nex <- nexml_read(f)
  g <- tempfile()
  nexml_write(nex, g)

  expect_true_or_null(nexml_validate(g))

  nex2 <- nexml_read(g)

  ## check the namespaces are added 
  expect_more_than(length(get_namespaces(nex2)), length(get_metadata(nex)))

  ## Check that the new abbreviations are added 

})
