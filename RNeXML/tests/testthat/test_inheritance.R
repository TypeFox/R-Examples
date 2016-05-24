context("inheritance")

## FIXME
## Should include expect_that tests, rather than just running without errors.  
## ADD test to show that toggling xml->s4->xml returns IDENTICAL objects, 
## Add tests to check values on some nodes/attributes...

test_that("we can perform simple conversions between NeXML XML and S4", {
  # basic example
  node <- newXMLNode("meta", 
                     attrs = c('xsi:type'="nex:LiteralMeta",
                               id="dict1",
                               property="cdao:has_tag",
                               datatype="xsd:boolean",
                               content="true"), 
                     suppressNamespaceWarning=TRUE)
  n2 <- newXMLNode("node", 
                   attrs = c(about="#n4",
                             label="n4", 
                             id = "n4"), 
                   .children = node)
  
  # check conversions to/from NeXML
  s4 <- as(n2, "node")
  xmlfroms4 <- as(s4, "XMLInternalNode")
  
##  expect_identical(n2, xmlfroms4) #cannot compare two external pointers
  expect_identical(saveXML(n2), saveXML(xmlfroms4))

})

# test_that("We can parse a complete NeXML file and toggle back and forth between XML and S4", {
test_that("Parse a complete NeXML file to a single otu", {
  doc <- xmlParse(system.file("examples", "trees.xml", package="RNeXML"))
  root <- xmlRoot(doc)

  otu <- as(root[["otus"]][[1]], "otu")
  expect_that(otu, is_a("otu"))
  as(otu, "XMLInternalNode")
})


doc <- xmlParse(system.file("examples", "trees.xml", package="RNeXML"))
root <- xmlRoot(doc)



test_that("Parse a complete NeXML file to trees", {
  trees <- as(root[["trees"]], "trees")
  expect_that(trees, is_a("trees"))
  as(trees, "XMLInternalNode")
})

test_that("Parse a complete NeXML file to many otus", {
  otus <- as(root[["otus"]], "otus")
  expect_that(otus, is_a("otus"))
  tt <- as(otus, "XMLInternalNode")
  expect_that(tt, is_a("XMLInternalNode"))
})

test_that("Parse a complete NeXML file to xmlinternalnode", {
  parsed <- as(root, "nexml")
  expect_that(parsed, is_a("nexml"))
  serialized <- as(parsed, "XMLInternalNode")
  expect_that(serialized, is_a("XMLInternalNode"))
})

test_that("Check that values are correct in the otu class element", {
  otu <- as(root[["otus"]][[1]], "otu")
  expect_that(otu@id[[1]], equals("t1"))
  expect_that(otu@label[[1]], equals("species 1"))
  expect_that(otu@meta, is_a("list"))
  expect_that(otu@about, is_identical_to(character(0)))
})

test_that("Check that values are correct in the trees class element", {
  trees <- as(root[["trees"]], "trees")
  expect_that(trees@tree, is_a("ListOftree"))
  expect_that(trees@otus[[1]], equals("tax1"))
  expect_that(trees@id[[1]], equals("Trees"))
  expect_that(trees@label[[1]], equals("TreesBlockFromXML"))
  expect_that(trees@meta, is_a("list"))
  expect_that(trees@about, is_identical_to(character(0)))
})

test_that("Check that values are correct in the otus class element", {
  otus <- as(root[["otus"]], "otus")
  expect_that(otus@otu, is_a("ListOfotu"))
  expect_that(otus@id[[1]], equals("tax1"))
  expect_that(otus@label[[1]], equals("RootTaxaBlock"))
  expect_that(otus@meta, is_a(class=c("list","ListOfmeta")))
  expect_that(otus@about, is_identical_to(character(0)))
})
