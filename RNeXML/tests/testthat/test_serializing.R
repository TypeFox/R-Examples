context("serializing")


## More tests at lower-level serializing from S4 to XML in inheritance.R

test_that("We can serialize ape to S4 RNeXML into valid NeXML",{
  data(bird.orders)


  nexml <- as(bird.orders, "nexml") 

  as(nexml, "XMLInternalNode")
  ###  Higher level API tests
  nexml_write(bird.orders, file="test.xml")
  expect_true_or_null(nexml_validate("test.xml"))

 ##  Clean up
  unlink("test.xml")

  })


test_that("We can serialize parsed NeXML to S4 RNeXML into valid NeXML",{
  root <- xmlRoot(xmlParse(system.file("examples", "trees.xml", package="RNeXML")))
  tree <- as(root, "nexml")
  nexml_write(tree, file="test.xml")

  ## validate
  expect_true_or_null(nexml_validate("test.xml"))

  ##  Clean up
  unlink("test.xml")

  })



#root <- xmlRoot(xmlParse(system.file("examples", "trees.xml", package="RNeXML")))
#tree <- as(root, "nexml")
#tree@trees[[1]]@tree[[1]]@node[[4]]@meta
#as(root[["trees"]][["tree"]][[4]][["meta"]], "meta")


