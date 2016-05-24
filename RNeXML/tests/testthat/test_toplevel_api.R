context("top level API")


test_that("read.nexml works", {
  ## The short version using an RNeXML API

  f <- system.file("examples", "trees.xml", package="RNeXML")
  nex <- read.nexml(f) # check alias
  expect_is(nex, "nexml")
})


test_that("write.nexml works (from ape::phylo)", {
  ## The short version using an RNeXML API

  data(bird.orders)
  nexml_write(bird.orders, file="example.xml")
  write.nexml(bird.orders, file="example.xml") # check alias too

## Check that that example is valid NeXML
  expect_true_or_null(nexml_validate("example.xml"))
  expect_is(nexml_read("example.xml", "nexml"), "nexml")

  unlink("example.xml") # cleanup

})


test_that("write.nexml can write multiple trees at once ", {
  f <- system.file("examples", "trees.xml", package="RNeXML")
  nex <- nexml_read(f)
  trees <- get_trees(nex)

  ##  We can write a listOfmultiPhylo if the argument is named
  nexml_write(trees = trees, file="example.xml")
  expect_true_or_null(nexml_validate("example.xml"))

  # we can write a multiPhylo (or phylo) by attempting coercion on the first argument instead:  
  nexml_write(trees[[1]], file="example.xml")
  expect_true_or_null(nexml_validate("example.xml"))


  unlink("example.xml") # cleanup

})



test_that("We can get the right level of lists of trees ", {

  f <- system.file("examples", "trees.xml", package="RNeXML")
  nex <- nexml_read(f)

  ## identical methods, Collapses length-1 lists 
  phy <- as(nex, "phylo")  ##
  phy2 <- get_trees(nex)
  phy3 <- nexml_get(nex, "trees")
  expect_identical(phy, phy2)
  expect_identical(phy3, phy2)

  ## Doesn't collapse the length-1 lists, returns list of multiPhylo always: 
  phy <- as(nex, "multiPhyloList")  ##
  phy2 <- get_trees_list(nex)
  phy3 <- nexml_get(nex, "trees_list")
  expect_identical(phy, phy2)
  expect_identical(phy3, phy2)

  ## Collapse to multiPhylo   
  phy <- as(nex, "multiPhylo")  ##
  phy2 <- get_trees(nex) # same because there are two trees in the same `trees` node.  
  expect_identical(phy, phy2)
  phy3 <- nexml_get(nex, "flat_trees") ## FIXME SOMETHING WRONG!
  expect_identical(phy3, phy2)
})



