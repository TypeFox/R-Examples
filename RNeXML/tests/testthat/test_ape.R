context("ape")


test_that("From ape::phylo to RNeXML::nexml object", {
   data(bird.orders)
   expect_is(as(bird.orders, "nexml"), class="nexml") 
})

test_that("We can go from various orderings of ape::phylo to RNeXML::nexml", {
  data(bird.orders)
  nexml <- as(bird.orders, "nexml")
  phy <- as(nexml, "phylo")

  ## Demonstrate that we now have a phylo object
  p <- plot(phy)
  expect_that(plot(phy), is_a("list"))
  expect_that(phy, is_a("phylo"))
})

test_that("From nexml to multiPhylo", {


  # part of base testing, could be replaced with higher level, but why 
  f <- system.file("examples", "trees.xml", package="RNeXML")
  doc <- xmlParse(f)
  root <- xmlRoot(doc)
  nexml <- as(root, "nexml")  ## parse the XML into S4

  ## APE TEST:  Coerce the S4 into phylo S3 object
  expect_warning(phy <- as(nexml, "phylo"), "Multiple trees found, Returning multiPhylo object")
  
  expect_is(phy, "multiPhylo")

})

## This unit test is really not testing ape functions but just the higher-level nexml_write function...
test_that("We can serialize the various versions of the ape format", {
  data(bird.orders)
  nexml <- as(bird.orders, "nexml")
  nexml_write(nexml, file = "test.xml")
  unlink("test.xml")
})





test_that("We can read and write NeXML to phylo and back without edge.lengths", {
           s <- "owls(((Strix_aluco,Asio_otus),Athene_noctua),Tyto_alba);"
           cat(s, file = "ex.tre", sep = "\n")
           owls <- read.tree("ex.tre")
           nexml_write(owls, file = "ex.xml")
           owls2 <- as(nexml_read("ex.xml"), "phylo")
           expect_equal(owls, owls2)
## FIXME what? 
           unlink("ex.tre")
           unlink("ex.xml")
})



test_that("Rooted trees remain rooted on conversions", {
          expect_true(is.rooted(bird.orders))
          expect_true(is.rooted(as(as(bird.orders, "nexml"), "phylo")))
          write.nexml(bird.orders, file = "tmp.xml")
          expect_true(is.rooted(as(read.nexml("tmp.xml"), "phylo")))
          unlink("tmp.xml")
})

phy <- unroot(bird.orders)
test_that("Unrooted trees remain unrooted on conversions", {
  expect_false(is.rooted(phy))
  expect_false(is.rooted(as(as(phy, "nexml"), "phylo")))
  write.nexml(phy, file = "tmp.xml")
  expect_false(is.rooted(as(read.nexml("tmp.xml"), "phylo")))
  unlink("tmp.xml")
})

test_that("We can convert trees with only some edge lengths into ape::phylo", {
          f <- system.file("examples", "some_missing_branchlengths.xml", package="RNeXML")
          expect_warning(a <- as(read.nexml(f), "phylo"), "Multiple trees found, Returning multiPhylo object")
          # We can parse it, goodness knows what anyone will do with it.  Better to hack off the branch lengths or convert to 0, but that's for the user.   
})

