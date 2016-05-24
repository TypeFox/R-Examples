context("simmap")


## Make a simmap tree 



test_that("we can coerce an ape::phylo tree with a 
          phytools:simmap extension into nexml", {
  skip_if_not_installed("phytools")
  library("phytools")
  set.seed(10) 
  tree <- rbdtree(b = log(50), d = 0, Tmax = .5)
  Q <- matrix(c(-2, 1, 1, 1, -2 ,1 ,1, 1, -2), 3, 3) 
  rownames(Q) <- colnames(Q) <- c("A", "B", "C") 
  ## Note that state symbols must be integers! factors will be converted
  mtree <- sim.history(tree, Q) 

  cols <- c("red", "blue", "green")
  names(cols) <- rownames(Q)
  nex <- simmap_to_nexml(mtree) 
  expect_is(nex, "nexml")

  phy <- nexml_to_simmap(nex) 

  orig <- plotSimmap(mtree,cols,ftype="off")
  roundtrip <- plotSimmap(phy,cols,ftype="off")

  # checks that the edge mappings are correct 
  expect_equal(mtree$maps, phy$maps)

  
  orig <- as.integer(as.factor(mtree$states[sort(names(mtree$states))]))
  converted <- as.integer(phy$states[sort(names(phy$states))])
  # checks that we got the states slot correct 
  expect_equal(converted, orig)
})



