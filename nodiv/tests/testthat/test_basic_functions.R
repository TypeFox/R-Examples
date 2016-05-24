context("basic functionality")

data(coquettes)

test_that("Node traversion",{
  expect_equal(basal_node(coquettes), 25)
  expect_equal(Descendants(28, coquettes), c(29,31)) 
  expect_equal(Parent(39, coquettes), 31)
  expect_equal(Sister(27, coquettes), 21)
})

test_that("nodenumbers()", {
  expect_equal(nodenumbers(coquettes), nodenumbers(coquettes$phylo))
  expect_equal(sum(nodenumbers(coquettes)), 828) 
  expect_equal(range(nodenumbers(coquettes)), c(25, 47)) 
})

test_that("nodes()", {
  expect_warning(nodes(coquettes), "no node labels")
  expect_equal(nodes(coquettes, TRUE), nodenumbers(coquettes))
})

test_that("Nspecies() and Nsites()", {
  expect_equal(Nspecies(coquettes), 24)
  expect_equal(Nsites(coquettes), 154)
})

test_that("sites()",{
  sit <- sites(coquettes)
  expect_equal(length(sit), 154)
  expect_equal(head(sit), c("5", "7", "14", "17", "18", "26"))
  expect_equal(tail(sit), c("699", "700", "704", "705", "715", "719"))
})

test_that("species()", {
  spe <- species(coquettes)
  expect_equal(length(spe), 24)
  expect_equal(spe[4], "Chalcostigma_herrani")
  expect_equal(class(spe), "character")
})

test_that("richness()", {
  ric <- richness(coquettes)
  expect_equal(class(ric), "numeric")
  expect_equal(length(ric), 154)
  expect_named(ric, sites(coquettes))
  expect_equivalent(head(ric), c(4,4,1,2,3,3))
  expect_equivalent(tail(ric), c(2,3,3,4,2,3))
  expect_equivalent(richness(coquettes, c("700", "142")), c(3, 3))
  expect_equivalent(richness(coquettes, c(5, 123)), c(3, 2))
})

test_that("occupancy()", {
  occup <- occupancy(coquettes)
  expect_equal(class(occup), "numeric")
  expect_equal(length(occup), 24)
  expect_named(occup, species(coquettes))
  expect_equal(sum(head(occup)), 182)
  expect_equal(sum(tail(occup)), 30)  
  expect_equal(sum(occupancy(coquettes, c("Heliangelus_strophianus", "Discosura_popelairii"))), 7)
  expect_equal(sum(occupancy(coquettes, c(5, 23))), 3)
})

