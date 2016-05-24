context("utils")

test_that("find_exclusive_clade", {
  phy <- read.tree("pinales.tre")
  exclude <- c(grep("Pinus_", phy$tip.label, value=TRUE),
               grep("Larix_", phy$tip.label, value=TRUE))
  sp <- find_exclusive_clade("Abies_recurvata", exclude, phy)

  expect_that(intersect(sp$species, exclude), equals(character(0)))
  expect_that(sp$node, equals(188))
  desc <- phy$tip.label[get_descendants(sp$node, phy, tips.only=TRUE)]
  expect_that(sp$species, equals(desc))
})
