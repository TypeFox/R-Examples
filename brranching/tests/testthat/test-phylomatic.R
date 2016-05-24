context("phylomatic")

taxa <- c("Poa annua", "Phlox diffusa", "Helianthus annuus")

test_that("phylomatic - GET (default) method works", {
  skip_on_cran()

  tree <- phylomatic(taxa = taxa)

  expect_is(taxa, "character")
  expect_is(tree, "phylo")
  expect_is(tree, "phylomatic")
})


test_that("phylomatic - POST method works", {
  skip_on_cran()

  tree <- phylomatic(taxa = taxa, get = 'POST')

	expect_is(tree, "phylo")
	expect_is(tree, "phylomatic")
})

test_that("phylomatic - stored tree", {
  skip_on_cran()

  tree <- phylomatic(taxa, storedtree = 'smith2011')

  expect_is(tree, "phylo")
  expect_is(tree, "phylomatic")
})

test_that("phylomatic - nexml output format", {
  skip_on_cran()

  out <- phylomatic(taxa, outformat = "nexml")

  expect_is(out, "phylomatic")
  expect_is(out[1], "character")
  expect_true(grepl("nexml", out))
})

test_that("phylomatic - treeuri param", {
  skip_on_cran()

  spp <- c("Abies_nordmanniana", "Abies_bornmuelleriana", "Abies_cilicica", "Abies_cephalonica",
    "Abies_numidica", "Abies_pinsapo", "Abies_alba")
  url <- "http://datadryad.org/bitstream/handle/10255/dryad.8791/final_tree.tre?sequence=1"
  tree <- phylomatic(taxa = spp, treeuri = url)

  expect_is(tree, "phylo")
  expect_is(tree, "phylomatic")
  expect_equal(tree$tip.label, spp)
})


test_that("phylomatic fails as expected", {
  skip_on_cran()

  # fails when no taxa given
	expect_error(phylomatic(), "argument \"taxa\" is missing")

	# fails when no taxnames FALSE and improper name strings passed
	expect_error(phylomatic(taxa, taxnames = FALSE), "No taxa in common")

	# fails when get isn't in allowed set
	expect_error(phylomatic(taxa, get = "STUFF"), "get must be one of 'POST' or 'GET'")

	# fails when too many taxa passed and get='GET'
	library("taxize")
	spp <- taxize::names_list("species", 200)
	expect_error(phylomatic(spp, get = "GET"), "\\(414\\) Request-URI Too Long")
})
