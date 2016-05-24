#
# --- Test setAs-Methods.R ---
#

### Get all the test files
if (Sys.getenv("RCMDCHECK") == FALSE) {
    pth <- file.path(getwd(), "..", "inst", "nexmlfiles")
} else {
    pth <- system.file(package="phylobase", "nexmlfiles")
}

# create ape::phylo version of a simple tree for testing
nwk <- "((t1:0.1,t2:0.2)n7:0.7,(t3:0.3,(t4:0.4,t5:0.5)n9:0.9)n8:0.8)n6:0.6;"
tr <- ape::read.tree(text=nwk)

# create analogous phylo4 object with a full complement of valid slots
ancestor <- as.integer(c(6,7,7,6,8,0,8,9,9))
descendant <- as.integer(c(7,1,2,8,3,6,9,4,5))
edge <- cbind(ancestor, descendant)
nid.tip <- 1:5
nid.int <- 6:9
lab.tip <- paste("t", nid.tip, sep="")
lab.int <- paste("n", nid.int, sep="")
elen <- descendant/10
elab <- paste("e", ancestor, descendant, sep="-")
phy <- phylo4(x=edge, tip.label=lab.tip, node.label=lab.int,
    edge.length=elen, edge.label=elab)

# create altered version such that each slot is out of order with
# respect to all others; methods should be able to handle this
phy.alt <- phy
phy.alt@label <- rev(phy@label)
phy.alt@edge <- phy@edge[c(6:9, 1:5), ]
phy.alt@edge.length <- phy@edge.length[c(7:9, 1:6)]
phy.alt@edge.label <- phy@edge.label[c(8:9, 1:7)]

## NeXML files
compFile <- file.path(pth, "comp_analysis.xml")
stopifnot(file.exists(compFile))

#-----------------------------------------------------------------------

context("setAs methods")

test_that("phylo to phylo4", {
    # simple case
    as.phy <- as(tr, "phylo4")
    expect_true(class(as.phy)=="phylo4")
    expect_equal(tr$edge, unname(edges(as.phy, drop.root=TRUE)))
    expect_equal(tr$tip.label, unname(tipLabels(as.phy)))
    expect_equal(tr$node.label, unname(nodeLabels(as.phy)))
    # TODO: ape keeps the root edge length in $root.edge
    #expect_equal(tr$edge.length, unname(edgeLength(as.phy)))
    expect_equal("preorder", edgeOrder(as.phy))

    ## test preservation of order attribute
    as.phy <- as(reorder(tr, "cladewise"), "phylo4")
    expect_equal("preorder", edgeOrder(as.phy))
    as.phy <- as(reorder(tr, "pruningwise"), "phylo4")
    expect_equal("postorder", edgeOrder(as.phy))

    ## test phylo import when only 2 tips
    tr2 <- ape::drop.tip(tr, 3:ape::Ntip(tr))
    expect_equal(nTips(as(tr2, "phylo4")), 2)
    expect_equal(nNodes(as(tr2, "phylo4")), 1)

    ## simple roundtrip test
    phy <- as(tr, "phylo4")
    expect_equal(tr, as(phy, "phylo"))
})

# note: this method mostly just wraps phylo->phylo4 coercion (tested
# above) and phylo4d("phylo4") method (tested in runit.class-phylo4d.R)
test_that("phylo to phylo4d", {
    expect_equal(as(tr, "phylo4d"), phylo4d(tr))
    phyd <- as(tr, "phylo4d")
    expect_true(class(phyd)=="phylo4d")
    # simple roundtrip test
    phyd <- as(tr, "phylo4d")
    expect_warning(phyo <- as(phyd, "phylo"))
    expect_equal(tr, phyo)
})

## test.multiPhylo.As.multiPhylo4 <- function() {
## }

## test.multiPhylo4.As.multiPhylo <- function() {
## }

test_that("nexml to phylo4", {
    nxml <- RNeXML::nexml_read(compFile)
    phy4 <- as(nxml, "phylo4")
    expect_true(all(tipLabels(phy4) %in% paste("taxon", 1:10, sep="_")))
    expect_equal(nEdges(phy4), 19)
})

test_that("nexml to phylo4d", {
    nxml <- RNeXML::nexml_read(compFile)
    phy4d <- as(nxml, "phylo4d")
    nxmldt <- RNeXML::get_characters(nxml)
    phy4d2 <- phylo4d(get_trees(nxml), nxmldt[sample(1:nrow(nxmldt)), ])
    expect_true(all(tipLabels(phy4d) %in% paste("taxon", 1:10, sep="_")))
    expect_equal(nEdges(phy4d), 19)
    expect_equal(phy4d, phy4d2)
    expect_equal(ncol(tdata(phy4d, "tip")), 2)
    expect_true(all(names(tdata(phy4d, "tip")) %in% c("log.snout.vent.length", "reef.dwelling")))
})

test_that("phylo4 to phylo", {
  ## phylo tree in unknown order
  expect_equal(suppressWarnings(as(phy, "phylo")), tr)
  # ...now check for warning for unknown order
  expect_warning(as(phy, "phylo"))

  # phylo tree in cladewise order
  tr.cladewise <- reorder(tr, "cladewise")
  phy.c <- as(tr.cladewise, "phylo4")
  expect_equal(as(phy.c, "phylo"), tr.cladewise)

  # phylo tree in pruningwise order
  tr.pruningwise <- reorder(tr, "pruningwise")
  phy.p <- as(tr.pruningwise, "phylo4")
  expect_equal(suppressWarnings(as(phy.p, "phylo")), tr.pruningwise)

  # after transforming the jumbled tree to phylo and back, edge matrix
  # and edge slots should still be in the original order, but node slots
  # should be back in nodeId order
  phy.r <- reorder(phy.alt)
  phy.roundtrip.r <- reorder(as(suppressWarnings(as(phy.alt, "phylo")), "phylo4"))
  expect_equal(edges(phy.roundtrip.r), edges(phy.r))
  expect_equal(edgeLength(phy.roundtrip.r), edgeLength(phy.r))
  expect_equal(labels(phy.roundtrip.r), labels(phy.r))
})

## this coerce method is defined implicitly
test_that("phylo to phylo4d", {
  ## phylo tree in unknown order
  phyd <- as(tr, "phylo4d")
  tipData(phyd) <- data.frame(x=1:5, row.names=tipLabels(phyd))
  expect_equal(suppressWarnings(as(phyd, "phylo")), tr)
  ## ...now check for warning for unknown order 
  expect_warning(as(phyd, "phylo"))

  ## phylo tree in cladewise order
  tr.cladewise <- reorder(tr, "cladewise")
  phyd <- as(tr.cladewise, "phylo4d")
  tipData(phyd) <- data.frame(x=1:5, row.names=tipLabels(phyd))
  expect_equal(suppressWarnings(as(phyd, "phylo")), tr.cladewise)
  ## ...now check for warning for dropping data
  expect_warning(as(phyd, "phylo"))

  ## phylo tree in pruningwise order
  tr.pruningwise <- reorder(tr, "pruningwise")
  phyd <- as(tr.pruningwise, "phylo4d")
  tipData(phyd) <- data.frame(x=1:5, row.names=tipLabels(phyd))
  expect_equal(suppressWarnings(as(phyd, "phylo")), tr.pruningwise)
})

##test.phylo4.As.phylog <- function() {
##}

test_that("phylo4 to data.frame", {
  phy.show <- phylobase:::.phylo4ToDataFrame(phy.alt, "pretty")
  expect_equal(phy.show$label, c(lab.tip, lab.int))
  expect_equal(phy.show$node, c(nid.tip, nid.int))
  expect_equal(phy.show$ancestor, ancestor[match(c(nid.tip, nid.int),
    descendant)])
  expect_equal(phy.show$edge.length, sort(elen))
  expect_equal(phy.show$node.type, factor(unname(nodeType(phy))))
})

## core functionality is already tested in test..phylo4ToDataFrame()
test_that("phylo4 to data.frame", {
    ## rooted tree
    expect_true(is.data.frame(as(phy, "data.frame")))

    ## unrooted tree
    tru <- ape::unroot(tr)
    phyu <- as(tru, "phylo4")
    # should probably check that this coercion results in something
    # *correct*, not just that it produces a data.frame
    expect_true(is.data.frame(as(phyu, "data.frame")))
})
