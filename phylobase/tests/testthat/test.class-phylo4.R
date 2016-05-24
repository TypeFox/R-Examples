#
# --- Test class-phylo4.R ---
#

### Get all the test files
if (Sys.getenv("RCMDCHECK") == FALSE) {
    pth <- file.path(getwd(), "..", "inst", "nexmlfiles")
} else {
    pth <- system.file(package="phylobase", "nexmlfiles")
}

## NeXML files
compFile <- file.path(pth, "comp_analysis.xml")
stopifnot(file.exists(compFile))

op <- phylobase.options()

context("test phylo4 class")

test_that("building from matrix works", {
    edge <- structure(c(6L, 7L, 8L, 8L, 9L, 9L, 7L, 6L, 7L, 8L, 1L, 9L,
                        2L, 3L, 4L, 5L), .Dim = c(8, 2))
    edge.length <- c(0.2, 0.5, 0.2, 0.15, 0.1, 0.1, 0.7, 1)
    tip.label <- paste("t", 1:5, sep="")
    node.label <- paste("n", 1:4, sep="")
    edge.label <- paste("e", 1:8, sep="")
    order <- "preorder"
    annote <- list(x="annotation")
    phy <- phylo4(edge, edge.length=edge.length, tip.label=tip.label,
                  node.label=node.label, edge.label=edge.label, order=order,
                  annote=annote)

    ## test each slot
    expect_equal(edge, unname(edges(phy)))
    expect_equal(edge.length, unname(edgeLength(phy)))
    expect_equal(4L, nNodes(phy))
    expect_equal(tip.label, unname(tipLabels(phy)))
    expect_equal(node.label, unname(nodeLabels(phy)))
    expect_equal(edge.label, unname(edgeLabels(phy)))
    expect_equal(order, edgeOrder(phy))
    expect_equal(annote, phy@annote)

    ## test improper cases
    ## expect_error(phylo4(edge, edge.length=999))  # recycling is allowed? FM (20140506: yes)
    expect_error(phylo4(edge, tip.label=999))
    expect_error(phylo4(edge, node.label=999))
    ## expect_error(phylo4(edge, edge.label=999))  # recycling is allowed? FM (20140506: yes)
    expect_error(phylo4(edge, order="invalid order"))
    expect_error(phylo4(edge, annote="invalid annotation"))
})


## note: this method mostly just wraps phylo->phylo4 coercion, which is
## tested more thoroughly in runit.setAs-methods.R; focus here is on
## annote and check.node.labels arguments         

test_that("phylo4 can be built from phylo (tests on what's not done in setAs tests)", {
    tr <- ape::read.tree(text="(((t1:0.2,(t2:0.1,t3:0.1):0.15):0.5,t4:0.7):0.2,t5:1):0.4;")
    
    ##
    ## annote
    ##

    annote <- list(x="annotation")
    phy <- phylo4(tr, annote=annote)
    expect_equal(annote, phy@annote)

    ##
    ## check.node.labels
    ##

    # case 0: no node labels
    phy <- phylo4(tr)
    expect_true(!hasNodeLabels(phy))

    # case 1: keep unique character labels
    tr$node.label <- paste("n", 1:4, sep="")
    phy <- phylo4(tr, check.node.labels="keep")
    expect_equal(tr$node.label, unname(nodeLabels(phy)))
    # keeping node labels should be the default
    expect_equal(phy, phylo4(tr))

    # case 2: keep unique number-like character labels
    tr$node.label <- as.character(1:4)
    phy <- phylo4(tr, check.node.labels="keep")
    expect_equal(tr$node.label, unname(nodeLabels(phy)))

    # case 3: keep unique numeric labels, but convert to character
    tr$node.label <- as.numeric(1:4)
    phy <- phylo4(tr, check.node.labels="keep")
    expect_equal(as.character(tr$node.label), unname(nodeLabels(phy)))

    # case 4: must drop non-unique labels
    tr$node.label <- rep("x", 4)
    ## with options allow.duplicated.labels="fail"
    phylobase.options(allow.duplicated.labels="fail")
    expect_error(phylo4(tr))
    expect_error(phylo4(tr, check.node.labels="keep"))
    phylobase.options(op)
    ## test dropping node labels
    phy <- phylo4(tr, check.node.labels="drop")
    expect_true(!hasNodeLabels(phy))
    ## with options allow.duplicated.labels="ok"
    phylobase.options(allow.duplicated.labels="ok")
    phy <- phylo4(tr)
    expect_equal(unname(nodeLabels(phy)), tr$node.label)
    phy <- phylo4(tr, check.node.labels="keep")
    expect_equal(unname(nodeLabels(phy)), tr$node.label)
    phy <- phylo4(tr, check.node.labels="drop")
    expect_true(!hasNodeLabels(phy))
    phylobase.options(op)
})

test_that("nexml to phylo4", {
    nxml <- RNeXML::nexml_read(compFile)
    phy4 <- phylo4(nxml)
    expect_true(all(tipLabels(phy4) %in% paste("taxon", 1:10, sep="_")))
    expect_equal(nEdges(phy4), 19)
})
