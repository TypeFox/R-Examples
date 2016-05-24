#
# --- Test treestruc.R functions ---
#

context("tree structures")

test_that("hasPoly", {
    # construct simple polytomy
    owls <- ape::read.tree(text =
        "((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3);")
    owls$edge <- matrix(c(4,4,4,1,2,3), ncol=2)
    owls$Nnode <- 1
    owls$edge.length <- owls$edge.length[-4]
    tr <- as(owls, "phylo4")
    expect_true(hasPoly(tr))
    # test against empty tree
    expect_true(!hasPoly(new("phylo4")))
})


test_that("hasSingle", {
    # test against empty tree
    expect_true(!hasSingle(new("phylo4")))
})

test_that("hasRetic", {
    # test against empty tree
    expect_true(!hasRetic(new("phylo4")))
})

