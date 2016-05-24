##
## --- Test methods-phylo4.R ---
##

# create ape::phylo version of a simple tree for testing
nwk <- "((t1:0.1,t2:0.2)n7:0.7,(t3:0.3,(t4:0.4,t5:0.5)n9:0.9)n8:0.8)n6:0.6;"
tr <- read.tree(text=nwk)

# create analogous phylo4 object with a full complement of valid slots
ancestor <- as.integer(c(6,7,7,6,8,0,8,9,9))
descendant <- as.integer(c(7,1,2,8,3,6,9,4,5))
edge <- cbind(ancestor, descendant)
nid.tip <- 1:5
nid.int <- 6:9
nid.all <- c(nid.tip, nid.int)
lab.tip <- paste("t", nid.tip, sep="")
lab.int <- paste("n", nid.int, sep="")
lab.all <- c(lab.tip, lab.int)
eid <- paste(ancestor, descendant, sep="-")
elen <- descendant/10
elab <- paste("e", eid, sep="")
phy <- phylo4(x=edge, tip.label=lab.tip, node.label=lab.int,
    edge.length=elen, edge.label=elab)

# create altered version such that each slot is out of order with
# respect to all others; methods should be able to handle this
phy.alt <- phy
phy.alt@label <- rev(phy@label)
phy.alt@edge <- phy@edge[c(6:9, 1:5), ]
phy.alt@edge.length <- phy@edge.length[c(7:9, 1:6)]
phy.alt@edge.label <- phy@edge.label[c(8:9, 1:7)]

# update test targets for edge-related slots
ancestor <- ancestor[c(6:9, 1:5)]
descendant <- descendant[c(6:9, 1:5)]
edge <- cbind(ancestor, descendant)
eid <- eid[c(6:9, 1:5)]
elen <- elen[c(6:9, 1:5)]
elab <- elab[c(6:9, 1:5)]

op <- phylobase.options()
#-----------------------------------------------------------------------

context("nTips, depthTips, nNodes, nodeType")

test_that("nTips works correctly",
  expect_that(nTips(phy.alt), equals(length(nid.tip)))
)

test_that("depthTips works when there are edge lengths", {
    edgeLengthVec <- c(1.2, 1.8, 1.8, 2.1, 2.3)
    names(edgeLengthVec) <- tipLabels(phy.alt)
    expect_that(depthTips(phy.alt), equals(edgeLengthVec))
})

test_that("depthTips works when there are no edge lengths", {
    tmpPhy <- phy.alt
    edgeLength(tmpPhy) <- NA
    expect_true(is.null(depthTips(tmpPhy)))
})

test_that("nTips works on ape objects",
          ## nTips phylo
          expect_equal(nTips(tr), 5))

test.nEdges.phylo4 <- function() {
  expect_identical(nEdges(phy.alt), nrow(edge))
}

test_that("nNodes works as expected",
          expect_equal(nNodes(phy.alt), length(nid.int)))

test_that("nodeType works as expected",
          expect_identical(nodeType(phy.alt),
                           setNames(c(rep("tip", length(nid.tip)),
                                      "root",
                                      rep("internal", length(nid.int)-1)),
                                    c(nid.tip, nid.int))))

context("nodeId")
test_that("nodeId works without arguments",
          expect_identical(nodeId(phy.alt), c(nid.tip, nid.int)))
test_that("nodeId works with argument all",
          expect_identical(nodeId(phy.alt, "all"), c(nid.tip, nid.int)))
test_that("nodeId works with argument tip",
          expect_identical(nodeId(phy.alt, "tip"), nid.tip))
test_that("nodeId works with argument internal",
          expect_identical(nodeId(phy.alt, "internal"), nid.int))
test_that("nodeId works woth argument root",
          expect_identical(nodeId(phy.alt, "root"), nid.int[1]))


context("nodeDepth")
allDepths <- c(1.2, 1.8, 1.8, 2.1, 2.3, 0.9, 1.0, 1.2, 1.6)
names(allDepths) <- names(getNode(phy.alt))
test_that("nodeDepth works without arguments", {
    expect_equal(nodeDepth(phy.alt), allDepths)
})

test_that("nodeDepth works with numeric argument", {
    expect_equal(nodeDepth(phy.alt, 1), allDepths[1])
})

test_that("nodeDepth works with character argument", {
    expect_equal(nodeDepth(phy.alt, "t1"), allDepths[1])
})

test_that("nodeDepth works with no branch length", {
    tmpPhy <- phy.alt
    edgeLength(tmpPhy) <- NA
    expect_true(is.null(nodeDepth(tmpPhy)))
})

############################################################################
## nodeHeight                                                             ##
############################################################################

context("nodeHeight")

tmp_nd_hgt_tree <- tempfile()
cat("(((A:1,B:1):2,(C:1,D:1):2):4,((E:10,F:1):2,(G:3,H:7):2):4);",
    file = tmp_nd_hgt_tree)
nd_hgt_tree <- readNewick(file = tmp_nd_hgt_tree)
unlink(tmp_nd_hgt_tree)

test_that("nodeHeight with 1 node", {
              expect_equal(nodeHeight(nd_hgt_tree, MRCA(nd_hgt_tree, c("A", "D")), "all_tip"),
                           setNames(c(3, 3, 3, 3), c("A", "B", "C", "D")))
              expect_equal(nodeHeight(nd_hgt_tree, MRCA(nd_hgt_tree, c("E", "H")), "min_tip"),
                           c("F" = 3))
              expect_equal(nodeHeight(nd_hgt_tree, MRCA(nd_hgt_tree, c("E", "H")), "max_tip"),
                           c("E" = 12))
              expect_equal(nodeHeight(nd_hgt_tree, MRCA(nd_hgt_tree, c("A", "D")), "root"),
                           4)
          })

test_that("nodeHeight with several nodes", {
              expect_equal(nodeHeight(nd_hgt_tree, c(
                  MRCA(nd_hgt_tree, c("A", "D")),
                  MRCA(nd_hgt_tree, c("A", "B"))),
                                      "all_tip"),
                           list("10" = setNames(c(3, 3, 3, 3), c("A", "B", "C", "D")),
                                "11" = c("A" = 1, "B" = 1)))

              expect_equal(nodeHeight(nd_hgt_tree, c(
                  MRCA(nd_hgt_tree, c("E", "H")),
                  MRCA(nd_hgt_tree, c("E", "F"))),
                                      "min_tip"),
                           list("13" = c("F" = 3),
                                "14" = c("F" = 1)))

              expect_equal(nodeHeight(nd_hgt_tree, c(
                  MRCA(nd_hgt_tree, c("E", "H")),
                  MRCA(nd_hgt_tree, c("E", "F"))), "max_tip"),
                           list("13" = c("E" = 12),
                                "14" = c("E" = 10)))

              expect_equal(nodeHeight(nd_hgt_tree, c(
                  MRCA(nd_hgt_tree, c("A", "D")),
                  MRCA(nd_hgt_tree, c("E", "F"))),
                                      "root"),
                           c("10" = 4, "14" = 6))
          })


test_that("nodeHeight for tips", {
              res <- as.list(rep(0, nTips(nd_hgt_tree)))
              for (i in seq_len(nTips(nd_hgt_tree))) names(res[[i]]) <- LETTERS[i]
              names(res) <- seq_len(nTips(nd_hgt_tree))

              expect_equal(nodeHeight(nd_hgt_tree, nodeId(nd_hgt_tree, "tip"), "all_tip"),
                           res)
              expect_equal(nodeHeight(nd_hgt_tree, nodeId(nd_hgt_tree, "tip"), "min_tip"),
                           res)
              expect_equal(nodeHeight(nd_hgt_tree, nodeId(nd_hgt_tree, "tip"), "max_tip"),
                           res)
          })

test_that("nodeHeight for mix of tips and internal nodes", {
              expect_equal(nodeHeight(nd_hgt_tree, c(1, 10), "all_tip"),
                           list("1" = c("A" = 0),
                                "10" = c("A" = 3, "B" = 3, "C" = 3, "D" = 3)))
              expect_equal(nodeHeight(nd_hgt_tree, c(1, 14), "min_tip"),
                           list("1" = c("A" = 0),
                                "14" = c("F" = 1)))
              expect_equal(nodeHeight(nd_hgt_tree, c(1, 14), "max_tip"),
                           list("1" = c("A" = 0),
                                "14" = c("E" = 10)))
              expect_equal(nodeHeight(nd_hgt_tree, c(5, 14), "root"),
                           c("5" = 16, "14" = 6))
          })


############################################################################
## edges                                                                  ##
############################################################################

context("edges")
test_that("edges works",  expect_identical(edges(phy.alt), edge))
test_that("edges work with drop.root=TRUE option",
          expect_identical(edges(phy.alt, drop.root=TRUE),
                           edge[edge[,1] != 0,]))

context("edge order")
test_that("edgeOrder works as expected", {
    expect_identical(edgeOrder(phy.alt), "unknown")
    expect_identical(edgeOrder(reorder(phy.alt, "preorder")), "preorder")
    expect_identical(edgeOrder(reorder(phy.alt, "postorder")), "postorder")
})

context("edgeId")
test_that("edgeId works with no argument",
          expect_identical(edgeId(phy.alt), eid))
test_that("edgeId works with argument all",
          expect_identical(edgeId(phy.alt, "all"), eid))
test_that("edgeId works with argument tip",
          expect_identical(edgeId(phy.alt, "tip"), eid[descendant %in% nid.tip]))
test_that("edgeId works with argument internal",
  expect_identical(edgeId(phy.alt, "internal"), eid[!descendant %in% nid.tip]))
test_that("edgeId works with argument root",
          expect_identical(edgeId(phy.alt, "root"), eid[ancestor == 0]))

context("hasEdgeLength")
test_that("hasEdgeLength works when edge lengths are present",
          expect_true(hasEdgeLength(phy.alt)))
test_that("hasEdgeLength works when no edge lengths are present", {
    phy.alt@edge.length <- NA_real_
    expect_true(!hasEdgeLength(phy.alt))
})


context("edgeLength")
test_that("default works (all edge lengths)",
          expect_identical(edgeLength(phy.alt), setNames(elen, eid)))
test_that("one edge length, by label",
    expect_equal(edgeLength(phy.alt, "t1"), c(`7-1`=0.1)))
test_that("one edge length, by node ID",
          expect_equal(edgeLength(phy.alt, 1), c(`7-1`=0.1)))
test_that("non-existent edge, by label", {
    ans <- structure(NA_real_, .Names = NA_character_)
    expect_equal(suppressWarnings(edgeLength(phy.alt, "xxx")), ans)
})
test_that("non-existent edge, by number", {
          ans <- structure(NA_real_, .Names = NA_character_)
          expect_equal(suppressWarnings(edgeLength(phy.alt, 999)), ans)
})
test_that("wrong number of edge lengths", {
    phy.tmp1 <- phy.alt
    phy.tmp1@edge.length <- phy.alt@edge.length[-1]
    expect_true(nzchar(checkPhylo4(phy.tmp1)))
    phy.tmp1 <- phy.alt
    phy.tmp1@edge.length <- c(phy.alt@edge.length, 1)
    expect_true(nzchar(checkPhylo4(phy.tmp1)))
})
test_that("negative edge lengths", {
    phy.tmp1 <- phy.alt
    phy.tmp1@edge.length[3] <- -1
    expect_true(nzchar(checkPhylo4(phy.tmp1)))
})
test_that("edge incorrectly labeled", {
    phy.tmp1 <- phy.alt
    names(phy.tmp1@edge.length)[1] <- "9-10"
    expect_true(nzchar(checkPhylo4(phy.tmp1)))
})

context("edgeLength <-")
emptyVec <- numeric()
attributes(emptyVec) <- list(names=character(0))
test_that("dropping all should produce empty slot", {
          edgeLength(phy.alt) <- numeric()
          expect_identical(edgeLength(phy.alt), setNames(rep(NA_real_, 9), edgeId(phy.alt, "all")))
          expect_identical(phy.alt@edge.length, emptyVec)
          edgeLength(phy.alt) <- NA_real_
          expect_identical(edgeLength(phy.alt), setNames(rep(NA_real_, 9), edgeId(phy.alt, "all")))
          expect_identical(phy.alt@edge.length, emptyVec)
})
test_that("vector with reversed names, get matched by default (complete replacement)", {
    edgeLength(phy.alt) <- numeric()
    revElen <- setNames(elen, rev(eid))
    edgeLength(phy.alt) <- revElen
    expect_identical(edgeLength(phy.alt), revElen[edgeId(phy.alt, "all")])
})
test_that("vector with reversed names, but specify no matching (complete replacement)", {
    edgeLength(phy.alt) <- numeric()
    revElen <- setNames(elen, rev(eid))
    edgeLength(phy.alt, use.names=FALSE) <- revElen
    elen1 <- elen
    expect_identical(edgeLength(phy.alt), setNames(elen1, edgeId(phy.alt, "all")))
})
test_that("vector with no names, should match to edgeId order (complete replacement)", {
    edgeLength(phy.alt) <- numeric()
    edgeLength(phy.alt) <- elen
    elen2 <- elen
    expect_identical(edgeLength(phy.alt), setNames(elen2, edgeId(phy.alt, "all")))
})
test_that("recycling applies if fewer the nEdges elements are supplied, \
          (duplicate edge length are okay), (complete replacement)", {
              edgeLength(phy.alt) <- 1
              expect_identical(edgeLength(phy.alt), setNames(rep(1, 9), edgeId(phy.alt, "all")))
})
edgeLength(phy.alt) <- elen
test_that("replace an edge length using numeric index (partial replacement)", {
          edgeLength(phy.alt)[9] <- 83
          expect_identical(edgeLength(phy.alt), setNames(c(elen[1:8], 83), edgeId(phy.alt, "all")))
})
test_that("back again, now using character index (partial replacement)", {
    edgeLength(phy.alt)["8-3"] <- 0.3
    elen3 <- elen
    expect_identical(edgeLength(phy.alt), setNames(elen3, edgeId(phy.alt, "all")))
})
test_that("error to add length for edges that don't exist (partial replacement)", {
    expect_error(edgeLength(phy.alt)["fake"] <- 999)
    expect_error(edgeLength(phy.alt)[999] <- 999)
})
test_that("NAs permitted only for root edge (or for *all* edges)", {
    edgeLength(phy.alt)[edgeId(phy.alt, "root")] <- NA
    expect_identical(edgeLength(phy.alt), setNames(c(NA, elen[2:9]), edgeId(phy.alt, "all")))
    edgeLength(phy.alt) <- elen
    expect_error(edgeLength(phy.alt)["8-3"] <- NA)
})


## TODO sumEdgeLength.phylo4 ## function(phy, node)

context("isRooted")
test_that("isRooted works as expected",
          expect_true(isRooted(phy.alt)))

context("rootNode")
test_that("rootNode works as expected",
          expect_identical(rootNode(phy.alt), getNode(phy, nid.int[1])))

context("rootNode <-")
test_that("rootNode <- is not yet implemented",
          expect_error(rootNode(phy.alt) <- 7))

context("labels")
test_that("labels works as expected with no argument",
  expect_identical(labels(phy.alt),
                   setNames(c(lab.tip, lab.int), c(nid.tip, nid.int))))
test_that("labels works as expected with argument all",
          expect_identical(labels(phy.alt, "all"),
                           setNames(c(lab.tip, lab.int), c(nid.tip, nid.int))))
test_that("labels works as expected with argument tip",
          expect_identical(labels(phy.alt, "tip"), setNames(lab.tip, nid.tip)))
test_that("labels works as expected with argument internal",
          expect_identical(labels(phy.alt, "internal"), setNames(lab.int, nid.int)))


context("labels <-")
test_that("dropping all should produce default tip labels, no internal labels", {
    labels(phy.alt) <- character()
    expect_identical(labels(phy.alt),
                     setNames(c(paste("T", 1:5, sep=""), rep(NA, 4)), nid.all))
})

## #
## # complete replacement
## #

## with names, not used
test_that("vector with reversed names, but names not used (all) - complete replacement", {
    labels(phy.alt) <- character()
    labels(phy.alt) <- setNames(lab.all, rev(nid.all))
    expect_identical(labels(phy.alt), setNames(lab.all, nid.all))
})
test_that("vector with reversed names, but names not used (tips) - complete replacement", {
    labels(phy.alt) <- character()
    labels(phy.alt, "tip") <- setNames(lab.tip, rev(nid.tip))
    expect_identical(tipLabels(phy.alt), setNames(lab.tip, nid.tip))
})
test_that("vector with reversed names, but names not used (internal) - complete replacement", {
    labels(phy.alt) <- character()
    labels(phy.alt, "internal") <- setNames(lab.int, rev(nid.int))
    expect_identical(nodeLabels(phy.alt), setNames(lab.int, nid.int))
})

## with names, used
test_that("vector with reversed names, but names used (all) - complete replacement", {
  labels(phy.alt) <- character()
  labels(phy.alt, use.names=TRUE) <- setNames(lab.all, rev(nid.all))
  expect_identical(labels(phy.alt), setNames(rev(lab.all), nid.all))
})
test_that("vector with reversed names, but names used (tips) - complete replacement", {
    labels(phy.alt) <- character()
    labels(phy.alt, "tip", use.names=TRUE) <- setNames(lab.tip, rev(nid.tip))
    expect_identical(tipLabels(phy.alt), setNames(rev(lab.tip), nid.tip))
})
test_that("vector with reversed names, but names used (internal) - complete replacement", {
    labels(phy.alt) <- character()
    labels(phy.alt, "internal", use.names=TRUE) <- setNames(lab.int, rev(nid.int))
    expect_identical(nodeLabels(phy.alt), setNames(rev(lab.int), nid.int))
})
## no names
test_that("vector with no names, should match to nodeId order (all) - complete replacement", {
  labels(phy.alt) <- character()
  labels(phy.alt) <- lab.all
  expect_identical(labels(phy.alt), setNames(lab.all, nid.all))
})
test_that("vector with no names, should match to nodeId order (all) - complete replacement", {
  labels(phy.alt) <- character()
  labels(phy.alt, type="tip") <- lab.tip
  expect_identical(tipLabels(phy.alt), setNames(lab.tip, nid.tip))
})
test_that("vector with no names, should match to nodeId order (all) - complete replacement", {
  labels(phy.alt) <- character()
  labels(phy.alt, type="internal") <- lab.int
  expect_identical(nodeLabels(phy.alt), setNames(lab.int, nid.int))
})

## partial replacement
labels(phy.alt) <- lab.all
test_that("replace a tip using numeric index", {
    labels(phy.alt)[5] <- "t5a"
    expect_identical(tipLabels(phy.alt), setNames(c(lab.tip[1:4], "t5a"), nid.tip))
})
test_that("and back again, now using character index", {
    labels(phy.alt)["5"] <- "t5"
    expect_identical(labels(phy.alt), setNames(lab.all, nid.all))
})
test_that("replace an internal node using numeric index", {
    labels(phy.alt)[9] <- "n9a"
    expect_identical(nodeLabels(phy.alt), setNames(c(lab.int[1:3], "n9a"), nid.int))
})
test_that("and back again, now using character index", {
    labels(phy.alt)["9"] <- "n9"
    expect_identical(labels(phy.alt), setNames(lab.all, nid.all))
})
test_that("error to produce duplicate tip or internal label", {
    phylobase.options(allow.duplicated.labels="fail")
    expect_error(labels(phy.alt)[1] <- "t2")
    expect_error(labels(phy.alt)[6] <- "n7")
})
test_that("no error in allow.duplicated.labels is ok", {
    phylobase.options(allow.duplicated.labels="ok")
    labels(phy.alt)[1] <- "t2"
    labels(phy.alt)[6] <- "n7"
    expect_identical(tipLabels(phy.alt), setNames(c("t2", "t2", "t3", "t4", "t5"), nid.tip))
    expect_identical(nodeLabels(phy.alt), setNames(c("n7", "n7", "n8", "n9"), nid.int))
})
test_that("error to add labels for nodes that don't exist", {
    expect_error(labels(phy.alt)["fake"] <- "xxx")
    expect_error(labels(phy.alt)[999] <- "xxx")
})

context("nodeLabels")
test_that("nodeLabels works as expected",
          expect_identical(nodeLabels(phy.alt), setNames(lab.int, nid.int)))

context("hasNodeLabels")
test_that("hasNodeLabels works as expected", {
    expect_true(hasNodeLabels(phy.alt))
    nodeLabels(phy.alt) <- NA_character_
    expect_true(!hasNodeLabels(phy.alt))
})

context("nodeLabels <-")
test_that("dropping all should produce no internal labels", {
    nodeLabels(phy.alt) <- character()
    expect_true(!any(nid.int %in% names(phy.alt@label)))
    expect_identical(nodeLabels(phy.alt), setNames(rep(NA_character_, 4), nid.int))
})
labels(phy.alt) <- lab.all
test_that("replace an internal node using numeric index", {
    nodeLabels(phy.alt)[4] <- "n9a"
    expect_identical(nodeLabels(phy.alt), setNames(c(lab.int[1:3], "n9a"), nid.int))
})
test_that("and back again, now using character index", {
    nodeLabels(phy.alt)["9"] <- "n9"
    expect_identical(labels(phy.alt), setNames(lab.all, nid.all))
})
test_that("error to produce duplicate internal label", {
    phylobase.options(allow.duplicated.labels="fail")
    expect_error(nodeLabels(phy.alt)["6"] <- "n7")
})
test_that("duplicated labels work as expected", {
    phylobase.options(op)
    phylobase.options(allow.duplicated.labels="ok")
    nodeLabels(phy.alt)["6"] <- "n7"
    expect_identical(nodeLabels(phy.alt), setNames(c("n7", "n7", "n8", "n9"), nid.int))
    expect_true(hasDuplicatedLabels(phy.alt))
    ## NAs are not considered duplicated
    nodeLabels(phy.alt)[1:2] <- NA
    expect_true(!hasDuplicatedLabels(phy.alt))
    phylobase.options(op)
    ## error to add labels for nodes that don't exist
    expect_error(nodeLabels(phy.alt)["fake"] <- "xxx")
    expect_error(nodeLabels(phy.alt)[999] <- "xxx")
})

context("tipLabels")
test_that("tipLabels works as expected",
          expect_identical(tipLabels(phy.alt), setNames(lab.tip, nid.tip)))

context("tipLabels <-")
test_that("dropping all tip labels should produce default labels", {
    tipLabels(phy.alt) <- character()
    expect_identical(tipLabels(phy.alt), setNames(paste("T", 1:5, sep=""), nid.tip))
})
labels(phy.alt) <- lab.all
test_that("replace a tip using numeric index", {
    tipLabels(phy.alt)[5] <- "t5a"
    expect_identical(tipLabels(phy.alt), setNames(c(lab.tip[1:4], "t5a"), nid.tip))
})
test_that("and back again, now using character index", {
    tipLabels(phy.alt)["5"] <- "t5"
    expect_identical(labels(phy.alt), setNames(lab.all, nid.all))
})
test_that("error to produce duplicate tip or internal label", {
    phylobase.options(allow.duplicated.labels="fail")
    expect_error(tipLabels(phy.alt)[1] <- "t2")
})
test_that("duplicated labels works as expected on tips", {
    phylobase.options(op)
    phylobase.options(allow.duplicated.labels="ok")
    tipLabels(phy.alt)[1] <- "t2"
    expect_identical(tipLabels(phy.alt), setNames(c("t2", "t2", "t3", "t4", "t5"), nid.tip))
    expect_true(hasDuplicatedLabels(phy.alt))
    tipLabels(phy.alt)[1:2] <- NA
    expect_true(!hasDuplicatedLabels(phy.alt))
    phylobase.options(op)
})
test_that("error to add labels for nodes that don't exist", {
    expect_error(tipLabels(phy.alt)["fake"] <- "xxx")
    expect_error(tipLabels(phy.alt)[999] <- "xxx")
})
test_that("hasEdgeLabels works as expected", {
    expect_true(hasEdgeLabels(phy.alt))
    phy.alt@edge.label <- NA_character_
    expect_true(!hasEdgeLabels(phy.alt))
})

context("edgeLabels")
test_that("edgeLabels works as expected", {
    expect_identical(edgeLabels(phy.alt), setNames(elab, eid))
})
test_that("edgeLabels returns named vector of NAs if edge labels are missing or NA", {
    phy.alt@edge.label <- NA_character_
    expect_identical(edgeLabels(phy.alt), setNames(rep(NA_character_, 9), eid))
    phy.alt@edge.label <- character()
    expect_identical(edgeLabels(phy.alt), setNames(rep(NA_character_, 9), eid))
})
test_that("if only some labels exists, should fill in NA for the others", {
    phy.alt@edge.label <- setNames(elab[-1], eid[-1])
    expect_identical(edgeLabels(phy.alt), setNames(c(NA, elab[-1]), eid))
})


context("edgeLabels <-")
test_that(" dropping all should produce empty slot", {
    edgeLabels(phy.alt) <- character()
    expect_identical(edgeLabels(phy.alt), setNames(rep(NA_character_, 9), eid))
})
test_that("vector with reversed names, which always get matched - complete replacement", {
    edgeLabels(phy.alt) <- character()
    edgeLabels(phy.alt) <- setNames(elab, rev(eid))
    expect_identical(edgeLabels(phy.alt), setNames(rev(elab), eid))
})
test_that("vector with no names, should match to edgeId order - complete replacement", {
    edgeLabels(phy.alt) <- character()
    edgeLabels(phy.alt) <- elab
    expect_identical(edgeLabels(phy.alt), setNames(elab, eid))
})
test_that("recycling applies if fewer the nEdges elements are supplied\\
           (duplicate edge labels are okay) - complete replacement.", {
               edgeLabels(phy.alt) <- "x"
               expect_identical(edgeLabels(phy.alt), setNames(rep("x", 9), eid))
           })
edgeLabels(phy.alt) <- elab
test_that("replace an edge label using numeric index - partial replacement", {
  edgeLabels(phy.alt)[9] <- "e8-3a"
  expect_identical(edgeLabels(phy.alt), setNames(c(elab[1:8], "e8-3a"), eid))
})
test_that("and back again, now using character index", {
    edgeLabels(phy.alt)["8-3"] <- "e8-3"
    expect_identical(edgeLabels(phy.alt), setNames(elab, eid))
})
test_that("error to add labels for edges that don't exist", {
    expect_error(edgeLabels(phy.alt)["fake"] <- "xxx")
    expect_error(edgeLabels(phy.alt)[999] <- "xxx")
})

## this is also the print method
## this mostly just wraps .phylo4ToDataFrame, which is tested elsewhere
## test.show.phylo4 <- function() {
## }
## test.names.phylo4 <- function() {
##   #TODO?
## }
## test.head.phylo4 <- function() {
##   #TODO?
## }
## test.tail.phylo4 <- function() {
##   #TODO?
## }

context("summary")
test_that("summary works as expected", {
    phy.sum <- summary(phy.alt, quiet=TRUE)
    expect_identical(phy.sum$name, "phy.alt")
    expect_identical(phy.sum$nb.tips, length(nid.tip))
    expect_identical(phy.sum$nb.nodes, length(nid.int))
    expect_identical(phy.sum$mean.el, mean(elen))
    expect_identical(phy.sum$var.el, var(elen))
    expect_identical(phy.sum$sumry.el, summary(elen))
})
test_that("summary works as expected when root edge as no length", {
    ## now make root edge length NA
    edgeLength(phy.alt)[edgeId(phy.alt, "root")] <- NA
    phy.sum2 <- summary(phy.alt, quiet=TRUE)
    expect_identical(phy.sum2$mean.el, mean(edgeLength(phy.alt), na.rm=TRUE))
    expect_identical(phy.sum2$var.el, var(edgeLength(phy.alt), na.rm=TRUE))
    expect_identical(phy.sum2$sumry.el, summary(stats::na.omit(edgeLength(phy.alt))))
})
test_that("now remove edge lengths altogether", {
    phy.alt@edge.length[] <- NA
    phy.sum3 <- summary(phy.alt, quiet=TRUE)
    expect_true(is.null(phy.sum3$mean.el))
    expect_true(is.null(phy.sum3$var.el))
    expect_true(is.null(phy.sum3$sumry.el))
})

## not an exported function -- called internally by reorder("phylo4")
## test.orderIndex <- function() {
## }

## test.reorder.phylo4 <- function() {
##   ## TODO
## }

context("isUltrametric")
test_that("isUltrametric works as expected", {
    expect_true(!isUltrametric(phy.alt))
    tmpPhy <- as(rcoal(10), "phylo4")
    expect_true(isUltrametric(tmpPhy))
    tmpPhy <- phy.alt
    edgeLength(tmpPhy) <- NA
    expect_error(isUltrametric(tmpPhy))
})

phylobase.options(op)
