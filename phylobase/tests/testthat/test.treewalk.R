#
# --- Test treewalk.R ---
#

# Create sample phylo4 tree for testing
tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;")
phytr <- as(tr, "phylo4")

# create phylo4 object with a full complement of valid slots
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

#-----------------------------------------------------------------------



## Note: we're not explicitly testing missing="warn" condition below;
## however, if "OK" and "fail" both work as expected, then so must "warn"

#test.getNode <- function() {

context("getNode")
test_that("getNode works when nodes provided only has valid characters", {
    expect_that(getNode(phytr, "spA"), equals(c(spA=1)))
    expect_that(getNode(phytr, c("spA", "spC")), equals(c(spA=1, spC=3)))
})

test_that("getNode works when nodes provided only has valid integers", {
    ans <- 4
    names(ans) <- "spD"
    expect_that(getNode(phytr, 4), equals(ans))
    ans <- c(4,6)
    names(ans) <- c("spD", NA)
    expect_that(getNode(phytr, c(4,6)), equals(ans))
})

test_that("getNode works when node includes only missing characters (names), but missing=OK", {
    ans <- rep(NA_integer_, 2)  # return values should be NA
    names(ans) <- rep(NA, 2)  # return values should have NA names
    expect_that(getNode(phytr, c("xxx", "yyy"), missing="OK"), equals(ans))
    # now missing = "fail"
    expect_error(getNode(phytr, c("xxx", "yyy"), missing="fail"))
})

test_that("getNode works wehn node includes only missing numbers (IDs), but missing=OK", {
    ans <- rep(NA_integer_, 3)  # return values should be NA
    names(ans) <- rep(NA, 3)  # return values should have NA names
    expect_that(getNode(phytr, c(-9, 0, 50), missing="OK"), equals(ans))
    # now missing = "fail"
    expect_error(getNode(phytr, c(-9, 0, 50), missing="fail"))
})

test_that("getNode works when node includes NAs, but missing = \"OK\"", {
    expect_true(is.na(getNode(phytr, NA_integer_, missing="OK")))
    expect_true(is.na(getNode(phytr, NA_character_, missing="OK")))
})

test_that("getNode works when node includes mixture of valid values and NAs", {
    ans <- c(2, NA)
    names(ans) <- c("spB", NA)
    expect_that(getNode(phytr, c("spB", NA), missing="OK"), equals(ans))
    expect_that(getNode(phytr, c(2, NA), missing="OK"), equals(ans))
})

test_that("getNode throws exception when node is neither integer-like nor character",
    expect_error(getNode(phytr, 1.5)))

test_that("getNode works even when a tip is labeled as \"0\"", {
    phyTmp <- phytr
    tipLabels(phyTmp)[1] <- "0"
    ans <- 1
    names(ans) <- "0"
    expect_that(getNode(phyTmp, "0"), equals(ans))
})

## TODO context("ancestor function")

## TODO context("children function")


context("descendants function")
phytr <- phylo4(read.tree(text="((t3,t4),(t1,(t2,t5)));"))

test_that("descendants() works with tips", {
    expect_identical(descendants(phytr, 5), setNames(5L, "t5"))
    expect_identical(descendants(phytr, 5, "tips"), setNames(5L, "t5"))
    expect_identical(descendants(phytr, 5, "children"),
                     setNames(integer(0), character(0)))
    expect_identical(descendants(phytr, 5, "all"), setNames(5L, "t5"))
    expect_identical(descendants(phytr, 5, "ALL"), setNames(5L, "t5"))
})

test_that("descendants() works when provided with a vector of nodes", {
              expect_identical(descendants(phytr, 5:7),
                               list("5" = c(t5 = 5L),
                                    "6" = c(t3 = 1L, t4 = 2L, t1 = 3L, t2 = 4L, t5 = 5L),
                                    "7" = c(t3 = 1L, t4 = 2L)))
              expect_identical(descendants(phytr, 5:7, "tips"),
                               list("5" = c(t5 = 5L),
                                    "6" = c(t3 = 1L, t4 = 2L, t1 = 3L, t2 = 4L, t5 = 5L),
                                    "7" = c(t3 = 1L, t4 = 2L)))
              expect_identical(descendants(phytr, 5:7, "children"),
                               list("5" = setNames(integer(0), character(0)),
                                    "6" = setNames(c(7L, 8L), c(NA, NA)),
                                    "7" = c(t3 = 1L, t4 = 2L))
                               )
              expect_identical(descendants(phytr, 5:7, "ALL"),
                               list("5" = c(t5 = 5L),
                                    "6" = setNames(c(6L, 7L, 1L, 2L, 8L, 3L, 9L, 4L, 5L),
                                        c(NA, NA, "t3", "t4", NA, "t1", NA, "t2", "t5")),
                                    "7" = setNames(c(7L, 1L, 2L), c(NA, "t3", "t4")))
                               )
          })

test_that("descendants() works with internal nodes", {
    expect_identical(descendants(phytr, 8),
        setNames(c(3L, 4L, 5L), c("t1", "t2", "t5")))
    expect_identical(descendants(phytr, 8, "tips"),
        setNames(c(3L, 4L, 5L), c("t1", "t2", "t5")))
    expect_identical(descendants(phytr, 8, "children"),
        setNames(c(3L, 9L), c("t1", NA)))
    expect_identical(descendants(phytr, 8, "all"),
                     setNames(c(3L, 9L, 4L, 5L), c("t1", NA, "t2", "t5")))
    expect_identical(descendants(phytr, 8, "ALL"),
                     setNames(c(8L, 3L, 9L, 4L, 5L),
                              c(NA, "t1", NA, "t2", "t5")))
})

## TODO siblings  # function(phy, node, include.self=FALSE)
## TODO ancestors # function (phy, node, type=c("all","parent","ALL"))
## TODO MRCA    # function(phy, ...)
## TODO shortestPath # function(phy, node1, node2)

context("test on getEdge with nodes as descendants")
## function(phy, node, type=c("descendant", "ancestor"),
##     missing=c("warn", "OK", "fail"))

test_that("getEdge works when node only has valid descendants, as characters", {
    expect_identical(getEdge(phy.alt, "t1"), setNames("7-1", 1))
    expect_identical(getEdge(phy.alt, c("t1", "t3")),
                     setNames(c("7-1", "8-3"), c(1,3)))
})

test_that("getEdge works when node only has valid descendants, as integers", {
    expect_identical(getEdge(phy.alt, 1), setNames("7-1", 1))
    expect_identical(getEdge(phy.alt, c(1,3)),
                     setNames(c("7-1", "8-3"), c(1,3)))
})

test_that("node includes only missing characters (labels), missing=OK", {
    expect_identical(getEdge(phy.alt, c("x", "y", "z"), missing="OK"),
                     setNames(rep(NA, 3), rep(NA, 3)))
})

test_that("node includes only missing characters (labels), missing=fail", {
    expect_error(getEdge(phy.alt, c("x", "y", "z"), missing="fail"))
})

test_that("node includes only missing numbers (IDs), but missing=OK",
    expect_identical(getEdge(phy.alt, c(-9, 0, 50), missing="OK"),
                     setNames(rep(NA, 3), rep(NA, 3))))

test_that("node includes only missing numbers (IDs), but missing=fail",
    expect_error(getEdge(phy, c(-9, 0, 50), missing="fail")))

test_that("node includes NAs, but missing = OK", {
    expect_true(is.na(getEdge(phy, NA_integer_, missing="OK")))
    expect_true(is.na(getEdge(phy, NA_character_, missing="OK")))
})

test_that("node includes mixture of valid values and NAs", {
    expect_identical(getEdge(phy, c("t3", NA), missing="OK"),
                     setNames(c("8-3", NA), c(3, NA)))
    expect_identical(getEdge(phy, c(3, NA), missing="OK"),
                     setNames(c("8-3", NA), c(3, NA)))
})

test_that("node is neither integer-like nor character", {
    expect_error(getEdge(phy, 1.5))
})

context("test on getEdge with nodes as ancestors")

test_that("node only has valid ancestors, as characters", {
    expect_identical(getEdge(phy.alt, "n6", type="ancestor"),
                     setNames(c("6-7", "6-8"), c(6, 6)))
    expect_identical(getEdge(phy.alt, c("n6", "n8"), type="ancestor"),
                     setNames(c("6-7", "6-8", "8-9", "8-3"), c(6, 6, 8, 8)))
})

test_that("node only has valid ancestors, as integers", {
    expect_identical(getEdge(phy.alt, 6, type="ancestor"),
                     setNames(c("6-7", "6-8"), c(6, 6)))
    expect_identical(getEdge(phy.alt, c(6, 8), type="ancestor"),
                     setNames(c("6-7", "6-8", "8-9", "8-3"), c(6, 6, 8, 8)))
    })

test_that("node includes only missing characters (labels), but missing=OK", {
    expect_identical(getEdge(phy.alt, c("x", "y", "z"), type="ancestor",
                             missing="OK"), setNames(rep(NA, 3), rep(NA, 3)))
})

test_that("node includes only tips (labels), but missing=OK", {
    expect_identical(
        getEdge(phy.alt, c("t1", "t3"), type="ancestor", missing="OK"),
        setNames(rep(NA, 2), c(1, 3)))
})

test_that("node includes only tips (labels), now missing = fail", {
    expect_error(getEdge(phy.alt, c("x", "y", "z"), missing="fail"))
    expect_error(getEdge(phy.alt, c("t1", "t3"), type="ancestor",
                         missing="fail"))
})

test_that("node includes only missing numbers (IDs), but missing=OK", {
    expect_identical(
        getEdge(phy.alt, c(-9, 0, 50), type="ancestor", missing="OK"),
        setNames(rep(NA, 3), rep(NA, 3)))
})

test_that("node includes only tips (labels), but missing=OK", {
    expect_identical(
        getEdge(phy.alt, c(1, 3), type="ancestor", missing="OK"),
        setNames(rep(NA, 2), c(1, 3)))
})

test_that("node includes only tips (labels), but missing=fail", {
    expect_error(getEdge(phy.alt, c(-9, 0, 50), missing="fail"))
    expect_error(getEdge(phy.alt, c(1, 3), type="ancestor",
                         missing="fail"))
})

test_that("node includes NAs, but missing = OK", {
    expect_true(is.na(getEdge(phy.alt, NA_integer_, type="ancestor",
                              missing="OK")))
    expect_true(is.na(getEdge(phy.alt, NA_character_, type="ancestor",
                              missing="OK")))
})

test_that("node includes mixture of valid values and NAs", {
    expect_identical(
        getEdge(phy.alt, c("t3", "n8", NA), type="ancestor", missing="OK"),
        setNames(c(NA, "8-9", "8-3", NA), c(3, 8, 8, NA)))
    expect_identical(
        getEdge(phy.alt, c(3, 8, NA), type="ancestor", missing="OK"),
        setNames(c(NA, "8-9", "8-3", NA), c(3, 8, 8, NA)))
})
