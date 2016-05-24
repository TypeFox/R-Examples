#
# --- Test checkdata.R ---
#

if (Sys.getenv("RCMDCHECK") == FALSE) {
    pth <- file.path(getwd(), "..", "inst", "nexusfiles")
} else {
    pth <- system.file(package="phylobase", "nexusfiles")
}
## co1.nex -- typical output from MrBayes. Contains 2 identical trees, the first
## one having posterior probabilities as node labels
co1File <- file.path(pth, "co1.nex")

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
elen <- descendant/10
elab <- paste("e", ancestor, descendant, sep="-")
phy <- phylo4(x=edge, tip.label=lab.tip, node.label=lab.int,
    edge.length=elen, edge.label=elab)

op <- phylobase.options()

context("test phylo4 validator/phylobase.options()")

test_that("test polytomies", {
    phylobase.options(poly="fail")
    expect_error(readNexus(file=co1File, check.node.labels="drop"))
    phylobase.options(op)
})

test_that("test retic", {
    phylobase.options(retic="fail")
    edgeRetic <- rbind(edge, c(6, 3))
    expect_error(phy <- phylo4(x=edgeRetic))
    phylobase.options(op)
})

test_that("test multiroot", {
    phylobase.options(multiroot="fail")
    edgeMultiRoot <- rbind(edge, c(0, 7))
    expect_error(phy <- phylo4(x=edgeMultiRoot))
    phylobase.options(op)
})

test_that("test singleton", {
    phylobase.options(singleton="fail")
    edgeSingleton <- cbind(c(9,7,7,6,6,8,8,10,10,0), 1:10)
    expect_error(phylo4(x=edgeSingleton))
    phylobase.options(op)
})

## checkPhylo4Data <- function() {
## }

## formatData <- function() {
##     # function(phy, dt, type=c("tip", "internal", "all"),
##     #   match.data=TRUE, label.type=c("rownames", "column"),
##     #   label.column=1, missing.data=c("fail", "warn", "OK"),
##     #   extra.data=c("warn", "OK", "fail"), rownamesAsLabels=FALSE)
## }


