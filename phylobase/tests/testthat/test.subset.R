##
## --- Test subset.R ---
##

## create phylo4 object with a full complement of valid slots
ancestor <- as.integer(c(6,7,7,6,8,0,8,9,9))
descendant <- as.integer(c(7,1,2,8,3,6,9,4,5))
edge <- cbind(ancestor, descendant)
nid.tip <- 1:5
nid.int <- 6:9
nid.all <- c(nid.tip, nid.int)
lab.tip <- paste("t", nid.tip, sep="")
lab.int <- paste("n", nid.int, sep="")
elen <- descendant/10
elab <- paste("e", ancestor, descendant, sep="-")
phy <- phylo4(x=edge, tip.label=lab.tip, node.label=lab.int,
                          edge.length=elen, edge.label=elab)

## create altered version such that each slot is out of order with
## respect to all others; methods should be able to handle this
phy.alt <- phy
phy.alt@label <- rev(phy@label)
phy.alt@edge <- phy@edge[c(6:9, 1:5), ]
phy.alt@edge.length <- phy@edge.length[c(7:9, 1:6)]
phy.alt@edge.label <- phy@edge.label[c(8:9, 1:7)]

## now create phylo4d by adding data (with node IDs as row.names)
phyd.alt <- as(phy.alt, "phylo4d")
allDt <- data.frame(a=letters[nid.all], b=10*nid.all, row.names=nid.all)
tdata(phyd.alt, "all") <- allDt

## create altered version such that data slots are out of order with
## respect to all others; methods should be able to handle this
nid.tip.r <- c(2,5,4,3,1)
nid.int.r <- c(8,7,9,6)
nid.all.r <- c(nid.tip.r, nid.int.r)
phyd.alt@data <- phyd.alt@data[rank(nid.all.r), ]

#-----------------------------------------------------------------------

context("subset and friends")

## Also be testing "[" phylo4 methods here
test_that("subset on phylo4", {
    # subset 2 tips
    phy.sub2 <- subset(phy.alt, tips.include=c(2, 5))
    expect_equal(tipLabels(phy.sub2), setNames(c("t2", "t5"), c("1", "2")))
    expect_equal(nodeLabels(phy.sub2), setNames(c("n6"), c("3")))
    expect_equal(edgeLength(phy.sub2),
                 setNames(c(0.6, 0.9, 2.2), c("0-3", "3-1", "3-2")))
    expect_equal(subset(phy.alt, tips.exclude=c(1, 3, 4)), phy.sub2)
    expect_equal(subset(phy.alt, tips.include=c("t2", "t5")), phy.sub2)
    expect_equal(subset(phy.alt, tips.exclude=c("t1", "t3", "t4")), phy.sub2)
    # subset 4 tips
    phy.sub4 <- subset(phy.alt, tips.include=c(1, 2, 4, 5))
    expect_equal(tipLabels(phy.sub4),
                 setNames(c("t1", "t2", "t4", "t5"), c("1", "2", "3", "4")))
    expect_equal(nodeLabels(phy.sub4),
                 setNames(c("n6", "n7", "n9"), c("5", "6", "7")))
    expect_equal(edgeLength(phy.sub4),
                 setNames(c(0.6, 0.4, 0.5, 0.7, 0.1, 0.2, 1.7),
                          c("0-5", "7-3", "7-4", "5-6", "6-1", "6-2", "5-7")))
    expect_equal(subset(phy.alt, tips.exclude=3), phy.sub4)
    expect_equal(subset(phy.alt, tips.include=c("t1", "t2", "t4", "t5")),
                 phy.sub4)
    expect_equal(subset(phy.alt, tips.exclude="t3"), phy.sub4)
    # check variants that should all return the original object
    expect_equal(phy.alt, subset(phy.alt))
    expect_equal(phy.alt, subset(phy.alt, tipLabels(phy.alt)))
    expect_equal(phy.alt, subset(phy.alt, seq_len(nTips(phy.alt))))
    expect_equal(phy.alt, phy.alt[tipLabels(phy.alt)])
    expect_equal(phy.alt, phy.alt[seq_len(nTips(phy.alt))])
    expect_equal(phy.alt, phy.alt[TRUE])
    # error if only one valid tip requested
    expect_error(subset(phy, tips.include="t1"))
    expect_error(suppressWarnings(subset(phy, tips.include=c("t1", "t999"))))
    # error if zero valid tips requested
    expect_error(suppressWarnings(subset(phy, tips.include="t999")))
    # error if more than one subset criteria are supplied
    expect_error(subset(phyd, tips.include="t1", tips.exclude="t3"))
})

## Also testing "[" phylo4d methods here
##TODO get rid of some tests that are pretty much redundant with the
##above, and add tests focused more on tree data
test_that("subset on phylo4d", {
    ## subset 2 tips
    phyd.sub2 <- subset(phyd.alt, tips.include=c(2, 5))
    expect_equal(tipLabels(phyd.sub2), setNames(c("t2", "t5"), c("1", "2")))
    expect_equal(nodeLabels(phyd.sub2), setNames(c("n6"), c("3")))
    expect_equal(edgeLength(phyd.sub2),
                 setNames(c(0.6, 0.9, 2.2), c("0-3", "3-1", "3-2")))
    expect_equal(subset(phyd.alt, tips.exclude=c(1, 3, 4)), phyd.sub2)
    expect_equal(subset(phyd.alt, tips.include=c("t2", "t5")), phyd.sub2)
    expect_equal(subset(phyd.alt, tips.exclude=c("t1", "t3", "t4")), phyd.sub2)
    ## subset 4 tips
    phyd.sub4 <- subset(phyd.alt, tips.include=c(1, 2, 4, 5))
    expect_equal(tipLabels(phyd.sub4),
                 setNames(c("t1", "t2", "t4", "t5"), c("1", "2", "3", "4")))
    expect_equal(nodeLabels(phyd.sub4),
                 setNames(c("n6", "n7", "n9"), c("5", "6", "7")))
    expect_equal(edgeLength(phyd.sub4),
                 setNames(c(0.6, 0.4, 0.5, 0.7, 0.1, 0.2, 1.7),
                          c("0-5", "7-3", "7-4", "5-6", "6-1", "6-2", "5-7")))
    expect_equal(subset(phyd.alt, tips.exclude=3), phyd.sub4)
    expect_equal(subset(phyd.alt, tips.include=c("t1", "t2", "t4", "t5")),
                 phyd.sub4)
    expect_equal(subset(phyd.alt, tips.exclude="t3"), phyd.sub4)      
    ## check variants that should all return the original object
    expect_equal(phyd.alt, subset(phyd.alt))
    expect_equal(phyd.alt, subset(phyd.alt, tipLabels(phyd.alt)))
    expect_equal(phyd.alt, subset(phyd.alt, seq_len(nTips(phyd.alt))))
    expect_equal(phyd.alt, phyd.alt[tipLabels(phyd.alt)])
    expect_equal(phyd.alt, phyd.alt[seq_len(nTips(phyd.alt))])
    expect_equal(phyd.alt, phyd.alt[TRUE])
    ## error if only one valid tip requested
    expect_error(subset(phyd.alt, tips.include="t1"))
    expect_error(suppressWarnings(subset(phyd.alt, tips.include=c("t1", "t999"))))
    ## error if zero valid tips requested
    expect_error(suppressWarnings(subset(phyd.alt, tips.include="t999")))
    # subset tips that include an NA value
    ##TODO uncomment this after tdata is working right with scrambled order
    ##    tdata(phyd.alt)["t5", "a"] <- NA
    ##    tdata(phyd.sub2)["t5", "a"] <- NA
    ##    expect_equal(phyd.sub2, subset(phyd.alt, tips.include=c(2, 5)))
})

test_that("subset on extractTree", {
    # extract phylo4 from itself
    expect_equal(phy.alt, extractTree(phy.alt))
    # extract phylo4 from phylo4d
    expect_equal(phy.alt, extractTree(phyd.alt))
})
