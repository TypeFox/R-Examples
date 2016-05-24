#
# --- Test formatData.R ---
#

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

# create altered version such that each slot is out of order with
# respect to all others; methods should be able to handle this
phy.alt <- phy
phy.alt@label <- rev(phy@label)
phy.alt@edge <- phy@edge[c(6:9, 1:5), ]
phy.alt@edge.length <- phy@edge.length[c(7:9, 1:6)]
phy.alt@edge.label <- phy@edge.label[c(8:9, 1:7)]

# create data to add to phylo4 to create phylo4d, but with data rows out
# of order
set.seed(1)
nid.tip.r <- sample(nid.tip)
nid.int.r <- sample(nid.int)
nid.all.r <- sample(c(nid.tip, nid.int))
allDt <- data.frame(a=letters[nid.all.r], b=10*nid.all.r)
tipDt <- data.frame(c=letters[nid.tip.r], d=10*nid.tip.r)
nodDt <- data.frame(c=letters[nid.int.r], e=10*nid.int.r)
## set row.names as numeric node IDs (may be changed in tests below)
row.names(allDt) <- nid.all.r
row.names(tipDt) <- nid.tip.r
row.names(nodDt) <- nid.int.r

#-----------------------------------------------------------------------

context("test formatData")

## function(phy, dt, type=c("tip", "internal", "all"),
##   match.data=TRUE, rownamesAsLabels=FALSE,
##   label.type=c("rownames", "column"), label.column=1,
##   missing.data=c("fail", "warn", "OK"),
##   extra.data=c("warn", "OK", "fail"), keep.all=TRUE

test_that("works with data.frame", {
    ## vector data coerced to data.frame (colname dt)
    expect_equal(phylobase:::formatData(phy.alt, 1:5),
                 phylobase:::formatData(phy.alt, data.frame(dt=1:5)))
})

test_that("works with lists of vector", {
    ## list of vector data coerced to data.frame (colnames as given)
    expect_equal(phylobase:::formatData(phy.alt, list(a=1:5, b=6:10)),
                 phylobase:::formatData(phy.alt, data.frame(a=1:5, b=6:10)))
})

test_that("works factors", {
    ## factor data coerced to data.frame (colname dt)
    expect_equal(phylobase:::formatData(phy.alt, factor(letters[1:5])),
                 phylobase:::formatData(phy.alt, data.frame(dt=letters[1:5])))
})

test_that("works with data.frame and 2 columns", {
    ## matrix data coerced to data.frame (colnames V1, V2)
    expect_equal(phylobase:::formatData(phy.alt, matrix(1:10, ncol=2)),
                 phylobase:::formatData(phy.alt, data.frame(V1=1:5, V2=6:10)))
})

test_that("works with data.frame colname as given",  {
    ## matrix data coerced to data.frame (colname as given)
    expect_equal(phylobase:::formatData(phy.alt, matrix(1:10, ncol=2,
        dimnames=list(NULL, c("a", "b")))),
                 phylobase:::formatData(phy.alt, data.frame(a=1:5, b=6:10)))
})

test_that("fails with non-supported objects (i.e. a phylo4)", {
    ## error if dt is, say, a phylo4 object
    expect_error(phylobase:::formatData(phy.alt, phy.alt))
})

test_that("fails with column number is out of range", {
    ## error if column number is out of range
    expect_error(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip", match.data=FALSE,
        label.type="column", label.column=3))
})

test_that("fails with column name is wrong", {
    ## error if column name is wrong
    expect_error(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip", match.data=FALSE,
        label.type="column", label.column="foo"))
})


##
## matching options
##

test_that("matching options work as expected", {
    ## don't match (purely positional)
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=rev(nid.tip)), type="tip", match.data=FALSE),
        data.frame(a=c(1:5, rep(NA, 4)), row.names=nid.all))
    ## match on rownames (node numbers)
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=rev(nid.tip)), type="tip", match.data=TRUE),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=rev(nid.tip)), type="tip"), data.frame(a=c(5:1,
        rep(NA, 4)), row.names=nid.all))
    ## match on rownames (labels)
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=rev(lab.tip)), type="tip", match.data=TRUE),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## match on rownames (mixed node numbers and labels)
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])),
        type="tip", match.data=TRUE),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## but fails if rownamesAsLabels is TRUE
    expect_error(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])),
        type="tip", match.data=TRUE, rownamesAsLabels=TRUE))
})

##
##   label.type="column" and label.column=2
##

test_that("label.type=column works", {
    ## should ignore label (purely positional) and retain a label col
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip", match.data=FALSE,
        label.type="column", label.column=2),
        data.frame(a=c(1:5, rep(NA, 4)), lab=c(rev(nid.tip), rep(NA,
        4)), row.names=nid.all))
    ## match on label column (node numbers)
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip", match.data=TRUE,
        label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip",
        label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## match on label column (labels)
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(lab.tip)), type="tip", match.data=TRUE,
        label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=rev(lab.tip)), type="tip", match.data=TRUE,
        label.type="column", label.column="lab"),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## match on label column (mixed node numbers and labels)
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])), type="tip",
        match.data=TRUE, label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## but fails if rownamesAsLabels is TRUE
    expect_error(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        lab=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])),
        type="tip", match.data=TRUE, rownamesAsLabels=TRUE,
        label.type="column", label.column=2))
    ## try to match internal nodes when type='tips'
    expect_error(phylobase:::formatData(phy.alt, data.frame(a=1:5, row.names=4:8),
        type="tip"))
    ## and vice versa
    expect_error(phylobase:::formatData(phy.alt, data.frame(a=6:9, row.names=1:4),
        type="internal"))
})

##
## missing.data
##

test_that("behaves as expected with missing data", {
    ## force error conditions
    expect_error(phylobase:::formatData(phy.alt, data.frame(a=1:3), type="tip"))
    expect_error(phylobase:::formatData(phy.alt, data.frame(a=1:3), type="tip",
                                        missing.data="fail"))
    expect_warning(phylobase:::formatData(phy.alt, data.frame(a=1:3), type="tip",
                                        missing.data="warn"))
    ## missing data with matching
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=rev(nid.tip)[-1],
        row.names=rev(nid.tip)[-1]), type="tip", missing.data="OK"),
        data.frame(a=c(nid.tip[-5], rep(NA, 5))))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=rev(nid.int)[-1],
        row.names=rev(nid.int)[-1]), type="internal", missing.data="OK"),
        data.frame(a=c(rep(NA, 5), nid.int[-4], NA)))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=rev(nid.all)[-1],
        row.names=rev(nid.all)[-1]), type="all", missing.data="OK"),
        data.frame(a=c(nid.all[-9], NA)))
    ## missing data without matching
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=rev(nid.tip)[-1]),
        type="tip", match.data=FALSE, missing.data="OK"),
        data.frame(a=c(rev(nid.tip)[-1], rep(NA, 5))))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=rev(nid.int)[-1]),
        type="internal", match.data=FALSE, missing.data="OK"),
        data.frame(a=c(rep(NA, 5), rev(nid.int)[-1], NA)))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=rev(nid.all)[-1]),
        type="all", match.data=FALSE, missing.data="OK"),
        data.frame(a=c(rev(nid.all)[-1], NA)))
})

##
## extra.data
##

test_that("works as expected with extra data", {
    ## force error conditions
    expect_error(phylobase:::formatData(phy.alt, data.frame(a=1:3), type="tip",
        missing.data="fail"))
    expect_warning(phylobase:::formatData(phy.alt, data.frame(a=0:5, row.names=0:5),
        type="tip", missing="warn"), "not found in the tree")
    expect_warning(phylobase:::formatData(phy.alt, data.frame(a=0:5, row.names=0:5),
        type="tip"), "not found in the tree")
    ## extra data with matching
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=c(0L, rev(nid.tip)),
        row.names=c(0, rev(nid.tip))), type="tip", extra.data="OK"),
        data.frame(a=c(nid.tip, rep(NA, 4))))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=c(0L, rev(nid.int)),
        row.names=c(0, rev(nid.int))), type="internal", extra.data="OK"),
        data.frame(a=c(rep(NA, 5), nid.int)))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=c(0L, rev(nid.all)),
        row.names=c(0, rev(nid.all))), type="all", extra.data="OK"),
        data.frame(a=nid.all))
    ## extra data without matching
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:15),
        type="tip", match.data=FALSE, extra.data="OK"),
        data.frame(a=c(1:5, rep(NA, 4))))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:15),
        type="internal", match.data=FALSE, extra.data="OK"),
        data.frame(a=c(rep(NA, 5), 1:4)))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:15),
        type="all", match.data=FALSE, extra.data="OK"),
        data.frame(a=c(1:9)))
})

test_that("works as expected with both missing & extra data", {
    ## allow both extra.data and missing.data
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=0:3, row.names=0:3),
        type="tip", extra.data="OK", missing.data="OK"),
        data.frame(a=c(1:3, rep(NA, 6))))
})

##
## keep.all
##

test_that("keep.all works", {
    ## keep all rows
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=nid.tip), type="tip", keep.all=TRUE),
        data.frame(a=c(1:5, rep(NA, 4)), row.names=nid.all))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=nid.tip), type="tip"),
        data.frame(a=c(1:5, rep(NA, 4)), row.names=nid.all))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=6:9,
        row.names=nid.int), type="internal", keep.all=TRUE),
        data.frame(a=c(rep(NA, 5), 6:9), row.names=nid.all))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=6:9,
        row.names=nid.int), type="internal"),
        data.frame(a=c(rep(NA, 5), 6:9), row.names=nid.all))
    ## only keep 'type' rows
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=1:5,
        row.names=nid.tip), type="tip", keep.all=FALSE),
        data.frame(a=c(1:5), row.names=nid.tip))
    expect_equal(phylobase:::formatData(phy.alt, data.frame(a=6:9,
        row.names=nid.int), type="internal", keep.all=FALSE),
        data.frame(a=c(6:9), row.names=nid.int))
})

context("formatData with duplicated labels in object")

test_that("formatData works with duplicated labels", {
    ## Saving default options
    op <- phylobase.options()

    ## Changing default options
    phylobase.options(allow.duplicated.labels="ok")

    ## Creating phylo4 object with duplicated labels
    phy.dup <- phy.alt
    tipLabels(phy.dup)[2] <- tipLabels(phy.dup)[1]

    ## vector data coerced to data.frame (colname dt)
    expect_equal(phylobase:::formatData(phy.dup, 1:5),
        phylobase:::formatData(phy.dup, data.frame(dt=1:5)))
    ## list of vector data coerced to data.frame (colnames as given)
    expect_equal(phylobase:::formatData(phy.dup, list(a=1:5, b=6:10)),
        phylobase:::formatData(phy.dup, data.frame(a=1:5, b=6:10)))
    ## factor data coerced to data.frame (colname dt)
    expect_equal(phylobase:::formatData(phy.dup, factor(letters[1:5])),
        phylobase:::formatData(phy.dup, data.frame(dt=letters[1:5])))
    ## matrix data coerced to data.frame (colnames V1, V2)
    expect_equal(phylobase:::formatData(phy.dup, matrix(1:10, ncol=2)),
        phylobase:::formatData(phy.dup, data.frame(V1=1:5, V2=6:10)))
    ## matrix data coerced to data.frame (colname as given)
    expect_equal(phylobase:::formatData(phy.dup, matrix(1:10, ncol=2,
        dimnames=list(NULL, c("a", "b")))),
        phylobase:::formatData(phy.dup, data.frame(a=1:5, b=6:10)))
    ## error if dt is, say, a phylo4 object
    expect_error(phylobase:::formatData(phy.dup, phy.dup))

    #
    # matching options
    #

    ## don't match (purely positional)
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=rev(nid.tip)), type="tip", match.data=FALSE),
        data.frame(a=c(1:5, rep(NA, 4)), row.names=nid.all))
    ## match on rownames (node numbers)
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=rev(nid.tip)), type="tip", match.data=TRUE),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=rev(nid.tip)), type="tip"), data.frame(a=c(5:1,
        rep(NA, 4)), row.names=nid.all))
    ## match on rownames (labels)
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=c(1,3,4,5),
        row.names=rev(lab.tip[-2])), type="tip", match.data=TRUE),
        data.frame(a=c(5,5,4,3,1, rep(NA, 4)), row.names=nid.all))
    ## match on rownames (mixed node numbers and labels)
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=c(1,2,3,4,5),
        row.names=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])),
        type="tip", match.data=TRUE),
        data.frame(a=c(5,4,3,2,1, rep(NA, 4)), row.names=nid.all))
    ## but fails if rownamesAsLabels is TRUE
    expect_error(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])),
        type="tip", match.data=TRUE, rownamesAsLabels=TRUE))

    ##
    ##   label.type="column" and label.column=2
    ##

    ## should ignore label (purely positional) and retain a label col
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip", match.data=FALSE,
        label.type="column", label.column=2),
        data.frame(a=c(1:5, rep(NA, 4)), lab=c(rev(nid.tip), rep(NA,
        4)), row.names=nid.all))
    ## match on label column (node numbers)
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip", match.data=TRUE,
        label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        lab=rev(nid.tip)), type="tip",
        label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## match on label column (labels)
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:4,
        lab=rev(lab.tip[-2])), type="tip", match.data=TRUE,
        label.type="column", label.column=2),
        data.frame(a=as.integer(c(4, 4:1, rep(NA, 4))), row.names=nid.all))
    ## match on label column (mixed node numbers and labels)
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        lab=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])), type="tip",
        match.data=TRUE, label.type="column", label.column=2),
        data.frame(a=c(5:1, rep(NA, 4)), row.names=nid.all))
    ## but fails if rownamesAsLabels is TRUE
    expect_error(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        lab=c(rev(lab.tip)[1:3], rev(nid.tip)[4:5])),
        type="tip", match.data=TRUE, rownamesAsLabels=TRUE,
        label.type="column", label.column=2))

    ## try to match internal nodes when type='tips'
    expect_error(phylobase:::formatData(phy.dup, data.frame(a=1:5, row.names=4:8),
        type="tip"))
    ## and vice versa
    expect_error(phylobase:::formatData(phy.dup, data.frame(a=6:9, row.names=1:4),
        type="internal"))

    ##
    ## missing.data
    ##

    ## force error conditions
    expect_error(phylobase:::formatData(phy.dup, data.frame(a=1:3), type="tip"))
    expect_error(phylobase:::formatData(phy.dup, data.frame(a=1:3), type="tip",
        missing.data="fail"))
    expect_warning(phylobase:::formatData(phy.dup, data.frame(a=1:3), type="tip",
        missing.data="warn"))
    ## missing data with matching
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=rev(nid.tip)[-1],
        row.names=rev(nid.tip)[-1]), type="tip", missing.data="OK"),
        data.frame(a=c(nid.tip[-5], rep(NA, 5))))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=rev(nid.int)[-1],
        row.names=rev(nid.int)[-1]), type="internal", missing.data="OK"),
        data.frame(a=c(rep(NA, 5), nid.int[-4], NA)))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=rev(nid.all)[-1],
        row.names=rev(nid.all)[-1]), type="all", missing.data="OK"),
        data.frame(a=c(nid.all[-9], NA)))
    ## missing data without matching
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=rev(nid.tip)[-1]),
        type="tip", match.data=FALSE, missing.data="OK"),
        data.frame(a=c(rev(nid.tip)[-1], rep(NA, 5))))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=rev(nid.int)[-1]),
        type="internal", match.data=FALSE, missing.data="OK"),
        data.frame(a=c(rep(NA, 5), rev(nid.int)[-1], NA)))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=rev(nid.all)[-1]),
        type="all", match.data=FALSE, missing.data="OK"),
        data.frame(a=c(rev(nid.all)[-1], NA)))

    ##
    ## extra.data
    ##

    ## force error conditions
    expect_error(phylobase:::formatData(phy.dup, data.frame(a=1:3), type="tip",
        missing.data="fail"))
    expect_warning(phylobase:::formatData(phy.dup, data.frame(a=0:5, row.names=0:5),
        type="tip", missing="warn"))
    expect_warning(phylobase:::formatData(phy.dup, data.frame(a=0:5, row.names=0:5),
        type="tip"))
    ## extra data with matching
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=c(0L, rev(nid.tip)),
        row.names=c(0, rev(nid.tip))), type="tip", extra.data="OK"),
        data.frame(a=c(nid.tip, rep(NA, 4))))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=c(0L, rev(nid.int)),
        row.names=c(0, rev(nid.int))), type="internal", extra.data="OK"),
        data.frame(a=c(rep(NA, 5), nid.int)))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=c(0L, rev(nid.all)),
        row.names=c(0, rev(nid.all))), type="all", extra.data="OK"),
        data.frame(a=nid.all))
    ## extra data without matching
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:15),
        type="tip", match.data=FALSE, extra.data="OK"),
        data.frame(a=c(1:5, rep(NA, 4))))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:15),
        type="internal", match.data=FALSE, extra.data="OK"),
        data.frame(a=c(rep(NA, 5), 1:4)))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:15),
        type="all", match.data=FALSE, extra.data="OK"),
        data.frame(a=c(1:9)))

    ## allow both extra.data and missing.data
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=0:3, row.names=0:3),
        type="tip", extra.data="OK", missing.data="OK"),
        data.frame(a=c(1:3, rep(NA, 6))))

    ##
    ## keep.all
    ##

    ## keep all rows
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=nid.tip), type="tip", keep.all=TRUE),
        data.frame(a=c(1:5, rep(NA, 4)), row.names=nid.all))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=nid.tip), type="tip"),
        data.frame(a=c(1:5, rep(NA, 4)), row.names=nid.all))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=6:9,
        row.names=nid.int), type="internal", keep.all=TRUE),
        data.frame(a=c(rep(NA, 5), 6:9), row.names=nid.all))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=6:9,
        row.names=nid.int), type="internal"),
        data.frame(a=c(rep(NA, 5), 6:9), row.names=nid.all))
    ## only keep 'type' rows
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=1:5,
        row.names=nid.tip), type="tip", keep.all=FALSE),
        data.frame(a=c(1:5), row.names=nid.tip))
    expect_equal(phylobase:::formatData(phy.dup, data.frame(a=6:9,
        row.names=nid.int), type="internal", keep.all=FALSE),
        data.frame(a=c(6:9), row.names=nid.int))

    ## restoring default options
    phylobase.options(op)
})
