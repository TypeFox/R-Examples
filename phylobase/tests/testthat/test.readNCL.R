#
# --- Test readNCL.R ---
#

### Get all the test files
if (Sys.getenv("RCMDCHECK") == FALSE) {
    pth <- file.path(getwd(), "..", "inst", "nexusfiles")
} else {
    pth <- system.file(package="phylobase", "nexusfiles")
}

## co1.nex -- typical output from MrBayes. Contains 2 identical trees, the first
## one having posterior probabilities as node labels
co1File <- file.path(pth, "co1.nex")

## MultiLineTrees.nex -- 2 identical trees stored on several lines
multiLinesFile <- file.path(pth, "MultiLineTrees.nex")

## treeWithDiscreteData.nex -- Mesquite file with discrete data
treeDiscDt <- file.path(pth, "treeWithDiscreteData.nex")

## treeWithPolyExcludedData.nex -- Mesquite file with polymorphic and excluded
##  characters
treePolyDt <- file.path(pth, "treeWithPolyExcludedData.nex")

## treeWithContinuousData.nex -- Mesquite file with continuous characters
treeContDt <- file.path(pth, "treeWithContinuousData.nex")

## treeWithDiscAndContData.nex -- Mesquite file with both discrete and
##    continuous data
treeDiscCont <- file.path(pth, "treeWithDiscAndContData.nex")

## noStateLabels.nex -- Discrete characters with missing state labels
noStateLabels <- file.path(pth, "noStateLabels.nex")

## Newick trees
newick <- file.path(pth, "newick.tre")

## Test with trees that don't include all the taxa listed in TAXA block
treeSubset <- file.path(pth, "testSubsetTaxa.nex")

## Contains representation of data associated with continuous data
ExContDataFile <- file.path(pth, "ExContData.Rdata")


stopifnot(file.exists(co1File))
stopifnot(file.exists(treeDiscDt))
stopifnot(file.exists(multiLinesFile))
stopifnot(file.exists(treePolyDt))
stopifnot(file.exists(treeContDt))
stopifnot(file.exists(treeDiscCont))
stopifnot(file.exists(ExContDataFile))
stopifnot(file.exists(noStateLabels))
stopifnot(file.exists(treeSubset))

op <- phylobase.options()


## function (file, simplify=TRUE, type=c("all", "tree", "data"),
##   char.all=FALSE, polymorphic.convert=TRUE, levels.uniform=TRUE,
##   check.node.labels=c("keep", "drop", "asdata"))



## ########### CO1 -- MrBayes file -- tree only

## Tree properties
## Labels
labCo1 <- c("Cow", "Seal", "Carp", "Loach", "Frog", "Chicken", "Human",
            "Mouse", "Rat", "Whale", NA, NA, NA, NA, NA, NA, NA, NA)
names(labCo1) <- 1:18
## Edge lengths
eLco1 <- c(0.143336, 0.225087, 0.047441, 0.055934, 0.124549, 0.204809, 0.073060, 0.194575,
           0.171296, 0.222039, 0.237101, 0.546258, 0.533183, 0.154442, 0.134574, 0.113163,
           0.145592)
names(eLco1) <- c("11-1", "11-2", "11-12", "12-13", "13-14", "14-15", "15-16", "16-17", "17-3",
                  "17-4", "16-5", "15-6", "14-7", "13-18", "18-8", "18-9", "12-10")
## Node types
nTco1 <-  c("tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip",
            "tip", "internal", "internal", "internal", "internal", "internal",
            "internal", "internal", "internal")
names(nTco1) <- 1:18
## Label values
lVco1 <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.93, 0.88, 0.99, 1.00,
           0.76, 1.00, 1.00)
context("readNCL can deal with simple NEXUS files (tree only)")
test_that("file with 2 trees (warning normal)", {
    ## Read trees
    co1 <- suppressWarnings(readNCL(file=co1File, check.node.labels="asdata"))
    ## Tree 1
    co1Tree1 <- co1[[1]]
    edgeNm <- paste(edges(co1Tree1)[, "ancestor"], edges(co1Tree1)[, "descendant"], sep = "-")
    expect_equal(labels(co1Tree1), labCo1)     # check labels
    expect_equal(edgeLength(co1Tree1), eLco1[edgeNm])  # check edge lengths
    expect_equal(nodeType(co1Tree1), nTco1)    # check node types
    expect_equal(as(co1Tree1, "data.frame")$labelValues, lVco1) # check label value
    ## Tree 2
    co1Tree2 <- co1[[2]]
    expect_equal(labels(co1Tree2), labCo1)     # check labels
    expect_equal(edgeLength(co1Tree2), eLco1[edgeNm])  # check edge lengths
    expect_equal(nodeType(co1Tree2), nTco1)    # check node types
})

test_that("test option simplify", {
    ## Check option simplify
    co1 <- readNCL(file=co1File, check.node.labels="asdata", simplify=TRUE)
    edgeNm <- paste(edges(co1)[, "ancestor"], edges(co1)[, "descendant"], sep = "-")
    expect_equal(length(co1), as.integer(1))   # make sure there is only one tree
    expect_equal(labels(co1), labCo1)          # check labels
    expect_equal(edgeLength(co1), eLco1[edgeNm])       # check edge lengths
    expect_equal(nodeType(co1), nTco1)         # check node type
    expect_equal(as(co1, "data.frame")$labelValues, lVco1)  # check label values
})

test_that("test option check.node.labels", {
    ## Check option check.node.labels
    phylobase.options(allow.duplicated.labels="fail")
    expect_error(readNCL(file=co1File, check.node.labels="keep")) # fail because labels aren't unique
    phylobase.options(op)
    phylobase.options(allow.duplicated.labels="ok")
    co1 <- readNCL(file=co1File, check.node.labels="keep", simplify=TRUE)
    expect_equal(nodeLabels(co1),
                 setNames(c(NA, "0.93", "0.88", "0.99", "1.00", "0.76", "1.00", "1.00"),
                          11:18))
    phylobase.options(op)
    co1 <- readNCL(file=co1File, check.node.labels="drop", simplify=TRUE)
    edgeNm <- paste(edges(co1)[, "ancestor"], edges(co1)[, "descendant"], sep = "-")
    expect_equal(labels(co1), labCo1)          # check labels
    expect_equal(edgeLength(co1), eLco1[edgeNm])       # check edge lengths
    expect_equal(nodeType(co1), nTco1)         # check node type
    expect_equal(as(co1, "data.frame")$labelValues, NULL)  # check label values don't exist
})

test_that("labelled root", {
    tmp_file <- tempfile()
    cat("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;", file = tmp_file)
    ape_tree <- as(ape::read.tree(file = tmp_file), "phylo4")
    ph4_tree <- readNewick(file = tmp_file)
    expect_equal(tipLabels(ape_tree), tipLabels(ph4_tree))
    expect_equal(nodeLabels(ape_tree), nodeLabels(ph4_tree))
    expect_equal(sort(edgeLength(ape_tree)), sort(edgeLength(ape_tree)))
})

test_that("readNCL can handle multi line files", {
    ## ########### Mutli Lines -- tree only
    multiLines <- readNCL(file=multiLinesFile)
    ## load correct representation and make sure that the trees read
    ## match it
    ml <- rncl::read_nexus_phylo(file = multiLinesFile)
    ml1 <- as(ml[[1]], "phylo4")
    ml2 <- as(ml[[2]], "phylo4")
    expect_equal(tipLabels(multiLines[[1]]), tipLabels(ml1))
    expect_equal(tipLabels(multiLines[[2]]), tipLabels(ml2))
    expect_equivalent(sort(edgeLength(multiLines[[1]])), sort(edgeLength(ml1)))
    expect_equivalent(sort(edgeLength(multiLines[[2]])), sort(edgeLength(ml2)))
    expect_equal(nodeType(multiLines[[1]]), nodeType(ml1))
    expect_equal(nodeType(multiLines[[2]]), nodeType(ml2))
})

## ########### Tree + data -- file from Mesquite
context("readNCL can handle files with tree & data")
## tree properties
labTr <-  c("Myrmecocystussemirufus", "Myrmecocystusplacodops",
            "Myrmecocystusmendax", "Myrmecocystuskathjuli",
            "Myrmecocystuswheeleri", "Myrmecocystusmimicus",
            "Myrmecocystusdepilis", "Myrmecocystusromainei",
            "Myrmecocystusnequazcatl", "Myrmecocystusyuma",
            "Myrmecocystuskennedyi", "Myrmecocystuscreightoni",
            "Myrmecocystussnellingi", "Myrmecocystustenuinodis",
            "Myrmecocystustestaceus", "Myrmecocystusmexicanus",
            "Myrmecocystuscfnavajo", "Myrmecocystusnavajo",
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
names(labTr) <- 1:35
eTr <- c(NA, 1.699299, 12.300701, 0.894820, 0.836689, 10.569191, 4.524387, 6.044804,
         0.506099, 0.198842, 0.689044, 4.650818, 2.926053, 1.724765, 1.724765, 4.255993,
         1.083870, 1.083870, 0.802512, 2.027251, 2.708942, 2.708942, 0.284767, 4.451425,
         2.257581, 2.193845, 2.193845, 8.635503, 2.770378, 2.770378, 8.275077, 5.724923,
         2.855375, 2.869547, 2.869547)
names(eTr) <- c("0-19", "19-20", "20-15", "20-21", "21-22", "22-12", "22-23", "23-11", "23-24",
                "24-25", "25-26", "26-3", "26-27", "27-1", "27-2", "25-28", "28-4", "28-5",
                "24-29", "29-30", "30-6", "30-7", "29-31", "31-10", "31-32", "32-8", "32-9",
                "21-33", "33-13", "33-14", "19-34", "34-16", "34-35", "35-17", "35-18")
nTtr <- c("tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip",
          "tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip", "tip",
          "root", "internal", "internal", "internal", "internal", "internal",
          "internal", "internal", "internal", "internal", "internal",
          "internal", "internal", "internal", "internal", "internal",
          "internal")
names(nTtr) <- 1:35
## data to test against
dtTest1 <- data.frame(time = factor(c(2,1,0,0,0,0,2,0,2,0,0,0,0,1,1,1,0,1)),
                      subgenus = factor(c(2,1,0,0,0,0,2,0,2,0,0,0,0,1,1,2,0,1)))
row.names(dtTest1) <- c("Myrmecocystuscfnavajo","Myrmecocystuscreightoni",
                        "Myrmecocystusdepilis","Myrmecocystuskathjuli",
                        "Myrmecocystuskennedyi","Myrmecocystusmendax",
                        "Myrmecocystusmexicanus","Myrmecocystusmimicus",
                        "Myrmecocystusnavajo","Myrmecocystusnequazcatl",
                        "Myrmecocystusplacodops","Myrmecocystusromainei",
                        "Myrmecocystussemirufus","Myrmecocystussnellingi",
                        "Myrmecocystustenuinodis","Myrmecocystustestaceus",
                        "Myrmecocystuswheeleri","Myrmecocystusyuma")
dtTest2 <- dtTest1
levels(dtTest2$time) <- c("diurnal", "crepuscular", "nocturnal")
levels(dtTest2$subgenus) <- c("Endiodioctes", "Eremnocystus", "Myrmecocystus")
p4 <- "phylo4"
p4d <- "phylo4d"
attributes(p4) <- attributes(p4d) <- list(package="phylobase")

test_that("readNCL can deal with the tree only", {
    ## Tree only
    tr <- readNCL(file=treeDiscDt, type="tree")
    tr2 <- rncl::read_nexus_phylo(file = treeDiscDt)
    tr2 <- as(tr2, "phylo4")
    expect_equal(labels(tr), labTr)   # check labels
    expect_equal(nodeType(tr), nTtr)  # check node types
    expect_equal(class(tr), p4)       # check class
    expect_equal(edgeLength(tr), edgeLength(tr2)[names(edgeLength(tr))])
})

test_that("readNCL can deal with data only", {
    ## Data only
    dt1 <- readNCL(file=treeDiscDt, type="data", return.labels=FALSE,
                   levels.uniform=FALSE)
    expect_equal(dt1, dtTest1)
    dt2 <- readNCL(file=treeDiscDt, type="data", return.labels=TRUE,
                   levels.uniform=FALSE)
    expect_equal(dt2, dtTest2)
})

test_that("readNCL can deal with tree + data", {
    ## Tree + Data
    trDt1 <- readNCL(file=treeDiscDt, type="all", return.labels=FALSE,
                     levels.uniform=FALSE)
    expect_equal(labels(trDt1), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trDt1)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trDt1), nTtr)  # check node types
    expect_equal(class(trDt1), p4d)      # check class
    expect_equal(tdata(trDt1, type="tip")[rownames(dtTest1), ], dtTest1)
    trDt2 <- readNCL(file=treeDiscDt, type="all", return.labels=TRUE,
                       levels.uniform=FALSE)
    expect_equal(labels(trDt2), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trDt2)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trDt2), nTtr)  # check node types
    expect_equal(class(trDt2), p4d)      # check class
    expect_equal(tdata(trDt2, type="tip")[rownames(dtTest2), ], dtTest2)
})


## ########## Tree + Data -- Test for polymorphic.convert, levels.uniform and char.all
## data to test against
## dtTest 3 -- levels.uniform=FALSE, return.labels=FALSE, polymorphic.convert=FALSE
dtPoly1 <- data.frame(Test1=factor(c(0,0,1,1,0,NA,1,1,1,0,0,NA,1,1,NA,0,1, NA)),
                      Test2=factor(c(0,0,0,0,0,NA,0,1,0,1,1,"{0,1}",NA,0,NA,0,"{0,1}",1)),
                      Test3=factor(c(1,1,1,0,0,0,2,"{0,1,2}",0,NA,0,"{0,1}",0,1,0,0,"{0,1,2}",1)),
                      row.names=c("Myrmecocystussemirufus","Myrmecocystusplacodops",
                          "Myrmecocystusmendax","Myrmecocystuskathjuli",
                          "Myrmecocystuswheeleri","Myrmecocystusmimicus",
                          "Myrmecocystusdepilis","Myrmecocystusromainei",
                          "Myrmecocystusnequazcatl","Myrmecocystusyuma",
                          "Myrmecocystuskennedyi","Myrmecocystuscreightoni",
                          "Myrmecocystussnellingi","Myrmecocystustenuinodis",
                          "Myrmecocystustestaceus","Myrmecocystusmexicanus",
                          "Myrmecocystuscfnavajo","Myrmecocystusnavajo"))
## dtPoly2 -- levels.uniform=FALSE, return.labels=FALSE, polymorphic.convert=TRUE
dtPoly2 <- dtPoly1
dtPoly2[c(12,17),2] <- NA
dtPoly2[c(8,12,17),3] <- NA
dtPoly2$Test1 <- factor(dtPoly2$Test1)
dtPoly2$Test2 <- factor(dtPoly2$Test2)
dtPoly2$Test3 <- factor(dtPoly2$Test3)
## dtPoly3 -- levels.uniform=FALSE, return.labels=TRUE, polymorphic.convert=TRUE
dtPoly3 <- dtPoly2
levels(dtPoly3$Test1) <- c("test1A", "test1B")
levels(dtPoly3$Test2) <- c("test2A", "test2B")
levels(dtPoly3$Test3) <- c("test3A", "test3B", "test3C")
## dtPoly4 -- levels.uniform=FALSE, return.labels=TRUE, polymorphic.convert=FALSE
##    not yet implemented

## dtPoly5 -- levels.uniform=TRUE, return.labels=FALSE, polymorphic.convert=FALSE
dtPoly5 <- dtPoly1
levels(dtPoly5$Test1) <- levels(dtPoly5$Test2) <- levels(dtPoly5$Test3) <-
    union(levels(dtPoly1$Test1), c(levels(dtPoly1$Test2), levels(dtPoly1$Test3)))
## dtPoly6 -- levels.uniform=TRUE, return.labels=FALSE, polymorphic.convert=TRUE
dtPoly6 <- dtPoly2
levels(dtPoly6$Test1) <- levels(dtPoly6$Test2) <- levels(dtPoly6$Test3) <-
    union(levels(dtPoly2$Test1), c(levels(dtPoly2$Test2), levels(dtPoly2$Test3)))
## dtPoly7 -- levels.uniform=TRUE, return.labels=TRUE, polymorphic.convert=FALSE
##    not yet implemented

## dtPoly8 -- levels.uniform=TRUE, return.labels=TRUE, polymorphic.convert=TRUE
dtPoly8 <- dtPoly3
levels(dtPoly8$Test1) <- levels(dtPoly8$Test2) <- levels(dtPoly8$Test3) <-
    union(levels(dtPoly3$Test1), c(levels(dtPoly3$Test2), levels(dtPoly3$Test3)))
## dtPoly5F -- char.all=FALSE, levels.uniform=TRUE, return.labels=FALSE, polymorphic.convert=FALSE
dtPoly5F <- dtPoly1[, 1:2]
levels(dtPoly5F$Test1) <- levels(dtPoly5F$Test2) <-
    union(levels(dtPoly1$Test1), levels(dtPoly1$Test2))
## dtPoly6F -- char.all=FALSE, levels.uniform=TRUE, return.labels=FALSE, polymorphic.convert=TRUE
dtPoly6F <- dtPoly2[, 1:2]
levels(dtPoly6F$Test1) <- levels(dtPoly6F$Test2) <-
    union(levels(dtPoly2$Test1), levels(dtPoly2$Test2))
## dtPoly8F -- char.all=FALSE, levels.uniform=TRUE, return.labels=TRUE, polymorphic.convert=TRUE
dtPoly8F <- dtPoly3[, 1:2]
levels(dtPoly8F$Test1) <- levels(dtPoly8F$Test2) <-
    union(levels(dtPoly3$Test1), levels(dtPoly3$Test2))

test_that("char.all=TRUE, levels.uniform=FALSE, return.labels=FALSE, polymorphic.convert=FALSE", {
    trChr1 <- readNCL(file=treePolyDt, type="all", polymorphic.convert=FALSE,
                      levels.uniform=FALSE, char.all=TRUE, return.labels=FALSE)
    expect_equal(labels(trChr1), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trChr1)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trChr1), nTtr)  # check node types
    expect_equal(class(trChr1), p4d)      # check class
    expect_equal(tdata(trChr1, "tip"), dtPoly1[tipLabels(trChr1), ])
})

test_that("char.all=TRUE, levels.uniform=FALSE, return.labels=FALSE, polymorphic.convert=TRUE", {
    trChr2 <- readNCL(file=treePolyDt, type="all", polymorphic.convert=TRUE,
                      levels.uniform=FALSE, return.labels=FALSE, char.all=TRUE)
    expect_equal(labels(trChr2), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trChr2)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trChr2), nTtr)  # check node types
    expect_equal(class(trChr2), p4d)      # check class
    expect_equal(tdata(trChr2, "tip"), dtPoly2[tipLabels(trChr2), ])
})

test_that("char.all=TRUE, levels.uniform=FALSE, return.labels=TRUE, polymorphic.convert=TRUE", {
    trChr3 <- readNCL(file=treePolyDt, type="all", polymorphic.convert=TRUE,
                      levels.uniform=FALSE, char.all=TRUE, return.labels=TRUE)
    expect_equal(labels(trChr3), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trChr3)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trChr3), nTtr)  # check node types
    expect_equal(class(trChr3), p4d)      # check class
    expect_equal(tdata(trChr3, "tip"), dtPoly3[tipLabels(trChr3), ])
})

test_that("char.all=TRUE, levels.uniform=FALSE, return.labels=TRUE, polymorphic.convert=FALSE is not yet implemented", {
## trChr4 <-
    expect_error(readNCL(file=treePolyDt, type="all",
                         levels.uniform=FALSE,
                         return.labels=TRUE,
                         polymorphic.convert=FALSE))
})

test_that("char.all=TRUE, levels.uniform=TRUE, return.labels=FALSE, polymorphic.convert=FALSE", {
    trChr5 <- readNCL(file=treePolyDt, type="all", polymorphic.convert=FALSE,
                        levels.uniform=TRUE, char.all=TRUE, return.labels=FALSE)
    expect_equal(labels(trChr5), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trChr5)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trChr5), nTtr)  # check node types
    expect_equal(class(trChr5), p4d)      # check class
    expect_equal(tdata(trChr5, "tip"), dtPoly5[tipLabels(trChr5), ])
})

test_that("char.all=TRUE, levels.uniform=TRUE, return.labels=FALSE, polymorphic.convert=TRUE", {
    trChr6 <- readNCL(file=treePolyDt, type="all", polymorphic.convert=TRUE,
                      levels.uniform=TRUE, char.all=TRUE, return.labels=FALSE)
    expect_equal(labels(trChr6), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trChr6)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trChr6), nTtr)  # check node types
    expect_equal(class(trChr6), p4d)      # check class
    expect_equal(tdata(trChr6, "tip"), dtPoly6[tipLabels(trChr6), ])
})

test_that("char.all=TRUE, levels.uniform=TRUE, return.labels=TRUE, polymorphic.convert=FALSE is not yet implemented", {
    ## trChr7 <-
    expect_error(readNCL(file=treePolyDt, type="all", char.all=TRUE,
                         levels.uniform=TRUE,
                         return.labels=TRUE,
                         polymorphic.convert=FALSE))
})

test_that("char.all=TRUE, levels.uniform=TRUE, return.labels=TRUE, polymorphic.convert=TRUE", {
    trChr8 <- readNCL(file=treePolyDt, type="all", char.all=TRUE,
                      levels.uniform=TRUE,
                      return.labels=TRUE,
                      polymorphic.convert=TRUE)
    expect_equal(labels(trChr8), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trChr8)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trChr8), nTtr)  # check node types
    expect_equal(class(trChr8), p4d)      # check class
    expect_equal(tdata(trChr8, "tip"), dtPoly8[tipLabels(trChr8), ])
})

## -- with char.all=FALSE
test_that("char.all=FALSE, levels.uniform=FALSE, return.labels=FALSE, polymorphic.convert=FALSE", {
    trChr1F <- readNCL(file=treePolyDt, type="all", polymorphic.convert=FALSE,
                       levels.uniform=FALSE, char.all=FALSE, return.labels=FALSE)
    expect_equal(labels(trChr1F), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trChr1F)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trChr1F), nTtr)  # check node types
    expect_equal(class(trChr1F), p4d)      # check class
    expect_equal(tdata(trChr1F, "tip"), dtPoly1[tipLabels(trChr1F), 1:2])
})

test_that("char.all=FALSE, levels.uniform=FALSE, return.labels=FALSE, polymorphic.convert=TRUE", {
    trChr2F <- readNCL(file=treePolyDt, type="all", polymorphic.convert=TRUE,
                       levels.uniform=FALSE, return.labels=FALSE, char.all=FALSE)
    expect_equal(labels(trChr2F), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trChr2F)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trChr2F), nTtr)  # check node types
    expect_equal(class(trChr2F), p4d)      # check class
    expect_equal(tdata(trChr2F, "tip"), dtPoly2[tipLabels(trChr2F), 1:2])
})

test_that("char.all=FALSE, levels.uniform=FALSE, return.labels=TRUE, polymorphic.convert=TRUE", {
    trChr3F <- readNCL(file=treePolyDt, type="all", polymorphic.convert=TRUE,
                        levels.uniform=FALSE, char.all=FALSE, return.labels=TRUE)
    expect_equal(labels(trChr3F), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trChr3F)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trChr3F), nTtr)  # check node types
    expect_equal(class(trChr3F), p4d)      # check class
    expect_equal(tdata(trChr3F, "tip"), dtPoly3[tipLabels(trChr3F), 1:2])
})

test_that("char.all=FALSE, levels.uniform=FALSE, return.labels=TRUE, polymorphic.convert=FALSE is not yet implemented", {
    ## trChr4F <-
    expect_error(readNCL(file=treePolyDt, type="all",
                         levels.uniform=FALSE,
                         return.labels=TRUE,
                         polymorphic.convert=FALSE))
})

test_that("char.all=FALSE, levels.uniform=TRUE, return.labels=FALSE, polymorphic.convert=FALSE", {
    trChr5F <- readNCL(file=treePolyDt, type="all", polymorphic.convert=FALSE,
                       levels.uniform=TRUE, char.all=FALSE, return.labels=FALSE)
    expect_equal(labels(trChr5F), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trChr5F)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trChr5F), nTtr)  # check node types
    expect_equal(class(trChr5F), p4d)      # check class
    expect_equal(tdata(trChr5F, "tip"), dtPoly5F[tipLabels(trChr5F), ])
})

test_that("char.all=FALSE, levels.uniform=TRUE, return.labels=FALSE, polymorphic.convert=TRUE", {
    trChr6F <- readNCL(file=treePolyDt, type="all", polymorphic.convert=TRUE,
                        levels.uniform=TRUE, char.all=FALSE, return.labels=FALSE)
    expect_equal(labels(trChr6F), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trChr6F)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trChr6F), nTtr)  # check node types
    expect_equal(class(trChr6F), p4d)      # check class
    expect_equal(tdata(trChr6F, "tip"), dtPoly6F[tipLabels(trChr6F), ])
})

test_that("char.all=FALSE, levels.uniform=TRUE, return.labels=TRUE, polymorphic.convert=FALSE is not yet implemented", {
    ## trChr7F <-
    expect_error(readNCL(file=treePolyDt, type="all", char.all=FALSE,
                         levels.uniform=TRUE,
                         return.labels=TRUE,
                         polymorphic.convert=FALSE))
})

test_that("char.all=FALSE, levels.uniform=TRUE, return.labels=TRUE, polymorphic.convert=TRUE", {
    trChr8F <- readNCL(file=treePolyDt, type="all", char.all=FALSE,
                       levels.uniform=TRUE,
                       return.labels=TRUE,
                       polymorphic.convert=TRUE)
    expect_equal(labels(trChr8F), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trChr8F)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trChr8F), nTtr)  # check node types
    expect_equal(class(trChr8F), p4d)      # check class
    expect_equal(tdata(trChr8F, "tip"), dtPoly8F[tipLabels(trChr8F), ])
})

## ########## Tree + Data -- test with continuous Characters
test_that("test of readNCL with tree data, with continuous characters", {
    DtCont <- readNCL(file=treeContDt, type="data")
    trDtCont <- readNCL(file=treeContDt, type="all")
    load(ExContDataFile)
    expect_equal(DtCont, ExContData[rownames(DtCont), ])
    expect_equal(tdata(trDtCont, "tip"), ExContData[tipLabels(trDtCont), ])
    expect_equal(labels(trDtCont), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trDtCont)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trDtCont), nTtr)  # check node types
    expect_equal(class(trDtCont), p4d)      # check class
})


## ########## Tree + Data -- both types (Discrete & Continuous)
test_that("tree + data for both types (discrete & continuous)", {
    dtDiscCont <- readNCL(file=treeDiscCont, type="data", levels.uniform=FALSE)
    trDtDiscCont <- readNCL(file=treeDiscCont, type="all", levels.uniform=FALSE)
    load(ExContDataFile)
    dtDiscContTest <- cbind(ExContData, dtTest2[rownames(ExContData), ])
    expect_equal(dtDiscCont, dtDiscContTest[rownames(dtDiscCont), ])
    expect_equal(tdata(trDtDiscCont, "tip"), dtDiscContTest[tipLabels(trDtDiscCont), ])
    expect_equal(labels(trDtDiscCont), labTr)   # check labels
    expect_equivalent(sort(edgeLength(trDtDiscCont)), sort(eTr)) # check edge lengths
    expect_equal(nodeType(trDtDiscCont), nTtr)  # check node types
    expect_equal(class(trDtDiscCont), p4d)      # check class
})

## ########### Check for proper handling of missing files
test_that("readNCL can handle missing files", {
    expect_error(readNCL(file="foo.bar"), regexp="doesn't exist")
})

## ########### Check behavior in case of missing state labels
test_that("readNCL warns in case of missing state labels", {
    expect_warning(readNCL(file=noStateLabels, return.labels=TRUE),
                   regexp="state labels are missing")
})

test_that("readNCL warns in case of missing state labels", {
    expect_warning(dtNoSt <- readNCL(file=noStateLabels, type="data",
                                     return.labels=TRUE),
                   regexp="state labels are missing")
    expect_equal(dtNoSt$char1, factor(c(1,2,0,1)))
})

## ########### Newick files
context("test with Newick files")
## Tree representation
labNew <- c("a", "b", "c", NA, NA)
names(labNew) <- 1:5
eLnew <- c(NA, 1, 4, 2, 3)
names(eLnew) <- c("0-4", "4-1", "4-5", "5-2", "5-3")
nTnew <- c("tip", "tip", "tip", "root", "internal")
names(nTnew) <- 1:5

test_that("check.node.labels='drop' with readNCL", {
    newTr <- readNCL(file=newick, file.format="newick", check.node.labels="drop")
    expect_equal(labels(newTr), labNew)
    expect_equivalent(sort(edgeLength(newTr)), sort(eLnew))
    expect_equal(nodeType(newTr), nTnew)
})

test_that("check.node.labels='drop' with readNewick", {
    newTr <- readNewick(file=newick, check.node.labels="drop")
    expect_equal(labels(newTr), labNew)
    expect_equivalent(sort(edgeLength(newTr)), sort(eLnew))
    expect_equal(nodeType(newTr), nTnew)
})

test_that("check.node.labels='asdata' with readNCL", {
    newTr <- readNCL(file=newick, file.format="newick", check.node.labels="asdata")
    expect_equal(labels(newTr), labNew)
    expect_equal(tdata(newTr)$labelValues, factor(c(NA, NA, NA, "yy", "xx")))
})

test_that("check.node.labels='asdata' with readNewick", {
    newTr <- readNewick(file=newick, check.node.labels="asdata")
    expect_equal(labels(newTr), labNew)
    expect_equal(tdata(newTr)$labelValues, factor(c(NA, NA, NA, "yy", "xx")))
})

test_that("check.node.labels='keep' with readNCL", {
    labNew[4:5] <- c("yy", "xx")
    newTr <- readNCL(file=newick, file.format="newick", check.node.labels="keep")
    expect_equal(labels(newTr), labNew)
})

test_that("check.node.labels='keep' with readNewick", {
    labNew[4:5] <- c("yy", "xx")
    newTr <- readNewick(file=newick, check.node.labels="keep")
    expect_equal(labels(newTr), labNew)
})

### Test with files where trees don't include all taxa -------------------------
context("Trees that don't contain all the taxa listed in the TAXA block")

test_that("first tree is correct", {
              tr <- readNexus(file = treeSubset)
              expect_equivalent(rootNode(tr[[1]]), 6)
              expect_equivalent(rootNode(tr[[2]]), 6)
              expect_equivalent(rootNode(tr[[3]]), 7)
              expect_equivalent(tipLabels(tr[[1]]), c("porifera", "ctenophora", "cnidaria", "deuterostomia", "protostomia"))
              expect_equivalent(tipLabels(tr[[2]]), c("porifera", "ctenophora", "xeno", "deuterostomia", "protostomia"))
              expect_equivalent(tipLabels(tr[[3]]), c("deuterostomia", "protostomia", "porifera", "ctenophora", "cnidaria", "xeno"))
          }
)

### Test roundtrip with Myrmecus file ------------------------------------------

context("Compare output from rncl read file and phylobase")

test_that("output from rncl::read_nexus_phylo and readNexus match", {
            tr_ape <- rncl::read_nexus_phylo(file = treeDiscDt)
            tr_ph4 <- readNexus(file = treeDiscDt, type = "tree")
            tr_ape <- as(tr_ape, "phylo4")
            expect_equal(edges(tr_ape)[order(edges(tr_ape)[, 1]), ],
                         edges(tr_ph4)[order(edges(tr_ph4)[, 1]), ])
            expect_equal(edgeLength(tr_ape),
                         edgeLength(tr_ph4)[names(edgeLength(tr_ape))])
            expect_equal(labels(tr_ape), labels(tr_ph4))
})
