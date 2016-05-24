require(testthat)
require(pez)
require(picante)
data(phylocom)

test_that("PA comm dist", {
    cc <- comparative.comm(phylocom$phy,
                           ifelse(phylocom$sample, 1, 0),
                           phylocom$traits, warn=FALSE)
    distByPez  <- comm.dist(cc)
    expect_that(distByPez[69] , equals(1))
    expect_that(distByPez[256], equals(0.5))
    expect_that(distByPez[178], equals(1/3))
    
    distByHand <- 1-as.dist(crossprod(cc$comm)/(dim(cc$comm)[1] - crossprod(1-cc$comm)))
    expect_that(distByPez, is_equivalent_to(distByHand))
})

## TODO:  non-PA dist matrix / overlap warning

test_that("func dist", {
    cc <- comparative.comm(phylocom$phy, phylocom$sample, warn=FALSE)
    expect_error(traits.dist(cc))
    ## TODO:  non-error case
})

test_that("func phylo dist", {
    cc <- comparative.comm(phylocom$phy, phylocom$sample, phylocom$traits, warn=FALSE)
    expect_that(funct.phylo.dist(cc, 0.5, 2)[1:6], equals(c(0.19961604238199601169, 0.52213831399923527066, 0.60980426435637358207, 0.56108118570624387900, 0.57849689828881134535, 0.62785217564210749064)))
    aVec <- seq(0, 1, by = 0.2)
    pVec <- 0:4
    PDist <- phylo.dist(cc)
    FDist <- traits.dist(cc)
    for(a in aVec) {
        for(p in pVec) {
            distByHand <- (a*(PDist/max(PDist))^p +
                           (1-a)*(FDist/max(FDist))^p)^(1/p)
            expect_that(funct.phylo.dist(cc, a, p), equals(distByHand))
        }
    }
    PDist[] <- PDist/max(PDist)
    FDist[] <- FDist/max(FDist)
    expect_that(funct.phylo.dist(cc, 1, 2), is_equivalent_to(PDist))
    expect_that(funct.phylo.dist(cc, 0, 2), is_equivalent_to(FDist))
})
