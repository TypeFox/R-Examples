context("accuracy-of-stats")

test_that("Stats are accurately calculated", { 
    #setup
    g <- read.genepop("test.gen", ncode=3)
    dstats <- lapply(diff_stats(g), round, 3)
    dstats_loc <- dstats$per.loc
    dstats_glob <- dstats$global
    #stats per loci are accurate
    expect_equivalent(dstats_loc[,"D"], c(0.517, -0.056))
    expect_equivalent(dstats_loc[,"Gprime_st"], c(0.761, -0.111))
    expect_equivalent(dstats_loc[,"Gst"], c(0.337, -0.026))
    expect_equivalent(dstats_loc[,"Hs"], c(0.337, 0.526))
    expect_equivalent(dstats_loc[,"Ht"], c(0.508, .513))

    #gloabal stats are accurate
    expect_equivalent(dstats_glob["D_mean"], NA_real_ ) 
    expect_equivalent(dstats_glob["D_het"], 0.279) 
    expect_equivalent(dstats_glob["Gprime_st"], 0.472) 
    expect_equivalent(dstats_glob["Gst_est"], 0.155) 
    expect_equivalent(dstats_glob["Hs"], 0.432) 
    expect_equivalent(dstats_glob["Ht"], 0.511) 

    #stand alone functions
    dj <- lapply(D_Jost(g), round, 3)
    expect_equivalent(dj$per.locus, c(0.517, -0.056))
    expect_equivalent(dj$global.het, 0.279)
    expect_equivalent(dj$global.harm_mean, NA_real_)

    gh <- lapply(Gst_Hedrick(g), round, 3)
    expect_equivalent(gh$per.locus, c(0.761, -0.111))
    expect_equivalent(gh$global, 0.472)

    gn <- lapply(Gst_Nei(g), round, 3)
    expect_equivalent(gn$per.locus, c(0.337, -0.026))
    expect_equivalent(gn$global, 0.155)
})
