test_that("test_pv2sc", {
    pfm <- matrix(c(3, 5, 4, 2, 7, 0, 3, 4, 9, 1, 1, 3, 3, 6, 
        4, 1, 11, 0, 3, 0, 11, 0, 2, 1, 11, 0, 2, 1, 3, 3, 2, 
        6, 4, 1, 8, 1, 3, 4, 6, 1, 8, 5, 1, 0, 8, 1, 4, 1, 9, 
        0, 2, 3, 9, 5, 0, 0, 11, 0, 3, 0, 2, 7, 0, 5), nrow = 4, 
        dimnames = list(c("A", "C", "G", "T")))
    bg <- c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
    pvalue <- 0.00001
    type <- "PFM"
    score <- TFMpv2sc(pfm, pvalue, bg, type)
    expect_equal(score, 8.773708, tolerance=0.01)
})

