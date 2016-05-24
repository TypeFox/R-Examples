test_that("test_sc2pv", {
    pfm <- matrix(c(3, 5, 4, 2, 7, 0, 3, 4, 9, 1, 1, 3, 3, 6, 
        4, 1, 11, 0, 3, 0, 11, 0, 2, 1, 11, 0, 2, 1, 3, 3, 2, 
        6, 4, 1, 8, 1, 3, 4, 6, 1, 8, 5, 1, 0, 8, 1, 4, 1, 9, 
        0, 2, 3, 9, 5, 0, 0, 11, 0, 3, 0, 2, 7, 0, 5), nrow = 4, 
        dimnames = list(c("A", "C", "G", "T")))
    bg <- c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
    score <- 8.77
    type <- "PFM"
    pvalue <- TFMsc2pv(pfm, score, bg, type)
    expect_equal(pvalue, 0.00001007156, tolerance=1e-5)
})

