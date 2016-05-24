test_that("IJROI format is read correctly", {
    dataset <- file.path(system.file(package = "retistruct"), "extdata", "smi32")
    r <- retistruct.read.dataset(dataset)
    ## Test that points are read in correctly 
    P <- as.matrix(read.csv(file.path(dataset, "P.csv")))
    attr(P, "dimnames") <- NULL
    expect_that(r$P, equals(P))
    ## Test that optic disc is read in correctly 
    od <- rbind(c(354, 336),
                c(354, 347),
                c(359, 352),
                c(364, 347),
                c(364, 338),
                c(360, 333))
    colnames(od) <- c("x", "y")
    expect_that(r$Ss$OD, equals(od))
})
