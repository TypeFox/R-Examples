context('gCD')

test_that('gCD run', {

    #Exploratory
    nfact <- 3
    gCDresult <- gCD(holzinger, nfact)
    gCDresult.outlier <- gCD(holzinger.outlier, nfact)
    expect_is(gCDresult, 'gCD')
    expect_is(plot(gCDresult), 'trellis')
    expect_is(gCDresult.outlier, 'gCD')
    expect_is(plot(gCDresult.outlier), 'trellis')
    expect_equal(as.numeric(gCDresult$gCD[1:3]), c(1.646295e-04, 6.672447e-05, 1.033961e-04), tolerance=1e-5)

    #-------------------------------------------------------------------
    suppressMessages(model <- sem::specifyModel(file='sem-model/sem-model.txt', quiet=TRUE))
    gCDresult <- gCD(holzinger, model)
    gCDresult.outlier <- suppressWarnings(gCD(holzinger.outlier, model))
    expect_is(gCDresult, 'gCD')
    expect_is(plot(gCDresult), 'trellis')
    expect_is(gCDresult.outlier, 'gCD')
    expect_is(plot(gCDresult.outlier), 'trellis')
    expect_equal(as.numeric(gCDresult$gCD[1:3]), c(1.339090e-05, 4.191384e-06, 3.585573e-05), tolerance=1e-5)

    #-------------------------------------------------------------------
    #Confirmatory with lavaan
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MxdArit + OddWrds
    F3 =~ Boots + Gloves + Hatchts'

    gCDresult <- gCD(holzinger, model, orthogonal=TRUE)
    gCDresult.outlier <- gCD(holzinger.outlier, model, orthogonal=TRUE)
    expect_equal(as.numeric(gCDresult.outlier$gCD[1:3]), c(1.040270e-02, 1.746456e-05, 8.223349e-05),
                 tolerance = 1e-5)
    expect_is(gCDresult, 'gCD')
    expect_is(plot(gCDresult), 'trellis')
    expect_is(gCDresult.outlier, 'gCD')
    expect_is(plot(gCDresult.outlier), 'trellis')
})

test_that('gCD categorical', {
    set.seed(1234)
    cut <- rnorm(ncol(holzinger.outlier), 0 , .25)
    dat <- holzinger.outlier
    for(i in 1:length(cut))
        dat[,i] <- ifelse(holzinger.outlier[,i] < cut[i], 0, 1)

    dat <- as.data.frame(dat)
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MxdArit + OddWrds
    F3 =~ Boots + Gloves + Hatchts'

    gCDresult <- suppressWarnings(gCD(dat, model, orthogonal=TRUE, ordered=colnames(dat)))
    expect_is(gCDresult, 'gCD')
    expect_equal(as.numeric(head(gCDresult$gCD)), c(0.0007, 0.0005, 0.0069, 0.0005, 0.0012, 0.0038),
                 tolerance = 1e-2)

#     model <- mirt::mirt.model('F1 = 1-3
#                             F2 = 4-6
#                             F3 = 7-9')
#     gCDresult2 <- suppressWarnings(gCD(dat, model))

})

