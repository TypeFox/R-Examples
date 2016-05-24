context('GOF')

test_that('GOF run', {

    #Exploratory
    nfact <- 3
    GOFresult <- GOF(holzinger, nfact)
    GOFresult.outlier <- GOF(holzinger.outlier, nfact)
    expect_equal(as.numeric(GOFresult[1:3]), c(-0.7431551, -0.3786938, -1.2646524),
                 tolerance=1e-5)
    expect_is(GOFresult, 'GOF')
    expect_is(GOFresult.outlier, 'GOF')
    expect_is(plot(GOFresult), 'trellis')
    expect_is(plot(GOFresult.outlier), 'trellis')

    #-------------------------------------------------------------------
    suppressMessages(model <- sem::specifyModel(file='sem-model/sem-model.txt', quiet=TRUE))
    GOFresult <- GOF(holzinger, model)
    GOFresult.outlier <- GOF(holzinger.outlier, model)
    expect_equal(as.numeric(GOFresult[1:3]), c(-4.945440, -1.943843, -3.330193),
                 tolerance=1e-5)
    expect_is(GOFresult, 'GOF')
    expect_is(GOFresult.outlier, 'GOF')
    expect_is(plot(GOFresult), 'trellis')
    expect_is(plot(GOFresult.outlier), 'trellis')

    #-------------------------------------------------------------------
    #Confirmatory with lavaan
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MxdArit + OddWrds
    F3 =~ Boots + Gloves + Hatchts'

    GOFresult <- GOF(holzinger, model, orthogonal=TRUE)
    expect_equal(as.numeric(GOFresult[1:3]), c(-4.984276, -1.952360, -3.352713),
                 tolerance=1e-5)
    GOFresult <- GOF(holzinger.outlier, model, orthogonal=TRUE)
    expect_is(GOFresult, 'GOF')
    expect_is(GOFresult.outlier, 'GOF')
    expect_is(plot(GOFresult), 'trellis')
    expect_is(plot(GOFresult.outlier), 'trellis')
})

test_that('GOF categorical', {
    set.seed(1234)
    cut <- rnorm(ncol(holzinger.outlier), 0 , .25)
    dat <- holzinger.outlier
    for(i in 1:length(cut))
        dat[,i] <- ifelse(holzinger.outlier[,i] < cut[i], 0, 1)

    dat <- as.data.frame(dat)
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MxdArit + OddWrds
    F3 =~ Boots + Gloves + Hatchts'

    GOFresult <- GOF(dat, model, orthogonal=TRUE, ordered=colnames(dat))
    expect_is(GOFresult, 'GOF')
    expect_equal(as.numeric(head(GOFresult)),
                 c(2.392693, -8.217107, -2.637149, -17.303106, 4.298777, -11.997538),
                 tolerance = 1e-3)
})

test_that('GOF mirt', {
    data(LSAT7, package = 'mirt')
    dat <- mirt::expand.table(LSAT7)
    model <- mirt::mirt.model('F = 1-5')
    GOFresult <- GOF(dat, model)
    expect_equal(as.numeric(GOFresult[1:3]), c(0.028003080, 0.001699223, 0.034445923),
                 tolerance=1e-5)
    expect_is(GOFresult, 'GOF')
    expect_is(plot(GOFresult), 'trellis')
})