context('LD')

test_that('LD run', {

    #Exploratory
    nfact <- 3
    LDresult <- LD(holzinger, nfact)
    LDresult.outlier <- LD(holzinger.outlier, nfact)
    expect_equal(as.numeric(LDresult[1:3]), c(0.5813941, 0.4363624, 0.1987451),
                 tolerance=1e-5)
    expect_is(LDresult, 'LD')
    expect_is(LDresult.outlier, 'LD')
    expect_is(plot(LDresult), 'trellis')
    expect_is(plot(LDresult.outlier), 'trellis')

    #-------------------------------------------------------------------
    suppressMessages(model <- sem::specifyModel(file='sem-model/sem-model.txt', quiet=TRUE))
    LDresult <- LD(holzinger, model)
    LDresult.outlier <- LD(holzinger.outlier, model)
    expect_equal(as.numeric(LDresult[1:3]), c(-4.945440, -1.943843, -3.330193),
                 tolerance=1e-5)
    expect_is(LDresult, 'LD')
    expect_is(LDresult.outlier, 'LD')
    expect_is(plot(LDresult), 'trellis')
    expect_is(plot(LDresult.outlier), 'trellis')

    #-------------------------------------------------------------------
    #Confirmatory with lavaan
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MxdArit + OddWrds
    F3 =~ Boots + Gloves + Hatchts'

    LDresult <- LD(holzinger, model, orthogonal=TRUE)
    expect_equal(as.numeric(LDresult[1:3]), c(-21.65268, -16.08441, -28.68881),
                 tolerance=1e-5)
    LDresult.outlier <- LD(holzinger.outlier, model, orthogonal=TRUE)
    expect_is(LDresult, 'LD')
    expect_is(LDresult.outlier, 'LD')
    expect_is(plot(LDresult), 'trellis')
    expect_is(plot(LDresult.outlier), 'trellis')
})

test_that('LD categorical', {
    data(LSAT7, package = 'mirt')
    dat <- mirt::expand.table(LSAT7)
    model <- mirt::mirt.model('F = 1-5')
    LDresult <- LD(dat, model)
    expect_equal(as.numeric(LDresult[1:3]), c(-9.233593, -8.005580, -10.857399),
                 tolerance=1e-5)
    expect_is(LDresult, 'LD')
    expect_is(plot(LDresult), 'trellis')
})