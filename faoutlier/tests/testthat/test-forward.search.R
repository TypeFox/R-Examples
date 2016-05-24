context('forward.search')

setCluster()

test_that('forward.search run', {    
    set.seed(1234)
    
    #Exploratory
    nfact <- 3
    FS <- forward.search(holzinger, nfact, print.messages = FALSE)
    expect_is(FS, 'forward.search')
    expect_is(plot(FS), 'trellis')
    expect_equal(FS$GOF[c(1, length(FS$GOF))], c(3.717805, 26.293829), tolerance = 1e-5)
    
    #-------------------------------------------------------------------    
    suppressMessages(model <- sem::specifyModel(file='sem-model/sem-model.txt', quiet=TRUE))    
    FS.outlier <- suppressWarnings(forward.search(holzinger.outlier, model, print.messages = FALSE))
    expect_is(FS.outlier, 'forward.search')
    expect_is(plot(FS.outlier), 'trellis')
    expect_equal(FS.outlier$GOF[c(1, length(FS.outlier$GOF))], 
                 c(34.64542, 161.01480), tolerance = 1e-5)
    
    #---- lavaan
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MxdArit + OddWrds
    F3 =~ Boots + Gloves + Hatchts'
    FS.outlier <- forward.search(holzinger.outlier, model, print.messages = FALSE, n.subsets = 200,
                                 p.base = .7)
    expect_is(FS.outlier, 'forward.search')
    expect_is(plot(FS.outlier), 'trellis')
    expect_equal(FS.outlier$GOF[c(1, length(FS.outlier$GOF))], 
                 c(26.96743, 111.39844), tolerance = 1e-5)
    
})