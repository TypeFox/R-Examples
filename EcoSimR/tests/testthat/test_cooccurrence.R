library(EcoSimR)
context("Co-occurrence null model tests")
tmat <- testMatrix <- ranMatGen(aBetaCol=0.5,bBetaCol=0.5, aBetaRow=0.5,bBetaRow=0.5,numRows=30,numCols=30,
                                mFill=0.25,abun=0,emptyRow=FALSE,emptyCol=FALSE)$m

tmatW0 <- tmat
tmatW0[10,] <- rep(0,dim(tmat)[2])

realData <- as.matrix(dataWiFinches[,2:20])

test_that("vector_sample algorithm works",{
  expect_true(is.vector(vector_sample(speciesData=rbinom(10,1,0.5),weights=runif(10))))
})

test_that("sim1 algorithm works",{
  expect_true(is.matrix(sim1(tmat)))
  expect_true(is.matrix(sim1(realData)))
})

test_that("sim2 algorithm works",{
  expect_true(is.matrix(sim2(tmat)))
  expect_true(is.matrix(sim2(realData)))
  ### Test that row sums are preserved
  expect_equal(apply(sim2(tmat),1,sum),apply(tmat,1,sum))
})

test_that("sim3 algorithm works",{
  expect_true(is.matrix(sim3(tmat)))
  expect_true(is.matrix(sim3(realData)))
  ### Test that row sums are preserved
  expect_equal(apply(sim3(tmat),2,sum),apply(tmat,2,sum))
})


test_that("sim4 algorithm works",{
  expect_true(is.matrix(sim4(tmat)))
  expect_true(is.matrix(sim4(realData)))
  ### Test that row sums are preserved
  expect_equal(apply(sim4(tmat),1,sum),apply(tmat,1,sum))
})

test_that("sim5 algorithm works",{
  expect_true(is.matrix(sim5(tmat)))
  expect_true(is.matrix(sim5(realData)))
  ### Test that row sums are preserved
  expect_equal(apply(sim5(tmat),2,sum),apply(tmat,2,sum))
})
  
test_that("sim6 algorithm works",{
  expect_true(is.matrix(sim6(tmat)))
  expect_true(is.matrix(sim6(realData)))
})

test_that("sim7 algorithm works",{
  expect_true(is.matrix(sim7(tmat)))
  expect_true(is.matrix(sim7(realData)))
})


test_that("sim8 algorithm works",{
  expect_true(is.matrix(sim8(tmat)))
  expect_true(is.matrix(sim8(realData)))
})

test_that("sim10 algorithm works",{
  mopts <- list()
  mopts[["speciesData"]] <- tmat
  mopts[["rowWeights"]] <- runif(dim(tmat)[1])
  mopts[["colWeights"]] <- runif(dim(tmat)[2])
  
  
  expect_true(is.matrix(do.call(sim10,mopts)))
  mopts[["speciesData"]] <- realData
  mopts[["rowWeights"]] <- runif(dim(realData)[1])
  mopts[["colWeights"]] <- runif(dim(realData)[2])
  
  expect_true(is.matrix(do.call(sim10,mopts)))
  expect_true(is.matrix(sim10(tmat)))
})

test_that("species_combo metric works",{
  expect_equal(species_combo(realData),10)
  expect_true(species_combo(tmat) > 0)
})


test_that("checker metric works",{
  expect_equal(checker(realData),91)
  expect_true(is.numeric(checker(tmat)))
})


test_that("c_score metric works",{
  expect_equal(round(c_score(realData)), 4)
  expect_true(is.numeric(c_score(tmat)))
  expect_true(is.numeric(c_score(tmatW0)))
  
})

test_that("c_score_var metric works",{
  expect_equal(round(c_score_var(realData)), 64)
  expect_true(is.numeric(c_score_var(tmat)))
  expect_true(is.numeric(c_score_var(tmatW0)))
  
})

test_that("c_score_skew metric works",{
  expect_equal(round(c_score_skew(realData)), 4)
  expect_true(is.numeric(c_score_skew(tmat)))
  expect_true(is.numeric(c_score_skew(tmatW0)))
  
})

test_that("v_ratio metric works",{
  expect_equal(round(v_ratio(realData)), 1)
  expect_true(is.numeric(v_ratio(tmat)))
  expect_true(is.numeric(v_ratio(tmatW0)))
})

test_that("All combinations of co-occurrence models work",{
  

  
  expect_is(cooc_null_model(dataWiFinches, algo="sim1",metric="species_combo",nRep=10),"coocnullmod")
  expect_is(cooc_null_model(dataWiFinches, algo="sim2",metric="species_combo",nRep=10),"coocnullmod")
  expect_is(cooc_null_model(dataWiFinches, algo="sim3",metric="species_combo",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim4",metric="species_combo",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim5",metric="species_combo",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim6",metric="species_combo",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim7",metric="species_combo",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim8",metric="species_combo",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim10",metric="species_combo",nRep=10,algoOpts=list(rowWeights=runif(dim(dataWiFinches)[1]) ,colWeights=runif(dim(dataWiFinches)[2]-1))),"coocnullmod")
  
  expect_is(cooc_null_model(dataWiFinches, algo="sim1",metric="checker",nRep=10),"coocnullmod")
  expect_is(cooc_null_model(dataWiFinches, algo="sim2",metric="checker",nRep=10),"coocnullmod")
  expect_is(cooc_null_model(dataWiFinches, algo="sim3",metric="checker",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim4",metric="checker",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim5",metric="checker",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim6",metric="checker",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim7",metric="checker",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim8",metric="checker",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim10",metric="checker",nRep=10,algoOpts=list(rowWeights=runif(dim(dataWiFinches)[1]) ,colWeights=runif(dim(dataWiFinches)[2]-1))),"coocnullmod")
  
  expect_is(cooc_null_model(dataWiFinches, algo="sim1",metric="c_score",nRep=10),"coocnullmod")
  expect_is(cooc_null_model(dataWiFinches, algo="sim2",metric="c_score",nRep=10),"coocnullmod")
  expect_is(cooc_null_model(dataWiFinches, algo="sim3",metric="c_score",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim4",metric="c_score",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim5",metric="c_score",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim6",metric="c_score",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim7",metric="c_score",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim8",metric="c_score",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim10",metric="c_score",nRep=10,algoOpts=list(rowWeights=runif(dim(dataWiFinches)[1]) ,colWeights=runif(dim(dataWiFinches)[2]-1))),"coocnullmod")
  
  expect_is(cooc_null_model(dataWiFinches, algo="sim1",metric="c_score_var",nRep=10),"coocnullmod")
  expect_is(cooc_null_model(dataWiFinches, algo="sim2",metric="c_score_var",nRep=10),"coocnullmod")
  expect_is(cooc_null_model(dataWiFinches, algo="sim3",metric="c_score_var",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim4",metric="c_score_var",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim5",metric="c_score_var",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim6",metric="c_score_var",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim7",metric="c_score_var",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim8",metric="c_score_var",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim10",metric="c_score_var",nRep=10,algoOpts=list(rowWeights=runif(dim(dataWiFinches)[1]) ,colWeights=runif(dim(dataWiFinches)[2]-1))),"coocnullmod")
  
  expect_is(cooc_null_model(dataWiFinches, algo="sim1",metric="c_score_skew",nRep=10),"coocnullmod")
  expect_is(cooc_null_model(dataWiFinches, algo="sim2",metric="c_score_skew",nRep=10),"coocnullmod")
  expect_is(cooc_null_model(dataWiFinches, algo="sim3",metric="c_score_skew",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim4",metric="c_score_skew",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim5",metric="c_score_skew",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim6",metric="c_score_skew",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim7",metric="c_score_skew",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim8",metric="c_score_skew",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim10",metric="c_score_skew",nRep=10,algoOpts=list(rowWeights=runif(dim(dataWiFinches)[1]) ,colWeights=runif(dim(dataWiFinches)[2]-1))),"coocnullmod")
  
  expect_is(cooc_null_model(dataWiFinches, algo="sim1",metric="v_ratio",nRep=10),"coocnullmod")
  expect_is(cooc_null_model(dataWiFinches, algo="sim2",metric="v_ratio",nRep=10),"coocnullmod")
  expect_is(cooc_null_model(dataWiFinches, algo="sim3",metric="v_ratio",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim4",metric="v_ratio",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim5",metric="v_ratio",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim6",metric="v_ratio",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim7",metric="v_ratio",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim8",metric="v_ratio",nRep=10),"coocnullmod")  
  expect_is(cooc_null_model(dataWiFinches, algo="sim10",metric="v_ratio",nRep=10,algoOpts=list(rowWeights=runif(dim(dataWiFinches)[1]) ,colWeights=runif(dim(dataWiFinches)[2]-1))),"coocnullmod")
  
  smod <-   cooc_null_model(dataWiFinches, algo="sim8",metric="v_ratio",nRep=10)
  
  
  expect_output(summary(smod),"Metric:  v_ratio")
  expect_true(is.null(plot(smod,type="hist")))
  expect_true(is.null(plot(smod,type="cooc")))
  
})

test_that("all text data frames are handled proprely",{
  textMat <- cbind(letters[1:dim(tmat)[1]],tmat)
  expect_is(cooc_null_model(textMat, algo="sim8",metric="v_ratio",nRep=10),"coocnullmod")  
  
  
})

test_that("Reproduce model works",{
  smod <-   cooc_null_model(dataWiFinches, algo="sim8",metric="v_ratio",nRep=10,saveSeed = T)
  reproduce_model(smod)
  smod2 <- cooc_null_model(dataWiFinches, algo="sim8",metric="v_ratio",nRep=10)
  expect_equal(sum(smod2$Sim-smod$Sim),0)
  smod <-   cooc_null_model(dataWiFinches, algo="sim8",metric="v_ratio",nRep=10)
  expect_error(reproduce_model(smod))
  
  
})

test_that("testing for sim9_single",{
  expect_true(is.matrix(sim9_single()))
  expect_true(is.matrix(sim9_single(tmat)))
  expect_equal(sum(rowSums(tmat)-rowSums(sim9_single(tmat))),0)
  expect_equal(sum(colSums(tmat)-colSums(sim9_single(tmat))),0)
  
})

test_that("testing for sim9",{
  expect_is( cooc_null_model(dataWiFinches, algo="sim9",metric="species_combo",nReps=10,burn_in=10),"coocnullmod")
  expect_is( cooc_null_model(dataWiFinches, algo="sim9",metric="checker",nReps=10,burn_in=10),"coocnullmod") 
  expect_is( cooc_null_model(dataWiFinches, algo="sim9",metric="c_score",nReps=10,burn_in=10),"coocnullmod")
  expect_is( cooc_null_model(dataWiFinches, algo="sim9",metric="c_score_var",nReps=10,burn_in=10),"coocnullmod")
  expect_is( cooc_null_model(dataWiFinches, algo="sim9",metric="c_score_skew",nReps=10,burn_in=10),"coocnullmod")
    expect_is( cooc_null_model(dataWiFinches, algo="sim9",metric="v_ratio",nReps=10,burn_in=10),"coocnullmod")
  
  smod <- cooc_null_model(dataWiFinches, algo="sim9",metric="v_ratio",nReps=10,burn_in=10)
  expect_output(summary(smod),"Algorithm:  sim9")
  expect_true(is.null(plot(smod,type="hist")))
  expect_true(is.null(plot(smod,type="cooc")))
  expect_true(is.null(plot(smod,type="burn_in")))
  
  
})

