context("EWMA")

test_that("EWMA immediate stop",{
    expect_equal(c(spcadjust:::ARL_EWMA_Markovapprox(0.00000001,pobs=pnorm,lambda=0.5)),
                 1,tolerance=1e-3)
})

test_that("EWMA ARL c=3.09",{
 ## library(spc)
 ## xewma.arl(mu=0,l=1,c=3.09,sided="two")
 expect_equal(c(spcadjust:::ARL_EWMA_Markovapprox(3.09,pobs=pnorm,lambda=1)),
              499.6091,tolerance=1e-3)
})

test_that("EWMA ARL lambda=0.5",{
 lambda=0.5
 L=3
 c=L*sqrt(lambda/(2-lambda))
 ## library(spc)
 ## xewma.arl(mu=0,l=lambda,c=L,sided="two")
 expect_equal(c(spcadjust:::ARL_EWMA_Markovapprox(c,gridpoints=100,pobs=pnorm,lambda=lambda)),
              397.4608,tolerance=0.1)
})

test_that("EWMA ARL 5000",{
 ARL=500
 lambda=0.5
 c <- spcadjust:::calibrateARL_Markovapprox(ARL=ARL,f=ARL_EWMA_Markovapprox,pobs=pnorm,lambda=lambda)
 expect_equal(c(spcadjust:::ARL_EWMA_Markovapprox(c,gridpoints=100,pobs=pnorm,lambda=lambda)),500,tolerance=1e-1)

 ## library(spc)
 ## L <- xewma.crit(l=lambda,L0=ARL,sided="two")
 ## L*sqrt(lambda/(2-lambda))
 expect_less_than(abs(c-1.773076),1e-3)
})

test_that("EWMA n=200, lambda=0.5, hitprob 0.1",{
 hprob=0.1
 n=200
 lambda=0.5
 c <- spcadjust:::calibratehitprob_Markovapprox(hprob=hprob,n=n,f=hitprob_EWMA_Markovapprox,pobs=pnorm,lambda=lambda)
 expect_equal(c(spcadjust:::hitprob_EWMA_Markovapprox(c,n=n,gridpoints=100,pobs=pnorm,lambda=lambda)),
              0.1,tolerance=1e-5)
})
