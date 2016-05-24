Sys.setenv("R_TESTS" = "")
library(eggCounts)
library(testthat)
# simulate egg counts
set.seed(1)
count <- simData1s()
counts <- simData2s()

# test raw counts and calculated epg result same mcmc samples
set.seed(1)
t1 <- stan2mcmc(fec_stan(count[,"obs"], rawCounts=FALSE, CF=50, zeroInflation = TRUE, nburnin=10,nsamples=1e3,thinning=1))
set.seed(1)
t2 <- stan2mcmc(fec_stan(count[,"obs"]/50, rawCounts=TRUE, CF=rep(50,10), zeroInflation = TRUE, nburnin=10,nsamples=1e3,thinning=1))
expect_that(t1, equals(t2))

# test 1-sample model compiles okay
t3<-fec_stan(count[,"obs"],rawCounts=FALSE,CF=50,zeroInflation = TRUE,
         muPrior=list(priorDist = "normal",hyperpars=c(1000,100)),nsamples=1000,
         nburnin=500,thinning=10,nchain=2,ncore=1,verbose=TRUE)
expect_that(t3,is_a("stanfit"))

t4<-fec_stan(count[,"obs"],rawCounts=FALSE,CF=50,zeroInflation = FALSE,
             muPrior=list(priorDist = "normal",hyperpars=c(1000,100)),nsamples=1000,
             nburnin=500,thinning=10,nchain=2,ncore=1,verbose=TRUE)
expect_that(t4,is_a("stanfit"))

# test 2-sample model compiles okay
t5<-fecr_stan(counts[,"masterPre"],counts[,"masterPost"],rawCounts=TRUE,preCF=50,zeroInflation = TRUE,
             muPrior=list(priorDist = "normal",hyperpars=c(1000,100)),nsamples=1000,
             nburnin=500,thinning=10,nchain=2,ncore=1,verbose=TRUE)
expect_that(t5,is_a("stanfit"))

# test 2-sample default models compiles okay
t6<-fecr_stan(counts[,"masterPre"],counts[,"masterPost"],rawCounts=TRUE,paired=TRUE,zeroInflation = TRUE,nsamples=1000,
              nburnin=500,verbose=TRUE)
expect_that(t6,is_a("stanfit"))

t7<-fecr_stan(counts[,"masterPre"],counts[,"masterPost"],rawCounts=TRUE,paired=FALSE,zeroInflation = TRUE,nsamples=1000,
              nburnin=500,verbose=TRUE)
expect_that(t7,is_a("stanfit"))

t8<-fecr_stan(counts[,"masterPre"],counts[,"masterPost"],rawCounts=TRUE,paired=TRUE,zeroInflation = FALSE,nsamples=1000,
              nburnin=500,verbose=TRUE)
expect_that(t8,is_a("stanfit"))

t9<-fecr_stan(counts[,"masterPre"],counts[,"masterPost"],rawCounts=TRUE,paired=FALSE,zeroInflation = FALSE,nsamples=1000,
              nburnin=500,verbose=TRUE)
expect_that(t9,is_a("stanfit"))

# check incorrect inputs
expect_that(fecr_stan(counts[,"masterPre"],counts[,"masterPost"],rawCounts=FALSE,preCF=50),throws_error())
expect_that(fecr_stan(counts[,"obsPre"],counts[,"obsPost"],preCF=50.5), throws_error())

