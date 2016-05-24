context("MSM baseline covariates")

test_that("constant W in baseline equivalent to constant summary measure", {
  n <- 200
  niter <- 3
  for (i in 1:niter) {
    k <- runif(1, min=2, max=5)
    for (survivalOutcome in c(T, F)) {
      W <- k 
      A1 <- rbinom(n, 1, 0.5)
      A2 <- rbinom(n, 1, 0.5)
      Y1 <- rbinom(n, 1, 0.1)
      A3 <- rbinom(n, 1, 0.5)
      Y2 <- rexpit(-5 + W + A1 + A2 + A3 + rnorm(n))
      if (survivalOutcome) {
        Y2 <- as.numeric(Y2 | Y1)
        A3[Y1==1] <- NA
      }
      
      regimes <- array(dim=c(n, 3, 2))
      regimes[,,1] <- 1
      regimes[,,2] <- 0
      summary.measures <- array(c(k, k, 5, 6, k, k, 7, 8), dim=c(2, 2, 2)) 
      colnames(summary.measures) <- c("k", "junk")
      
      data <- data.frame(W, A1, A2, Y1, A3, Y2)
      x1 <- ltmleMSM(data, Anodes=c("A1", "A2", "A3"), Ynodes=c("Y1", "Y2"), regimes=regimes, summary.measures=summary.measures, working.msm="Y~-1 + k + junk", estimate.time=F, survivalOutcome=survivalOutcome, final.Ynodes=c("Y1", "Y2"))
      
      x2 <- ltmleMSM(data, Anodes=c("A1", "A2", "A3"), Ynodes=c("Y1", "Y2"), regimes=regimes, summary.measures=summary.measures, working.msm="Y~- 1+ W + junk", estimate.time=F, survivalOutcome=survivalOutcome, final.Ynodes=c("Y1", "Y2"))
      
      summ <- function(x) summary(x)$cmat[, 1:2] #estimate and std error
      
      expect_equal(summ(x1), summ(x2), tolerance=0.0001, scale=1, check.attributes=FALSE)
    }
  }
})

test_that("error if baseline and summary have the same name or covariate not found", {
  n <- 100
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  A1 <- rbinom(n, 1, 0.5)
  L <- rnorm(n)
  A2 <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)
  regimes <- array(dim=c(n, 2, 2))
  regimes[,,1] <- 1
  regimes[,,2] <- 0
  summary.measures <- array(1:4, dim=c(2, 2, 1)) 
  colnames(summary.measures) <- c("W2", "junk")
  
  data <- data.frame(W1, W2, A1, L, A2, Y)
  expect_error(ltmleMSM(data, Anodes=c("A1", "A2"), Lnodes="L", Ynodes="Y", regimes=regimes, summary.measures=summary.measures, working.msm="Y ~ W2 + junk", estimate.time=F), "Baseline covariate columns of data and columns of summary.measures may not have the same name")
  
  colnames(summary.measures) <- c("Z", "junk")
  expect_error(ltmleMSM(data, Anodes=c("A1", "A2"),  Lnodes="L", Ynodes="Y", regimes=regimes, summary.measures=summary.measures, working.msm="Y ~ W3 + junk", estimate.time=F), "All right hand side variables in working.msm must be either column names of summary.measures or column names of baseline covariates")
  
  colnames(summary.measures) <- c("L", "junk") #not a baseline covariate
  expect_error(ltmleMSM(data, Anodes=c("A1", "A2"),  Lnodes="L", Ynodes="Y", regimes=regimes, summary.measures=summary.measures, working.msm="Y ~ W3 + junk", estimate.time=F), "All right hand side variables in working.msm must be either column names of summary.measures or column names of baseline covariates")
  
})

test_that("ltmleMSM runs and returns beta with correct names for baseline only, summary only, or both", {
  data(sampleDataForLtmleMSM)
  Anodes <- grep("^A", names(sampleDataForLtmleMSM$data))
  Lnodes <- grep("^CD4", names(sampleDataForLtmleMSM$data))[-1]
  Ynodes <- grep("^Y", names(sampleDataForLtmleMSM$data))
  
  suppressWarnings(x <- ltmleMSM(data=sampleDataForLtmleMSM$data, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, regimes=sampleDataForLtmleMSM$regimes, summary.measures=sampleDataForLtmleMSM$summary.measures, final.Ynodes=Ynodes, working.msm="Y ~ age * male", estimate.time=FALSE, survivalOutcome=T)) #baseline only
  expect_equal(names(x$beta), c("(Intercept)", "age", "male", "age:male"))
  
  suppressWarnings(x <- ltmleMSM(data=sampleDataForLtmleMSM$data, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, regimes=sampleDataForLtmleMSM$regimes, summary.measures=sampleDataForLtmleMSM$summary.measures, final.Ynodes=Ynodes, working.msm="Y ~ age * switch.time", estimate.time=FALSE, survivalOutcome=T)) #baseline and summary
  expect_equal(names(x$beta), c("(Intercept)", "age", "switch.time", "age:switch.time"))
  
  suppressWarnings(x <- ltmleMSM(data=sampleDataForLtmleMSM$data, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, regimes=sampleDataForLtmleMSM$regimes, summary.measures=sampleDataForLtmleMSM$summary.measures, final.Ynodes=Ynodes, working.msm="Y ~ time + switch.time", estimate.time=FALSE, survivalOutcome=T)) #summary only
  expect_equal(names(x$beta), c("(Intercept)", "time", "switch.time"))
})
