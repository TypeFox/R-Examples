library(testthat)
library(rpf)
library(mirt)
#options(error = utils::recover)

context("mirt ICC")

myseed <- as.integer(runif(1) * 1e7)
#print(paste("set.seed =",myseed))
set.seed(myseed)

i.count <- 5
spec <- list()

test_that("3PL", {
  spec[1:i.count] <- rpf.drm()
  data <- rpf.sample(100, spec)
  data <- simplify2array(lapply(data, unclass)) - 1

  suppressWarnings(fit <- mirt(data, 1, rep('3PL',i.count), D=1,
                               verbose=FALSE, technical=list(NCYCLES=1)))

  for (ix in 1:i.count) {
    ii <- extract.item(fit, ix)
    expect_equal(c(t(probtrace(ii, c(-1,0,1)))),
                 c(rpf.prob(spec[[1]], c(ii@par[1:2], ii@par[3:4]), c(-1,0,1))))
  }
})

spec[1:i.count] <- rpf.grm(outcomes=3, multidimensional=TRUE)

data <- rpf.sample(100, spec)
data <- simplify2array(lapply(data, unclass)) - 1

test_that("GRM", {
  suppressWarnings(fit <- mirt(data, 1, rep('graded',i.count), D=1,
                               verbose=FALSE, technical=list(NCYCLES=1)))

  for (ix in 1:i.count) {
    ii <- extract.item(fit, ix)
    if (length(ii@par) < 3) {
                                        # mirt can lose a category if it is not represented in the data
                                        #print(data[,ix])
      next
    }
    expect_equal(c(t(probtrace(ii, c(-1,0,1)))),
                 c(rpf.prob(spec[[1]], ii@par[1:3], c(-1,0,1))))
  }
})

test_that("nominal", {
  spec[1:i.count] <- rpf.nrm(outcomes=3,
                             T.a=rbind(0, diag(2)),
                             T.c=rbind(0, diag(2)))

  suppressWarnings(fit <- mirt(data, 1, rep('nominal',i.count),
                               verbose=FALSE, D=1, technical=list(NCYCLES=1)))
  for (ix in 1:i.count) {
    ii <- extract.item(fit, ix)
    if (length(ii@par) < 7) next
    expect_equal(c(t(probtrace(ii, c(-1,0,1)))),
                 c(rpf.prob(spec[[1]], ii@par[c(1, 3:4, 6:7)], c(-1,0,1))))
    expect_equal(c(log(t(probtrace(ii, c(-1,0,1))))),
                 c(rpf.logprob(spec[[1]], ii@par[c(1, 3:4, 6:7)], c(-1,0,1))))
  }
})
