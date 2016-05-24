library(testthat)
library(rpf)
library(mvtnorm)

context("chen & thissen 1997")

mxSimplify2Array <- function(x, higher=FALSE) {
	if (higher) {
		stop("higher=TRUE is not implemented. Consider using simplify2array")
	}
  len <- sapply(x, length)
  biggest <- which(len == max(len))[1]
  out <- matrix(NA, nrow=max(len), ncol=length(x))
  for (iter in 1:length(x)) {
	  if (len[iter] == 0) next
    out[1:len[iter],iter] <- x[[iter]]
  }
  colnames(out) <- names(x)
  rownames(out) <- names(x[[biggest]])
  out
}

test_that("ct1997", {
  set.seed(1)
  
  spec <- list()
  spec[1:5] <- rpf.grm(factors=2)
  spec[6]   <- rpf.drm(factors=2)
  gen.param <- mxSimplify2Array(lapply(spec, rpf.rparam, version=1))
  colnames(gen.param) <- paste("i", 1:ncol(gen.param), sep="")
  gen.param[2,] <- c(0,0,.5,.5,1,1)
  
  theta <- rmvnorm(1000, c(0,0), diag(2))
  resp <- rpf.sample(t(theta), spec, gen.param)
  
  # hide latent factor that we don't know about
  tspec <- list()
  tspec[1:length(spec)] <- lapply(spec, rpf.modify, factors=1)
  
  grp <- list(spec=tspec, param=gen.param[-2,], mean=c(0), cov=diag(1), data=resp)
  
  got <- ChenThissen1997(grp)
  #cat(deparse(round(got$pval[!is.na(got$pval)], 2)))
  expect_equal(got$pval[!is.na(got$pval)], 
               c(-1.31, -0.35, 3.15, 5.7, -5.66, 0.54, 3.43, -5.45, -0.57,
                 7.5,  8.11, 1.78, 15.37, 2.79, 5.68), tolerance=.001)
  #cat(deparse(round(got$gamma[!is.na(got$gamma)], 3)))
  expect_equal(got$gamma[!is.na(got$gamma)],
               c(-0.064, -0.002, 0.058, 0.062, -0.197, 0.023, 0.055, -0.008,
                 -0.025, 0.158, 0.13, 0.081, 0.232, 0.02, 0.052), tolerance=.01)
})

test_that("ct1997 2tier", {
  set.seed(1)
  require(rpf)
	spec <- list()
	spec[1:5] <- rpf.drm(factors=3)
	gen.param <- sapply(spec, rpf.rparam)
  gen.param['a2', 1:2] <- 0
  gen.param['a3', 3] <- 0
  gen.param[c('a2','a3'), 4:5] <- 0
  colnames(gen.param) <- paste("i", 1:ncol(gen.param), sep="")
	resp <- rpf.sample(500, spec, gen.param)
	grp <- list(spec=spec, param=gen.param, mean=runif(3, 0, 1), cov=diag(runif(3,1,2)), data=resp)
	slow <- ChenThissen1997(grp, qpoints=13L, qwidth=4, .twotier=FALSE)
	fast <- ChenThissen1997(grp, qpoints=13L, qwidth=4, .twotier=TRUE)
  expect_equal(slow$raw[!is.na(slow$raw)],
               fast$raw[!is.na(fast$raw)], tolerance=.001)
  
  grp$data <- rpf.sample(200, spec, gen.param, mcar=.2)
  fast <- ChenThissen1997(grp, qpoints=13L, qwidth=4)
  got <- fast$pval[,'i1']
  names(got) <- NULL
#  cat(deparse(round(fast$pval[,'i1'],2)))
  expect_equal(got, c(1.49, 1.34, 2.79, -1.54), tolerance=.1)
})

mxSimplify2Array <- function(x) {
  par <- x
  len <- sapply(par, length)
  biggest <- which(len == max(len))[1]
  out <- matrix(NA, nrow=max(len), ncol=length(par))
  for (x in 1:length(par)) {
    out[1:len[x],x] <- par[[x]]
  }
  rownames(out) <- names(par[[biggest]])
  out
}

test_that("ct1997 permutations", {
  require(rpf)
  set.seed(1)
  grp <- list(spec=list())
  grp$spec[1:10] <- lapply(sample.int(6, 10, replace=TRUE), function(o) rpf.grm(outcomes=1+o))
  grp$param <- mxSimplify2Array(lapply(grp$spec, rpf.rparam))
  colnames(grp$param) <- paste("i", 1:10, sep="")
  grp$mean <- 0
  grp$cov <- diag(1)
  grp$free <- !is.na(grp$param) & grp$param != 0
  grp$data <- rpf.sample(500, grp=grp)
  grp$data <- grp$data[,colnames(grp$data)[sample.int(10)]]
  
  ChenThissen1997(grp, inames = colnames(grp$param)[sample.int(10, 9)])
  # will stop if something is wrong
  
  grp1 <- grp
  grp1$data <- grp$data[1:250,]
  grp2 <- grp
  grp2$data <- grp$data[251:500,]

  r1 <- ChenThissen1997(grp)
  r2 <- ChenThissen1997(grp1) + ChenThissen1997(grp2)
  expect_equal(r1$pval, r2$pval)
  expect_equal(r1$gamma, r2$gamma)
})

drawRandomProportion <- function(expected) {
  total <- sum(expected)
  prob <- expected / total
  sim <- rep(NA, length(expected))
  rowSim <- sample.int(length(expected), size=total, prob=prob, replace=TRUE)
  sim <- tabulate(rowSim, length(expected))
  sim
}

if (0) {
  spec <- list()
  spec[1:6] <- rpf.grm(factors=1)
  gen.param <- sapply(spec, rpf.rparam)
  grp <- list(spec=spec, param=gen.param, mean=c(0), cov=diag(1))

  pair <- c(1L,2L)
  numPeople <- 185
  E <- (numPeople * pairwiseExpected(grp, pair))

  ms <- pairwiseItemDistribution(grp, pair)
  hist(ms)
  quantile(ms, c(.95, .99))   # 25 38
  
  trials <- 10000
  ms <- rep(NA, trials)
  for (tx in 1:trials) {
    O <- drawPairwiseSample(grp, pair, numPeople)
    ms[tx] <- sum((c(E) - O)^2)
  }

  hist(ms)
  quantile(ms, c(.95, .99))   # 25 38
}

pearson.gof <- function(observed, expected, df) {
  x2 <- sum((observed - expected)^2/expected)
  if (missing(df)) {
    df <- (dim(observed)[1]-1) * (dim(observed)[2]-1)
  }
  pchisq(x2, df, lower.tail=FALSE)
}

if (0) {
  spec <- list()
  spec[1:6] <- rpf.grm(factors=1)
  gen.param <- sapply(spec, rpf.rparam)
  grp <- list(spec=spec, param=gen.param, mean=c(0), cov=diag(1))

  pair <- c(1L,2L)
  numPeople <- 500
  E <- (numPeople * pairwiseExpected(grp, pair))

  trials <- 50
  got <- expand.grid(trial=1:trials, method=c("pearson","rms"), pval=NA)
  for (rep in 1:trials) {
    O <- drawPairwiseSample(grp, pair, 50, qpts=11L, qwidth=4)  # doesn't work!
    O <- O * numPeople / 50
    got[got$trial==rep,'pval'] <- c(pearson.gof(O, E),
                                    pairwiseItemTest(grp, pair, O, qpts=11, qwidth=4))
#                                    ptw2011.gof.test(O, E)
   print(rep)
  }
  #
  require(ggplot2)
 mask <- got$method=="rms"
  #mask <- rep(TRUE, nrow(got))
 pval <- got[mask,'pval']
 got[mask,'pval'] <- 1 / (1+exp(-(logit(pval) - 2.8)))
  
  tbl <- expand.grid(alpha=seq(0,1,.01), method=c("pearson","rms"), emp=NA)
  tbl[tbl$method=="rms",'emp'] <-
    Vectorize(function(x) sum(got[got$method=="rms",'pval'] < x))(seq(0,1,.01))
  tbl[tbl$method=="pearson",'emp'] <-
    Vectorize(function(x) sum(got[got$method=="pearson",'pval'] < x))(seq(0,1,.01))
  tbl$emp <- tbl$emp / trials
  
  ggplot(tbl, aes(alpha, emp, color=method)) + geom_line() +
    geom_abline(slope=1, color="yellow")+ coord_fixed()
}

