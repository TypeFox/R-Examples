#options(error = browser)
library(testthat)
library(rpf)

context("ot2000")

test_that("simple case", {
  require(rpf)
  set.seed(1)

  spec <- list()
  spec[1:3] <- rpf.drm()

  gen.p <- matrix(c(1,0,0,1,
                    1,-1,0,1,
                    1,1,0,1), ncol=3, nrow=4)  # note byrow=FALSE
  gen.p[3:4,] <- logit(gen.p[3:4,])
  data <- rpf.sample(200, spec, gen.p)
                                        #  write.csv(data, "fit-test.csv", row.names=FALSE, quote=FALSE)

  param <- matrix(c(.71,.03,0, 1,
                    1.53,-.85, 0, 1,
                    1.1,.9,0, 1), ncol=3, nrow=4)
  param[3:4,] <- logit(param[3:4,])
  
  colnames(param) <- paste("i", 1:3, sep="")
  grp <- list(spec=spec, param=param, data=data, mean=0, cov=matrix(1,1,1))

  got <- SitemFit(grp, method="pearson")
  
  tbl <- round(t(sapply(got, function(row) c(stat=row$statistic, df=row$df, p=row$pval))),2)
  expect_equal(sum(tbl[,'df']), 9)  # not sure TODO
  expect_equal(sum(tbl[,'stat']), 22.55)
})

test_that("orlando-thissen-2000", {
  require(rpf)
  set.seed(7)
  grp <- list(spec=list())
  grp$spec[1:20] <- rpf.grm()
  grp$param <- sapply(grp$spec, rpf.rparam, version=1L)
  colnames(grp$param) <- paste("i", 1:20, sep="")
  grp$mean <- 0
  grp$cov <- diag(1)
  grp$free <- grp$param != 0
  grp$data <- rpf.sample(500, grp=grp)
  
  got <- SitemFit(grp, method="pearson", omit=1)
  expect_true(all(sapply(got, function(ii) is.null(ii$omitted))))
  stat <- sapply(got, function(x) x$statistic)
  names(stat) <- NULL
  Estat <- c(10.8, 13.38, 19.47, 15.34, 13.62, 6.52, 19.91, 8.19,  8.03, 11.73, 12.84,
             10.4, 5.08, 9.37, 10.37, 11.74, 20.67, 16.34,  13.5, 16.91)
  expect_equal(stat, Estat, tolerance=.01)
  
  E1orig <- structure(c(2, 3.99, 12.95, 21.82, 21.66, 25.29, 30.5, 32.29,  35.09, 24.41, 34.33, 15.74,
              18.61, 14.44, 10.79, 5.72, 3.82,  1.5, 0.23, 0.12, 0, 0.01, 0.05, 0.18, 0.34,
              0.71, 1.5, 2.71,  4.91, 5.59, 12.67, 9.26, 17.39, 21.56, 26.21, 23.28, 27.18,
              19.5,  5.77, 5.88), .Dim = c(20L, 2L))
  oexp <- got[[1]]$orig.expected
  dimnames(oexp) <- NULL
  expect_equal(oexp, E1orig, tolerance=.01)
  
  E1 <- structure(c(2, 3.99, 12.95, 21.82, 21.66, 25.29, 30.5, 32.29,  35.09, 24.41, 34.33,
                    15.74, 18.61, 14.44, 10.79, 5.72, 3.82,  1.85, 0, 0, 0, 0, 0, 0, 0, 1.3,
                    1.5, 2.71, 4.91, 5.59, 12.67,  9.26, 17.39, 21.56, 26.21, 23.28, 27.18,
                    19.5, 5.77, 5.88), .Dim = c(20L,  2L))
  mask <- !is.na(got[[1]]$expected)
  expect_equal(got[[1]]$expected[mask], E1[mask], tolerance=.01)
  
  expect_equal(got[[1]]$pval, -0.984, tolerance=.01)

  got <- SitemFit(grp, method="pearson", alt=TRUE)
  stat <- sapply(got, function(x) x$statistic)
  names(stat) <- NULL
  #cat(deparse(round(stat, 2)))
  Estat <- c(16.6, 13.68, 13.19, 14.03, 19.86, 11.02, 29.92, 6.95, 18.31,
             14.78, 10.06, 9.27, 5.67, 9.3, 14.75, 18.03, 20.18, 20.19, 15.89,  11.57)
  expect_equal(stat, Estat, tolerance=.01)
})

test_that("fit w/ mcar", {
  require(rpf)
  require(testthat)
  set.seed(7)
  grp <- list(spec=list(), qwidth=5, qpoints=31)
  grp$spec[1:20] <- rpf.grm()
  grp$param <- sapply(grp$spec, rpf.rparam, version=1L)
  colnames(grp$param) <- paste("i", 1:20, sep="")
  grp$free <- grp$param != 0
  grp$data <- rpf.sample(500, grp=grp, mcar=.1)

  got <- sumScoreEAPTest(omitMostMissing(grp, 3L))
  expect_equal(got$n, 101L)
  expect_equal(got$rms.p, -1.87, tolerance=.01)
  expect_equal(got$pearson.df, 16L)
  expect_equal(got$pearson.p, -1.46, tolerance=.01)

  grp1 <- grp
  grp1$data <- grp$data[1:250,]
  grp2 <- grp
  grp2$data <- grp$data[251:500,]
  e1 <- sumScoreEAPTest(grp1) + sumScoreEAPTest(grp2)
  e2 <- sumScoreEAPTest(grp)
  chk <- c('n','pearson.df', 'pearson.chisq', 'pearson.p')
  expect_equal(unlist(e1[chk]), unlist(e2[chk]))
  
  got <- SitemFit(grp, omit=2L)
  stat <- sapply(got, function(x) x$statistic)
  names(stat) <- NULL
  Estat <- c(6.68, 6.87, 9.6, 14.31, 11.74, 11.63, 15.05, 2.08, 6.16, 6.8,
             10.28, 2.37, 6.86, 6.47, 10.02, 11.48, 13.96, 9.96, 8.24, 11.77 )
  expect_equal(stat, Estat, tolerance=.01)
  
  e1 <- SitemFit(grp1) + SitemFit(grp2)
  e2 <- SitemFit(grp)
  expect_equal(sapply(e1, function(ii) ii$statistic),
               sapply(e2, function(ii) ii$statistic))
})

test_that("2tier fit", {
  set.seed(1)
  require(rpf)
  numItems <- 6
  spec <- list()
  spec[1:numItems] <- rpf.drm(factors=3)
  param <- sapply(spec, rpf.rparam, version=1)
  gsize <- numItems/3
  for (gx in 0:2) {
    if (gx != 1) {
      param['a2', seq(gx * gsize+1, (gx+1)*gsize)] <- 0
    }
    if (gx != 2) {
      param['a3', seq(gx * gsize+1, (gx+1)*gsize)] <- 0
    }
  }
  grp <- list(spec=spec, param=param, mean=runif(3, -1, 1), cov=diag(runif(3,.5,2)))
  grp$data <- rpf.sample(500, grp=grp)
  colnames(grp$param) <- colnames(grp$data)

  got <- SitemFit(grp, qwidth=2, qpoints=5L, .twotier=FALSE)
  tt <- SitemFit(grp, qwidth=2, qpoints=5L, .twotier=TRUE)
  expect_equal(sapply(got, function(x) x$statistic),
               sapply(tt, function(x) x$statistic), .001)

  got <- SitemFit(grp, qwidth=2, qpoints=5L, .twotier=FALSE, alt = TRUE)
  tt <- SitemFit(grp, qwidth=2, qpoints=5L, .twotier=TRUE, alt=TRUE)
  expect_equal(sapply(got, function(x) x$statistic),
               sapply(tt, function(x) x$statistic), .001)
})

if (0) {
  library(mirt)
  dat <- expand.table(LSAT6)
  model <- mirt.model('F = 1-5
                    CONSTRAIN = (1-5, a1)')
  (mod <- mirt(dat, model))
  itemfit(mod, X2=TRUE, S_X2.tables=TRUE)$E[[3]]
  itemfit(mod, X2=TRUE)
  
  library(rpf)
  spec <- list()
  spec[1:5] <- rpf.drm()
  param <- sapply(coef(mod)[-6], function(x) x)
  param[3:4,] <- rpf::logit(param[3:4,])
  dat2 <- as.data.frame(lapply(as.data.frame(dat), ordered, levels=0:1))
  grp <- list(spec=spec, param=param, mean=0, cov=diag(1), data=dat2)
  got <- SitemFit(grp, method="pearson", alt=FALSE)
}
