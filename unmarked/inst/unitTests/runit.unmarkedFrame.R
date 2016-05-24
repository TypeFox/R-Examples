test.emptyframe <- function() {
  checkException(umf <- unmarkedFrame())
}

test.frame <- function() {
  M <- 10
  J <- 3
  y <- matrix(rbinom(J * M, 1, 0.5), M, J)
  siteCovs <- data.frame(a = rnorm(M), b = factor(gl(2,5)))
  umf <- unmarkedFrame(y = y, siteCovs = siteCovs)
}


test.umfDS.args <- function() {
    y <- matrix(1, 1, 2)
    s <- "point"
    d <- 0:2
    uin <- "m"
    sc <- data.frame(1)
    oc <- matrix(1, 2)
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, obsCovs=oc,
        survey=s, dist.breaks=d, unitsIn=uin))
    umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s, dist.breaks=d,
        unitsIn=uin)
    checkException(obsCovs(umf) <- oc)
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s,
        dist.breaks=d))
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, dist.breaks=d,
        unitsIn=uin))
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s,
        unitsIn=uin))
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s,
        dist.breaks=0:3, unitsIn=uin))
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s,
        dist.breaks=1:3, unitsIn=uin))

    }



test.obsToY <- function() {
    y <- matrix(c(
        1, 0, 0,
        2, 1, 0,
        1, 0, 1,
        2, 1, 2,
        1, 0, 3,
        1, 1, 1), nrow=6, ncol=3, byrow=TRUE)
    oc <- matrix(c(
        1, 0,
        2, 1,
        1, 1,
        NA, 0,
        1, NA,
        NA, NA), nrow=6, ncol=2, byrow=TRUE)
    umf <- unmarkedFrameMPois(y = y, obsCovs = list(x=oc), type="double")
    o2y <- obsToY(umf)

    checkEquals(o2y, matrix(1, 2, 3))
    oc.na <- is.na(oc)
    oc.na %*% o2y

    }
