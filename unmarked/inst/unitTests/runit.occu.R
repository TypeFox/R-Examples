test.occu.fit.simple.1 <- function() {
  
  y <- matrix(rep(1,10),5,2)
  umf <- unmarkedFrameOccu(y = y)
  fm <- occu(~ 1 ~ 1, data = umf)

  occ <- fm['state']
  det <- fm['det']

  occ <- coef(backTransform(occ))
  checkEqualsNumeric(occ,1)

  det <- coef(backTransform(det))
  checkEqualsNumeric(det,1)

  bt <- backTransform(fm, type = 'state')
  checkEqualsNumeric(coef(bt), 1)

  bt <- backTransform(fm, type = 'det')
  checkEqualsNumeric(coef(bt), 1)

}

test.occu.fit.simple.0 <- function() {

  y <- matrix(rep(0,10),5,2)
  umf <- unmarkedFrameOccu(y = y)
  fm <- occu(~ 1 ~ 1, data = umf)

  occ <- fm['state']
  det <- fm['det']

  occ <- coef(backTransform(occ))
  checkEqualsNumeric(occ, 0, tolerance = 1e-4)

  det <- coef(backTransform(det))
  checkEqualsNumeric(det,0, tolerance = 1e-4)

  bt <- backTransform(fm, type = 'state')
  checkEqualsNumeric(coef(bt), 0, tolerance = 1e-4)

  bt <- backTransform(fm, type = 'det')
  checkEqualsNumeric(coef(bt), 0, tolerance = 1e-4)


}

test.occu.fit.covs <- function() {

  y <- matrix(rep(0:1,10),5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ x, data = umf)
  
  occ <- fm['state']
  det <- fm['det']

  checkException(occ <- coef(backTransform(occ)))

  checkEqualsNumeric(coef(occ), c(8.590737, 2.472220), tolerance = 1e-4)
  checkEqualsNumeric(coef(det), c(0.44457, -0.14706, 0.44103), tolerance = 1e-4)

  occ.lc <- linearComb(fm, type = 'state', c(1, 0.5))
  det.lc <- linearComb(fm, type = 'det', c(1, 0.3, -0.3))
    
  checkEqualsNumeric(coef(occ.lc), 9.826848, tol = 1e-4)
  checkEqualsNumeric(coef(det.lc), 0.2681477, tol = 1e-4)

  checkEqualsNumeric(coef(backTransform(occ.lc)), 1, tol = 1e-4)
  checkEqualsNumeric(coef(backTransform(det.lc)), 0.5666381, tol = 1e-4)

  checkException(backTransform(fm, type = "state"))
  checkException(backTransform(fm, type = "det"))

  fitted <- fitted(fm)
  checkEqualsNumeric(fitted, structure(c(0.5738, 0.5014, 0.4318, 0.38581, 0.50171, 0.53764, 
0.46563, 0.40283, 0.39986, 0.79928), .Dim = c(5L, 2L)), tol = 1e-5)

}

test.occu.fit.covs.0 <- function() {

  y <- matrix(rep(0,10),5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  checkException(fm <- occu(~ o1 + o2 ~ x, data = umf))

}

test.occu.fit.NA <- function() {

  y <- matrix(rep(0:1,10),5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  siteCovs[3,1] <- NA
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ x, data = umf)
  checkEquals(fm@sitesRemoved, 3)
  checkEqualsNumeric(coef(fm), c(8.70123, 4.58255, 0.66243, -0.22862, 0.58192), tol = 1e-5)

  obsCovs[10,2] <- NA
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ x, data = umf)
  checkEquals(fm@sitesRemoved, 3)
  checkEqualsNumeric(coef(fm), c(8.91289, 1.89291, -1.42471, 0.67011, -8.44608), tol = 1e-5)

}

## Add some checks here.
test.occu.offest <- function() {

  y <- matrix(rep(0:1,10),5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ offset(x), data = umf)
  checkEqualsNumeric(coef(fm),
                     structure(c(9.74361, 0.44327, -0.14683, 0.44085), .Names = c("psi(Int)", 
"p(Int)", "p(o1)", "p(o2)")), tol = 1e-5)
  fm <- occu(~ o1 + offset(o2) ~ offset(x), data = umf)
  checkEqualsNumeric(coef(fm), structure(c(8.59459, 0.97574, -0.3096), .Names = c("psi(Int)", 
"p(Int)", "p(o1)")), tol=1e-5)

}
