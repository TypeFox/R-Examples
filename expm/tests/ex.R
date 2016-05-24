library(expm)

source(system.file("test-tools.R", package= "expm"), keep.source=FALSE)

## Note that these results are achieved with the default
## settings order=8, method="Pade" -- accuracy could
## presumably be improved still further by some tuning
## of these settings.


## ----------------------------
## Test case 1 from Ward (1977)
## ----------------------------
T1 <- rbind(c(4, 2, 0),
            c(1, 4, 1),
            c(1, 1, 4))
(m1   <- expm(T1, method="Pade"))
(m1O  <- expm(T1, method="PadeO"))# very slightly different
(m1T  <- expm(T1, method="Taylor"))
(m1TO <- expm(T1, method="TaylorO"))
## True Eigenvalue Decomposition of T1
s2 <- sqrt(2)
eV1 <- matrix(c(s2,s2,s2,  -2,1,1,  2,-1,-1) / sqrt(6),
              3,3)
L1 <- diag(lm1 <- c(6, 3, 3))
stopifnot(
    all.equal(eV1 %*% L1, T1 %*% eV1, tolerance=1e-15)
    )
## However,  eV1 is not orthogonal, but of rank 2

if(FALSE) { ## require("Rmpfr")) { ## 200 bit precision version of that
  S2 <- sqrt(mpfr(2,200))
  E1 <- c(S2,S2,S2,  -2,1,1,  2,-1,-1) / sqrt(mpfr(6,200))
  dim(E1) <- c(3,3)
  print(E1 %*% L1)
  print(E1)
}

## "true" result
m1.t <- matrix(c(147.866622446369, 127.781085523181, 127.781085523182,
                 183.765138646367, 183.765138646366, 163.679601723179,
                 71.797032399996,  91.8825693231832, 111.968106246371), 3,3)
stopifnot(all.equal(m1.t, m1,  check.attributes=FALSE, tolerance = 1e-13),
          all.equal(m1.t, m1O, check.attributes=FALSE, tolerance = 1e-13),
          all.equal(m1.t,m1T,  check.attributes=FALSE, tolerance = 1e-13),
          all.equal(m1.t,m1TO, check.attributes=FALSE, tolerance = 1e-13),
          all.equal(m1.t, expm(T1,"Ward77"),    tolerance = 1e-13),
          all.equal(m1.t, expm(T1,"R_Pade"),    tolerance = 1e-13),
          all.equal(m1.t, expm(T1,"R_Ward77"),  tolerance = 1e-13))
## -- these agree with ward (1977, p608)
##
m1.2 <- try( expm(T1, "R_Eigen") ) ## 32-bit: gives an error from solve; 64-bit "ok"
if(!inherits(m1.2, "try-error")) {
    if(FALSE)## with libatlas R_Eigen is "sehr eigen"
    stopifnot(all.equal(m1.t, m1.2, check.attributes=FALSE))
    ## but it's less accurate:
    print(all.equal(m1.t, m1.2, check.attributes=FALSE, tolerance= 1e-12))
    ##-> rel.diff = 6.44e-10 / 6.2023e-10
    ##__ BUT  0.1228099
    ##__ with libatlas (ubuntu 12.04 libatlas-base-dev Version: 3.8.4-3build1)
}

##
## ----------------------------
## Test case 2 from Ward (1977)
## ----------------------------
T2 <- t(matrix(c(
                    29.87942128909879, .7815750847907159, -2.289519314033932,
                    .7815750847907159, 25.72656945571064,  8.680737820540137,
                    -2.289519314033932, 8.680737820540137,  34.39400925519054),
                  3, 3))
(m2 <- expm(T2, method="Pade"))
##                   [,1]               [,2]               [,3]
##[1,]   5496313853692357 -18231880972009844 -30475770808580828
##[2,] -18231880972009852  60605228702227024 101291842930256144
##[3,] -30475770808580840 101291842930256144 169294411240859072
## -- which agrees with Ward (1977) to 13 significant figures
(m2O  <- expm(T2, method="PadeO"))
(m2T  <- expm(T2,method="Taylor"))
(m2TO <- expm(T2,method="TaylorO"))

m2.t <- matrix(c(5496313853692216, -18231880972008932, -30475770808579672,
                 -18231880972008928, 60605228702222480, 101291842930249776,
                 -30475770808579672, 101291842930249808, 169294411240850528),
               3, 3)

## -- in this case a very similar degree of accuracy -- even Taylor is ok
stopifnot(all.equal(m2.t, m2, check.attributes=FALSE, tolerance = 1e-12),
          all.equal(m2.t, m2O,check.attributes=FALSE, tolerance = 1e-12),
          all.equal(m2.t,m2T, check.attributes=FALSE, tolerance = 1e-12),
          all.equal(m2.t,m2TO,check.attributes=FALSE, tolerance = 1e-12),
          all.equal(m2.t, expm(T2,"Ward77"),   tolerance = 1e-12),
          all.equal(m2.t, expm(T2,"R_Ward77"), tolerance = 1e-12),
          all.equal(m2.t, expm(T2,"R_Pade"),   tolerance = 1e-12),
          TRUE)

## ----------------------------
## Test case 3 from Ward (1977)
## ----------------------------
T3 <- t(matrix(c(
    -131, 19, 18,
    -390, 56, 54,
    -387, 57, 52), 3, 3))
(m3 <- expm(T3, method="Pade"))
##                    [,1]                [,2]                [,3]
##[1,] -1.5096441587713636 0.36787943910439874 0.13533528117301735
##[2,] -5.6325707997970271 1.47151775847745725 0.40600584351567010
##[3,] -4.9349383260294299 1.10363831731417195 0.54134112675653534
## -- agrees to 10dp with Ward (1977), p608.
(m3O  <- expm(T3, method="PadeO"))
(m3T  <- expm(T3,method="Taylor"))
(m3TO <- expm(T3,method="TaylorO"))

m3.t <- matrix(c(-1.50964415879218, -5.6325707998812, -4.934938326092,
                 0.367879439109187, 1.47151775849686, 1.10363831732856,
                 0.135335281175235, 0.406005843524598, 0.541341126763207),
               3,3)

stopifnot(all.equal(m3.t, m3,           check.attributes=FALSE, tolerance = 3e-11),
					#			  ^^^^^
					# 1.2455e-11 for libatlas (above)
          all.equal(m3.t, m3T,          check.attributes=FALSE, tolerance = 1e-11),
          all.equal(m3.t, m3O,          check.attributes=FALSE, tolerance = 1e-11),
          all.equal(m3.t, m3TO,         check.attributes=FALSE, tolerance = 1e-11),
          all.equal(m3.t, expm(T3,"R_Eigen"), tolerance = 1e-11),
          all.equal(m3.t, expm(T3,"Ward77"), tolerance = 1e-11),
          all.equal(m3.t, expm(T3,"R_Ward"), tolerance = 1e-11),
          all.equal(m3.t, expm(T3,"R_Pade"), tolerance = 1e-11),
          TRUE)
## -- in this case, a similar level of agreement with Ward (1977).

## ----------------------------
## Test case 4 from Ward (1977)
## ----------------------------
T4 <-
    array(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-10,
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
          dim = c(10, 10))
(m4   <- expm(T4, method="Pade"))
(m4O  <- expm(T4, method="PadeO"))
(m4T  <- expm(T4,method="Taylor"))
(m4TO <- expm(T4,method="TaylorO"))

stopifnot(all.equal(m4  [,10], 1/gamma(10:1), tolerance=1e-14),
          all.equal(m4O [,10], 1/gamma(10:1), tolerance=1e-14),
          all.equal(m4T [,10], 1/gamma(10:1), tolerance=1e-14),
          all.equal(m4TO[,10], 1/gamma(10:1), tolerance=1e-14),
	  all.equal(m4, m4O, check.attributes=FALSE, tolerance=5e-15),
	  all.equal(m4, m4T, check.attributes=FALSE, tolerance=5e-15),
	  all.equal(m4, m4TO,check.attributes=FALSE, tolerance=5e-15),
          all.equal(m4, expm(T4,"Ward77"), check.attributes=FALSE, tolerance = 1e-14),
          all.equal(m4, expm(T4,"R_Ward"), check.attributes=FALSE, tolerance = 1e-14),
          all.equal(m4, expm(T4,"R_Pade"), check.attributes=FALSE, tolerance = 1e-14),
          max(abs(m4 - expm(T4,"R_Eigen"))) < 1e-7)
## here expm(., EV ) is accurate only to 7 d.p., whereas
##      expm(.,Pade) is correct to at least 14 d.p.

### Test case with diagonalizable matrix with multiple Eigen values:
A4 <- rbind(c(-1, 3, -1),
            c(-3, 5, -1),
            c(-3, 3,  1))
Ea4 <- eigen(A4)
stopifnot(all.equal(Ea4$values, (lam4 <- c(2,2,1))))
## However, the eigen values don't show the nice property:
V4 <- Ea4$vectors
crossprod(V4)
## i.e., they are *not* orthogonal
## but still diagonalize:
stopifnot(all.equal(A4,
                    V4 %*% diag(x=lam4) %*% solve(V4)))
## whereas this diagonalizes *and* looks nice
W4 <- rbind(c(1, 1, -1),
	    c(1, 1,  0),
	    c(1, 0,  3))
(sW4 <- solve(W4))
stopifnot(all.equal(diag(x = c(1,2,2)),
		    solve(W4) %*% A4 %*% W4 ))
