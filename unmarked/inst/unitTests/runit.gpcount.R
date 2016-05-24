

test.gpcount.fit <- function()
{
    y <- matrix(c(0,0,0, 1,0,1, 2,2,2,
                3,2,3, 2,2,2, 1,1,1,
                NA,0,0, 0,0,0, 0,0,0,
                3,3,3, 3,1,3, 2,2,1,
                0,0,0, 0,0,0, 0,0,0), 5, 9, byrow=TRUE)
    siteCovs <- data.frame(x = c(0,2,-1,4,-1))
    obsCovs <- list(o1 = matrix(seq(-3, 3, length=length(y)), 5, 9))
    obsCovs$o1[5,4:6] <- NA
    yrSiteCovs <- list(yr=matrix(c('1','2','2'), 5, 3, byrow=TRUE))
    yrSiteCovs$yr[4,2] <- NA

    umf <- unmarkedFrameGPC(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, numPrimary=3)
    fm <- gpcount(~x, ~yr, ~o1, data = umf, K=23)
    checkEquals(fm@sitesRemoved, integer(0))
    checkEqualsNumeric(coef(fm),
        c(1.14754541, 0.44499137, -1.52079283, -0.08881542,  2.52037155, -0.10950615), tol = 1e-5)

    fm0c <- gpcount(~1, ~1, ~1, umf, K=20)
    fm0r <- gpcount(~1, ~1, ~1, umf, K=20, engine="R")
    checkEqualsNumeric(coef(fm0c), coef(fm0r))

}
