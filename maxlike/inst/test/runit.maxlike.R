test.maxlike.fit1 <- function() {

    data(MaungaWhau)
    elev <- raster(MaungaWhau$elev, 0, 61, 0, 87)
    precip <- raster(MaungaWhau$precip, 0, 61, 0, 87)
    xy <- MaungaWhau$xy

    # Stack them and make sure they are named
    ep <- stack(elev, precip)
    names(ep) <- c("elev", "precip")

    # Fit a model
    fm <- maxlike(~elev + I(elev^2) + precip, ep, xy)

    # Check estimates
    checkEqualsNumeric(coef(fm),
                       c(0.5366934, 2.4465578, -2.3575862, 2.1310296),
                       tol=1e-6)

    # Check variance-covariance matrix
    checkEqualsNumeric(vcov(fm), matrix(c(
                       0.05204765,  0.03300724, -0.03589617,  0.03561091,
                       0.03300724,  0.03921600, -0.03373490,  0.03307982,
                       -0.03589617, -0.03373490, 0.03618861, -0.03398646,
                       0.03561091,  0.03307982, -0.03398646,  0.04569138),
                       4, 4, byrow=TRUE), tol=1e-6)

    # Add missing values and refit
    elev2 <- elev
    elev2[c(1,5)] <- NA
    xy2 <- xy
    xy2[2,] <- NA
    ep2 <- stack(elev2, precip)
    names(ep2) <- c("elev", "precip")
    fm2 <- maxlike(~elev + I(elev^2) + precip, ep2, xy2)
    checkEqualsNumeric(fm2$pix.removed, c(1,5))
    checkEqualsNumeric(fm2$pts.removed, 2)

    checkEqualsNumeric(AIC(fm2), 16017.43, tol=1e-6)
    checkEqualsNumeric(logLik(fm2), -8004.717, tol=1e-6)
    checkEqualsNumeric(confint(fm2),
                       matrix(c(0.09770203,  0.992682,
                                2.06111183,  2.839246,
                                -2.73433869, -1.987336,
                                1.71726661,  2.558565),
                       4, 2, byrow=TRUE), tol=1e-6)

    # Test update
    fm2u <- update(fm2)
    checkEquals(fm2, fm2u)

    # refit using fixed values
    fix <- c(0,2,-2,2)
    checkException(fm3 <- update(fm2, fixed=fix))
    fix <- c(NA,NA,NA,NA)
    checkException(fm3 <- update(fm2, fixed=fix))
    fix <- c(0,NA,NA,NA)
    fm3 <- update(fm2, fixed=fix)
    checkEqualsNumeric(coef(fm3),
                       c(0.000000, 2.121397, -2.003473, 1.781255),
                       tol=1e-6)

    # Test cloglog
    fm4 <- update(fm2, link="cloglog", savedata=TRUE)
    checkEqualsNumeric(coef(fm4),
                       c(0.1765974, 2.0890600, -2.0136279,  1.7778590),
                       tol=1e-6)

    # Test predict
    fm5 <- update(fm2, savedata=TRUE)
    psi.hat <- predict(fm5)
    checkEqualsNumeric(cellStats(psi.hat, "mean"), 0.3538176, tol=1e-6)

    psi.hat <- predict(fm4)
    checkEqualsNumeric(cellStats(psi.hat, "mean"), 0.3817011, tol=1e-6)

}






