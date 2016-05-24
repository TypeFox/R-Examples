test.distsamp.covs <- function() {
    y <- matrix(rep(4:1, 10), 5, 2, byrow=TRUE)
    siteCovs <- data.frame(x = c(0, 2, 3, 4, 1))
    umf <- unmarkedFrameDS(y = y, siteCovs = siteCovs,
        dist.breaks=c(0, 5, 10)/1000, survey="line", tlength=rep(1, 5),
        unitsIn="km")
    fm <- distsamp(~ x ~ x, data = umf)

    lam <- fm['state']
    det <- fm['det']

    checkEqualsNumeric(coef(lam), c(1.4340999, -0.1102387), tolerance = 1e-4)
    checkEqualsNumeric(coef(det), c(-4.64686395, -0.09337832), tolerance = 1e-4)

    lam.lc <- linearComb(fm, type = 'state', c(1, 2))
    det.lc <- linearComb(fm, type = 'det', c(1, 2))

    checkEqualsNumeric(coef(lam.lc), 1.213623, tol = 1e-4)
    checkEqualsNumeric(coef(det.lc), -4.833621, tol = 1e-4)

    checkEqualsNumeric(coef(backTransform(lam.lc)), 3.365655, tol = 1e-4)
    checkEqualsNumeric(coef(backTransform(det.lc)), 0.007957658, tol = 1e-4)
    }



test.distsamp.line.keyfuns <- function()
{
    y <- structure(c(7, 7, 12, 9, 9, 11, 9, 5, 7, 6, 25, 26, 30, 26, 23,
        24, 20, 33, 26, 32, 5, 3, 8, 7, 1, 4, 4, 7, 7, 6, 3, 1, 1, 4,
        4, 4, 3, 6, 2, 3), .Dim = c(10L, 4L))
    umf <- unmarkedFrameDS(y = y, dist.breaks=c(0, 3, 15, 18, 20),
        survey="line", unitsIn="m", tlength=rep(100, nrow(y)))

    fm.halfnorm <- distsamp(~1~1, umf)
    D <- backTransform(fm.halfnorm, type="state")
    S <- backTransform(fm.halfnorm, type="det")
    checkEqualsNumeric(coef(D), 129.5509, tol=1e-4)
    checkEqualsNumeric(SE(D), 9.446125, tol=1e-4)
    checkEqualsNumeric(coef(S), 18.15386, tol=1e-4)
    checkEqualsNumeric(SE(S), 2.893362, tol=1e-4)

    fm.exp <- distsamp(~1~1, umf, keyfun="exp", starts=c(4, 0))
    D <- backTransform(fm.exp, type="state")
    S <- backTransform(fm.exp, type="det")
    checkEqualsNumeric(coef(D), 144.8802, tol=1e-4)
    checkEqualsNumeric(SE(D), 14.31655, tol=1e-4)
    checkEqualsNumeric(coef(S), 31.75738, tol=1e-4)
    checkEqualsNumeric(SE(S), 9.711254, tol=1e-4)

    fm.haz <- distsamp(~1~1, umf, keyfun="hazard", starts=c(4, 3, 1))
    D <- backTransform(fm.haz, type="state")
    Sh <- backTransform(fm.haz, type="det")
    Sc <- backTransform(fm.haz, type="scale")
    checkEqualsNumeric(coef(D), 137.0375, tol=1e-4)
    checkEqualsNumeric(SE(D), 16.82505, tol=1e-4)
    checkEqualsNumeric(coef(Sh), 15.90262, tol=1e-4)
    checkEqualsNumeric(SE(Sh), 5.099981, tol=1e-4)
    checkEqualsNumeric(coef(Sc), 0.8315524, tol=1e-4)
    checkEqualsNumeric(SE(Sc), 0.4753275, tol=1e-4)

    fm.unif <- distsamp(~1~1, umf, keyfun="uniform")
    D <- backTransform(fm.unif, type="state")
    checkEqualsNumeric(coef(D), 107.5000, tol=1e-4)

    checkEqualsNumeric(coef(fm.halfnorm),
                       coef(update(fm.halfnorm, engine="R")))
    checkEqualsNumeric(coef(fm.exp),
                       coef(update(fm.exp, engine="R")))
    checkEqualsNumeric(coef(fm.halfnorm),
                       coef(update(fm.halfnorm, engine="R")))
    checkEqualsNumeric(coef(fm.halfnorm),
                       coef(update(fm.halfnorm, engine="R")))

}



test.distsamp.point.keyfuns <- function()
{
    y <- structure(c(1, 0, 0, 0, 0, 0, 3, 1, 1, 0, 16, 15, 18, 14, 22,
        24, 12, 20, 20, 21, 10, 9, 9, 5, 6, 6, 6, 9, 5, 6, 6, 6, 4, 2,
        6, 3, 3, 3, 1, 4), .Dim = c(10L, 4L))

    umf <- unmarkedFrameDS(y = y, dist.breaks=c(0, 3, 15, 18, 20),
        survey="point", unitsIn="m", tlength=rep(100, 20))

    fm.halfnorm <- distsamp(~1~1, umf)
    D <- backTransform(fm.halfnorm, type="state")
    S <- backTransform(fm.halfnorm, type="det")
    checkEqualsNumeric(coef(D), 316.1711, tol=1e-4)
    checkEqualsNumeric(SE(D), 37.08797, tol=1e-4)
    checkEqualsNumeric(coef(S), 18.05958, tol=1e-4)
    checkEqualsNumeric(SE(S), 3.341798, tol=1e-4)

    fm.exp <- distsamp(~1~1, umf, keyfun="exp", starts=c(6, 0))
    D <- backTransform(fm.exp, type="state")
    S <- backTransform(fm.exp, type="det")
    checkEqualsNumeric(coef(D), 369.7526, tol=1e-4)
    checkEqualsNumeric(SE(D), 68.11901, tol=1e-4)
    checkEqualsNumeric(coef(S), 28.90848, tol=1e-4)
    checkEqualsNumeric(SE(S), 11.66219, tol=1e-4)

    fm.haz <- distsamp(~1~1, umf, keyfun="hazard", starts=c(5, 3, 1))
    D <- backTransform(fm.haz, type="state")
    Sh <- backTransform(fm.haz, type="det")
    Sc <- backTransform(fm.haz, type="scale")
    checkEqualsNumeric(coef(D), 266.3911, tol=1e-4)
    checkEqualsNumeric(SE(D), 20.45144, tol=1e-4)
    checkEqualsNumeric(coef(Sh), 18.69351, tol=1e-4)
    checkEqualsNumeric(SE(Sh), 0.8950444, tol=1e-4)
    checkEqualsNumeric(coef(Sc), 5.797366, tol=1e-4)
    checkEqualsNumeric(SE(Sc), 4.054381, tol=1e-4)

    fm.unif <- distsamp(~1~1, umf, keyfun="uniform")
    D <- backTransform(fm.unif, type="state")
    checkEqualsNumeric(coef(D), 236.3451, tol=1e-4)

    checkEqualsNumeric(coef(fm.halfnorm),
                       coef(update(fm.halfnorm, engine="R")))
    checkEqualsNumeric(coef(fm.exp),
                       coef(update(fm.exp, engine="R")))
    checkEqualsNumeric(coef(fm.halfnorm),
                       coef(update(fm.halfnorm, engine="R")))
    checkEqualsNumeric(coef(fm.halfnorm),
                       coef(update(fm.halfnorm, engine="R")))

}


