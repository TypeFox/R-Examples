#####-----------------------------------------------------------------------
## CEP estimate based on correlated bivariate normal distribution - mvnEll.R, hoyt.R
CEPCorrNormal <-
function(CEPlevel=0.5, ctr=c(0, 0), sigma=diag(length(ctr)), accuracy=FALSE) {
    p <- ncol(sigma)
    CorrNormal <- if(accuracy) {  # POA != POI
        ## quantile from offset circle probability - most general case
        qmvnEll(CEPlevel, mu=numeric(p), sigma=sigma, e=diag(p), x0=ctr)
    } else {                      # POA == POI
        if(p == 2L) {             # 2D -> exact Hoyt distribution
            HP <- getHoytParam(sigma)
            qHoyt(CEPlevel, qpar=HP$q, omega=HP$omega)
        } else {                  # 1D/3D case
            qmvnEll(CEPlevel, mu=numeric(p), sigma=sigma, e=diag(p), x0=numeric(p))
        }
    }

    setNames(CorrNormal, CEPlevel)
}

####-----------------------------------------------------------------------
## Grubbs estimates (Grubbs, 1964, p54, p55-56) - grubbs.R
## Grubbs-Patnaik CEP estimate based on Patnaik two-moment central chi^2 approximation
CEPGrubbsPatnaik <-
function(CEPlevel=0.5, ctr=c(0, 0), sigma=diag(length(ctr)), accuracy=FALSE) {
    GPP <- getGrubbsParam(sigma, ctr=ctr, accuracy=accuracy)
    GrubbsPatnaik <- qChisqGrubbs(CEPlevel, m=GPP$m, v=GPP$v, n=GPP$n, type="Patnaik")
    setNames(GrubbsPatnaik, CEPlevel)
}

## Grubbs-Pearson CEP estimate based on Pearson three-moment central chi^2 approximation
CEPGrubbsPearson <-
function(CEPlevel=0.5, ctr=c(0, 0), sigma=diag(length(ctr)), accuracy=FALSE) {
    GPP <- getGrubbsParam(sigma, ctr=ctr, accuracy=accuracy)
    GrubbsPearson <- qChisqGrubbs(CEPlevel, m=GPP$m, v=GPP$v, nPrime=GPP$nPrime, type="Pearson")
    setNames(GrubbsPearson, CEPlevel)
}

## Grubbs-Liu CEP estimate based on four-moment non-central chi^2
## approximation (Liu, Tang & Zhang, 2009)
CEPGrubbsLiu <-
function(CEPlevel=0.5, ctr=c(0, 0), sigma=diag(length(ctr)), accuracy=FALSE) {
    GPP <- getGrubbsParam(sigma, ctr=ctr, accuracy=accuracy)
    GrubbsLiu <- qChisqGrubbs(CEPlevel, m=GPP$m, v=GPP$v, muX=GPP$muX,
                              varX=GPP$varX, l=GPP$l, delta=GPP$delta, type="Liu")
    setNames(GrubbsLiu, CEPlevel)
}

#####-----------------------------------------------------------------------
## Rayleigh CEP estimate from Singh, 1992 - rayleigh.R, rice.R, maxwell.R
CEPRayleigh <-
function(CEPlevel=0.5, ctr=c(0, 0), sigma=diag(length(ctr)), accuracy=FALSE, doRob=FALSE, xy) {
    p <- ncol(sigma)
    Rayleigh <- if(p == 1L) {     # 1D
        if(accuracy) {            # POA != POI -> non-central F distribution
            RiceParam <- getRiceParam(xy, doRob=doRob)
            NU  <- RiceParam$nu
            SD  <- setNames(RiceParam$sigma["sigma"], NULL)
            N   <- length(xy)
            NCP <- 1*(NU/SD)^2    # N*f^2, N=1 because interested in individual shot, not mean
            SD * sqrt(qf(CEPlevel, df1=1, df2=length(xy)-1, ncp=NCP))
        } else {                  # POA = POI -> half normal / central F distribution
            SD <- getRayParam(xy, doRob=doRob)$sigma["sigma"]
            SD * sqrt(qchisq(CEPlevel, df=p))
            # SD * sqrt(qf(CEPlevel, df1=1, df2=length(xy)-1))
        }
    } else if(p == 2L) {          # 2D
        if(accuracy) {            # POA != POI -> Rice distribution
            RiceParam <- getRiceParam(xy, doRob=doRob)
            qRice(CEPlevel, nu=RiceParam$nu, sigma=RiceParam$sigma["sigma"])
        } else {                  # POA = POI -> Rayleigh
            RayParam <- getRayParam(xy, doRob=doRob)
            qRayleigh(CEPlevel, scale=RayParam$sigma["sigma"])
        }
    } else if(p == 3L) {          # 3D
        MaxParam <- getRayParam(xy, doRob=doRob)
        if(accuracy) {            # POA != POI -> offset sphere probability
            ## circular covariance matrix with estimated M-B param sigma
            sigMat <- diag(rep(MaxParam$sigma["sigma"]^2, p))
            qmvnEll(CEPlevel, sigma=sigMat, mu=numeric(p), x0=ctr, e=diag(p))
        } else {                  # POA = POI -> Maxwell-Boltzmann distribution
            qMaxwell(CEPlevel, sigma=MaxParam$sigma["sigma"])
        }
    } else {
        warning("Rayleigh CEP with accuracy=TRUE is only available for 2D/3D-data")
        rep(NA_real_, length(CEPlevel))
    }

    setNames(Rayleigh, CEPlevel)
}

#####-----------------------------------------------------------------------
## Krempasky CEP estimate from Krempasky, 2003
CEPKrempasky <-
function(CEPlevel=0.5, sigma=diag(2), accuracy=FALSE) {
    p <- ncol(sigma)
    ## only available for 2D 50% POA=POI case
    Krempasky50 <- if((p == 2L) && !accuracy) {
        ## estimated correlation, covariance and standard deviations
        rho    <- cov2cor(sigma)[1, 2]
        covXY  <- sigma[1, 2]
        sigmaX <- sqrt(sigma[1, 1])
        sigmaY <- sqrt(sigma[2, 2])

        ## rotation angle gamma and rotation matrix
        gamma <- atan((-2*covXY + sqrt(4*covXY^2 + (sigmaY^2 - sigmaX^2)^2)) /
                          (sigmaY^2 - sigmaX^2))

        A <- cbind(c(cos(gamma), sin(gamma)), c(-sin(gamma), cos(gamma)))

        ## covariance matrix, correlation and standard deviations of rotated data
        sigmaRot   <- t(A) %*% sigma %*% A     # covariance matrix
        rhoPrime   <- cov2cor(sigmaRot)[1, 2]  # correlation
        sigmaPrime <- sqrt(sigmaRot[1, 1])
        # isTRUE(all.equal(sigmaRot[1, 1], sigmaRot[2, 2])) # supposed to be equal

        CEP00 <- sigmaPrime*1.1774100225154746910115693264596996377473856893858
        C2    <- 0.3267132048600136726456919696354558579811249664099
        C4    <- 0.0568534980324428522403361717417556935145611584613
        ## CEP00 <- sigmaPrime*sqrt(2*log(2)) # sigmaPrime * qRayleigh(0.5, 1)
        ## C2    <- 0.5*(1 - log(2)/2)
        ## C4    <- C2*(log(2) - log(2)^2/4 - 0.5) + 3/8 - (9/16)*log(2) + (3/16)*log(2)^2 - (1/64)*log(2)^3 - C2^2*log(2)/2
        CEP00*(1 - 0.5*C2*rhoPrime^2 - 0.5*(C4 + 0.25*C2^2)*rhoPrime^4)
    } else {
        if(p != 2L) {
            warning("Krempasky CEP is only available for 2D-data")
        }

        if(accuracy) {
            warning("Krempasky CEP is only available for accuracy=FALSE")
        }
        rep(NA_real_, length(CEPlevel))
    }

    if(any(CEPlevel != 0.5)) {
        warning("Krempasky CEP is only available for CEPlevel 0.5")
    }

    Krempasky <- ifelse(CEPlevel == 0.5, Krempasky50, NA_real_)
    setNames(Krempasky, CEPlevel)
}

#####-----------------------------------------------------------------------
## Ignani estimate from van Ignani (2010)
CEPIgnani <-
function(CEPlevel=0.5, sigma=diag(2), accuracy=FALSE) {
    p  <- ncol(sigma)

    ## make sure eigenvalues >= 0 when very small
    ev      <- eigen(sigma)$values    # eigenvalues
    lambda  <- ev*sign(ev)
    alpha   <- sqrt(lambda[2] / lambda[1])
    beta    <- if(p == 3L) { sqrt(lambda[3] / lambda[2]) } else { 0 }
    betaVec <- c(1, beta, beta^2, beta^3)

    ## coefficients for polynomials in table 1
    con <- textConnection("coef   R0.5    R0.9   R0.95   R0.99
                           c11  0.6754  1.6494  1.9626  2.5686
                           c12 -0.1547  0.0332 -0.0906 -0.1150
                           c13  0.2616  1.3376  1.3214 -0.3475
                           c14  1.0489 -0.8445 -0.3994  1.3570
                           c21 -0.0208 -0.0588  0.0100  0.1479
                           c22  1.1739 -0.5605  0.2722  0.9950
                           c23  1.9540 -4.7170 -5.4821  1.3223
                           c24 -5.5678  5.7135  3.9732 -4.8917
                           c31  1.1009  0.3996  0.0700 -0.4285
                           c32 -2.6375  1.5739  0.0462 -1.9795
                           c33 -1.4838  5.3623  7.1658 -1.1104
                           c34  6.5837 -7.9347 -5.7194  6.9617
                           c41 -0.5821  0.1636  0.4092  0.7371
                           c42  1.5856 -1.0747 -0.1953  1.2495
                           c43 -0.0678 -1.7785 -3.0134 -0.2061
                           c44 -2.3324  3.2388  2.4661 -2.8968
                          ")
    coefDF <- read.table(con, header=TRUE)
    close(con)

    getR <- function(level) {
        column <- paste0("R", level)
        c1 <- coefDF[ 1:4,  column]
        c2 <- coefDF[ 5:8,  column]
        c3 <- coefDF[ 9:12, column]
        c4 <- coefDF[13:16, column]
        sqrt(lambda[1]) * (crossprod(c1, betaVec) +
                           crossprod(c2, betaVec)*alpha +
                           crossprod(c3, betaVec)*alpha^2 +
                           crossprod(c4, betaVec)*alpha^3)
    }

    Ignani <- if((p %in% c(2L, 3L)) && !accuracy) {
        vapply(CEPlevel, function(x) {
            if(x %in% c(0.5, 0.9, 0.95, 0.99)) { getR(x) } else { NA_real_ } }, numeric(1))
    } else {
        if(!(p %in% c(2L, 3L))) {
            warning("Ignani CEP is only available for 2D/3D-data")
        }

        if(accuracy) {
            warning("Ignani CEP is only available for accuracy=FALSE")
        }

        rep(NA_real_, length(CEPlevel))
    }

    if(!any(CEPlevel %in% c(0.5, 0.9, 0.95, 0.99))) {
        warning("Ignani CEP is only available for CEPlevel 0.5, 0.9, 0.95, 0.99")
    }

    setNames(Ignani, CEPlevel)
}

#####-----------------------------------------------------------------------
## RMSE-based estimate from van Diggelen (2007)
## CEP = sqrt(qchisq(0.5, 2)) * RMSEx
## RMSExy = sqrt(2) * RMSEx = sqrt(2) * (1/sqrt(qchisq(0.5, 2)))*CEP
CEPRMSE <-
function(CEPlevel=0.5, sigma=diag(2), accuracy=FALSE, xy) {
    p <- ncol(sigma)

    ## root mean squared error
    RMSE_AT <- sqrt(sum(colMeans(xy^2))) # non-centered data
    RMSE_AF <- sqrt(sum(diag(sigma)))    # centered data

    ## conversion factor from RMSExy or RMSExyz to CEP
    ## 2D: (1/sqrt(2)) * qRayleigh(CEPlevel, scale=1)
    ##     for 50%, this is just sqrt(log(2))
    ## 3D: (1/sqrt(3)) * qMaxwell(CEPlevel, sigma=1)
    RMSE2CEP <- (1/sqrt(p)) * sqrt(qchisq(CEPlevel, df=p))
    RMSE <- if(accuracy) {
        ## POA != POI -> essentially Rice (2D) / offset sphere (3D)
        RMSE2CEP * RMSE_AT
    } else {
        ## POA = POI -> essentially Rayleigh (2D) / Maxwell-Boltzmann(3D)
        RMSE2CEP * RMSE_AF
    }

    setNames(RMSE, CEPlevel)
}

#####-----------------------------------------------------------------------
## Ethridge CEP estimate from Ethridge (1983) after Puhek (1992)
CEPEthridge <-
function(CEPlevel=0.5, accuracy=FALSE, xy) {
    lnDTC <- if(accuracy) {              # log distance to group center (radius)
        rSqSum <- sqrt(rowSums(xy^2))    # log radii to origin = point of aim
        log(rSqSum)
    } else {
        log(getDistToCtr(xy))            # log radii to group center
    }

    mLnDTC   <- mean(lnDTC)              # mean log radius
    medLnDTC <- median(lnDTC)            # median log radius
    varLnDTC <- var(lnDTC)               # variance log radius

    ## weighted mean after Hogg (1967)
    ## sample kurtosis log radius
    kLnDTC <- mean((lnDTC - mLnDTC)^4) / mean((lnDTC - mLnDTC)^2)^2
    dHogg  <- pmax(1 + (0.03 * (kLnDTC-3)^3 * (lnDTC-medLnDTC)^2 / varLnDTC), 0.01)
    wHogg  <- (1/dHogg) / sum(1/dHogg)   # weighting factors
    uHogg  <- sum(wHogg * lnDTC)         # log median radius estimate
    Ethridge50 <- exp(uHogg)

    if(any(CEPlevel != 0.5)) {
        warning("Ethridge CEP estimate is only available for CEPlevel 0.5")
    }

    Ethridge <- ifelse(CEPlevel == 0.5, Ethridge50, NA_real_)
    setNames(Ethridge, CEPlevel)
}

#####-----------------------------------------------------------------------
## modified RAND-234 CEP estimate for 50% from Williams, 1997
## using the semi-major and semi-minor axes of the error ellipse (PCA)
CEPRAND <-
function(CEPlevel=0.5, ctr=c(0, 0), sigma=diag(length(ctr)), accuracy=FALSE) {
    p <- ncol(sigma)

    ## make sure eigenvalues >= 0 when very small
    ev     <- eigen(sigma)$values    # eigenvalues
    lambda <- ev*sign(ev)

    RAND50 <- if(p == 2L) {          # only available for 2D case
        RAND50MPI <- 0.5620*sqrt(lambda[1]) + 0.6152*sqrt(lambda[2])
        if(accuracy) {        # take systematic location bias into account
            bias <- sqrt(sum(ctr^2)) / RAND50MPI
            if(bias > 2.2) {
                warning(c("RAND location bias estimate is ",
                          round(bias, 2), " (> 2.2),\n",
                          "more than what RAND CEP should be considered for"))
            }

            ## cubic regression to take bias into account
            RAND50MPI * (1.0039 - 0.0528*bias + 0.4786*bias^2 - 0.0793*bias^3)
        } else {                         # ignore location bias
            RAND50MPI
        }                                # if(accuracy)
    } else {                       # 1D/3D case
        warning("RAND CEP estimate is only available for 2D-data")
        rep(NA_real_, length(CEPlevel))
    }

    if(any(CEPlevel != 0.5)) {
        warning("RAND CEP estimate is only available for CEPlevel 0.5")
    }

    RAND <- ifelse(CEPlevel == 0.5, RAND50, NA_real_)
    setNames(RAND, CEPlevel)
}

#####-----------------------------------------------------------------------
## Valstar CEP estimate for 50% from Williams, 1997
## using the semi-major and semi-minor axes of the error ellipse (PCA)
CEPValstar <-
function(CEPlevel=0.5, ctr=c(0, 0), sigma=diag(length(ctr)), accuracy=FALSE) {
    p <- ncol(sigma)
    aspRat <- sqrt(kappa(sigma, exact=TRUE))

    ## make sure eigenvalues >= 0 when very small
    ev     <- eigen(sigma)$values    # eigenvalues
    lambda <- ev*sign(ev)

    Valstar50 <- if(p == 2L) {       # only available for 2D case
        ValstarMPI <- if((1/aspRat) <= 0.369) {
            0.675*sqrt(lambda[1]) + sqrt(lambda[2])/(1.2*sqrt(lambda[1]))
        } else {
            0.5620*sqrt(lambda[1]) + 0.6152*sqrt(lambda[2])  # almost RAND
        }

        if(accuracy) {               # POA != POI
            Valstar <- sqrt(ValstarMPI^2 + sum(ctr^2))
        } else {                     # POA = POI
            ValstarMPI
        }                            # if(accuracy)
    } else {                         # 1D/3D case
        warning("Valstar CEP estimate is only available for 2D-data")
        rep(NA_real_, length(CEPlevel))
    }

    if(any(CEPlevel != 0.5)) {
        warning("Valstar CEP estimate is only available for CEPlevel 0.5")
    }

    Valstar <- ifelse(CEPlevel == 0.5, Valstar50, NA_real_)
    setNames(Valstar, CEPlevel)
}

# ## Siouris, GM. 1993. Appendix A
# CEPSiouris <- function(sigma) {
#     p <- ncol(sigma)
#     aspRat <- sqrt(kappa(sigma, exact=TRUE))
#
#     ## make sure eigenvalues >= 0 when very small
#     ev     <- eigen(sigma)$values    # eigenvalues
#     lambda <- ev*sign(ev)
#
#     if(p == 2L) {                # 2D
#         theta <- 0.5 * atan(2*sqrt(lambda[1]*lambda[2]) / (lambda[1] - lambda[2]))
#         0.589*(sqrt(lambda[1]*cos(theta)^2 + lambda[2]*sin(theta)^2) +
#                sqrt(lambda[1]*sin(theta)^2 + lambda[2]*cos(theta)^2))
#     } else if(p == 3L) {         # 3D
#         varTotal <- sum(lambda)
#         V <- 2*sum(lambda^2) / varTotal^2
#         sqrt(varTotal*(1-(V/9)^3))
#     }
# }
#
# ## Shultz, ME. 1963. Circular error probability of a quantity affected by a bias
# CEPShultz <- function(ctr, sigma) {
#     p <- ncol(sigma)
#     aspRat <- sqrt(kappa(sigma, exact=TRUE))
#
#     ## make sure eigenvalues >= 0 when very small
#     ev     <- eigen(sigma)$values    # eigenvalues
#     lambda <- ev*sign(ev)
#
#     muH <- sqrt(sum(ctr^2))
#     sdC <- mean(sqrt(lambda))
#     CE90 <- 2.1272*sdC + 0.1674*muH + 0.3623*(muH^2/sdC) - 0.055*(muH^3/sdC)
# }
#
# ## Ager, TP. 2004. An Analysis of Metric Accuracy Definitions
# ## and Methods of Computation NIMA InnoVision
# CEPAger <- function(ctr, sigma) {
#     p <- ncol(sigma)
#     aspRat <- sqrt(kappa(sigma, exact=TRUE))
#
#     ## make sure eigenvalues >= 0 when very small
#     ev     <- eigen(sigma)$values    # eigenvalues
#     lambda <- ev*sign(ev)
#
#     muH <- sqrt(sum(ctr^2))
#     sdC <- mean(sqrt(lambda))
#
#     CE90 <- if(muH/sdC <= 0.1) {
#         2.1460*sdC
#     } else if(muH/sdC <= 3) {
#         2.1272*sdC + 0.1674*muH + 0.3623*(muH^2/sdC) - 0.055*(muH^3/sdC)
#     } else {
#         0.9860*muH + 1.4548*sdC
#     }
# }
