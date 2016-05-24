getHitProb <-
function(xy, r=1, unit="unit", dstTarget=100, conversion="m2cm",
         accuracy=FALSE, type="CorrNormal", doRob=FALSE) {
    UseMethod("getHitProb")
}

getHitProb.data.frame <-
function(xy, r=1, unit="unit", dstTarget=100, conversion="m2cm",
         accuracy=FALSE, type="CorrNormal", doRob=FALSE) {
    xy <- getXYmat(xy, xyTopLeft=FALSE, relPOA=FALSE)
    NextMethod("getHitProb")
}

getHitProb.default <-
function(xy, r=1, unit="unit", dstTarget=100, conversion="m2cm",
         accuracy=FALSE, type="CorrNormal", doRob=FALSE) {
    xy <- as.matrix(xy)
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(!is.numeric(r))  { stop("r must be numeric") }
    r <- r[r > 0]
    if(length(r) == 0L) { stop("r must be > 0") }

    unit <- match.arg(unit,
                      choices=c("unit", "m", "cm", "mm", "yd", "ft", "in",
                                "deg", "MOA", "SMOA", "rad", "mrad", "mil"))

    type <- match.arg(type,
                      choices=c("CorrNormal", "GrubbsPearson", "GrubbsPatnaik",
                                "GrubbsLiu", "Rayleigh"), several.ok=TRUE)

    p <- ncol(xy)

    ## check if we can do robust estimation if so required
    haveRob <- if(nrow(xy) < 4L) {
        if(doRob) { warning("We need >= 4 points for robust estimations") }
        FALSE
    } else {
        TRUE
    }                                    # if(nrow(xy) < 4L)

    #####-----------------------------------------------------------------------
    ## some basic calculations used later
    if(doRob && haveRob) {        # center
        rob   <- robustbase::covMcd(xy, cor=FALSE)
        ctr   <- rob$center              # robust estimate: group center
        sigma <- rob$cov                 # robust estimate: covariance matrix
    } else {
        ctr   <- colMeans(xy)
        sigma <- cov(xy)
    }

    ## error ellipse characteristics -> radii = sqrt of eigenvalues
    ## aspect ratio of ellipse = sqrt of kappa condition index
    aspRat <- sqrt(kappa(sigma, exact=TRUE))
    flat   <- 1 - (1/aspRat)             # flattening

    ## infer (x,y)-coord units from conversion
    unitXY <- getUnits(conversion, first=FALSE)

    ## convert r to unit of (x,y)-coordinates
    rNew <- if(unit == "unit") {         # keep unit
        r                                # new r = r
    } else if(unit %in% c("deg", "MOA", "SMOA", "rad", "mrad", "mil")) {
        fromMOA(r, dst=dstTarget, conversion=conversion, type=unit)
    } else {                             # absolute size unit
        r2rNew <- getConvFac(paste0(unit, "2", unitXY))
        r2rNew * r
    }

    #####-----------------------------------------------------------------------
    ## estimate based on correlated bivariate normal distribution
    CorrNormal <- if("CorrNormal" %in% type) {
        if(p == 1L) {                        # 1D
            SD <- sqrt(sigma[1, 1])
            if(accuracy) {                   # offset interval +/- rNew
                pnorm(rNew, mean=ctr, sd=SD) - pnorm(-rNew, mean=ctr, sd=SD)
            } else {                         # centered interval +/- rNew
                pnorm(rNew, mean=0,   sd=SD) - pnorm(-rNew, mean=0,   sd=SD)
            }
        } else if(p == 2L) {                 # 2D
            if(accuracy) {                   # offset circle probability -> mvnEll.R
                pmvnEll(r=rNew, sigma=sigma, mu=numeric(p), x0=ctr, e=diag(p))
            } else {                         # Hoyt distribution -> hoyt.R
                HP <- getHoytParam(sigma)
                pHoyt(rNew, qpar=HP$q, omega=HP$omega)
            }
        } else {                             # 3D -> mvnEll.R
            if(accuracy) {                   # offset sphere probability
                pmvnEll(r=rNew, sigma=sigma, mu=numeric(p), x0=ctr,        e=diag(p))
            } else {                         # centered sphere probability
                pmvnEll(r=rNew, sigma=sigma, mu=numeric(p), x0=numeric(p), e=diag(p))
            }
        }
    } else { rep(NA_real_, length(r)) }

    #####-----------------------------------------------------------------------
    ## Grubbs-Pearson CEP estimate based on Pearson three-moment central
    ## chi^2 approximation (Grubbs, 1964, p55-56)
    GPP <- getGrubbsParam(sigma, ctr=ctr, accuracy=accuracy)
    GrubbsPearson <- if("GrubbsPearson" %in% type) {
        pChisqGrubbs(rNew, m=GPP$m, v=GPP$v, nPrime=GPP$nPrime, type="Pearson")
    } else { rep(NA_real_, length(r)) }

    #####-----------------------------------------------------------------------
    ## Grubbs-Patnaik CEP estimate based on Patnaik two-moment central
    ## chi^2 approximation (Grubbs, 1964, p54)
    GrubbsPatnaik <- if("GrubbsPatnaik" %in% type) {
        GrubbsPatnaik <- pChisqGrubbs(rNew, m=GPP$m, v=GPP$v, n=GPP$n, type="Patnaik")
    } else { rep(NA_real_, length(r)) }

    #####-----------------------------------------------------------------------
    ## Grubbs-Liu CEP estimate based on four-moment non-central chi^2
    ## approximation (Liu, Tang & Zhang, 2009)
    GrubbsLiu <- if("GrubbsLiu" %in% type) {
        pChisqGrubbs(rNew, m=GPP$m, v=GPP$v, muX=GPP$muX, varX=GPP$varX,
                           l=GPP$l, delta=GPP$delta, type="Liu")
    } else { rep(NA_real_, length(r)) }

    #####-----------------------------------------------------------------------
    ## Rayleigh estimate from Singh
    Rayleigh <- if("Rayleigh" %in% type) {
        if(aspRat > 2) {
            warning(c("Aspect ratio of error ellipse is ",
                      round(aspRat, 2) , " (> 2),\n",
                      "probably more than what Rayleigh distribution should be considered for"))
        }

        if(p == 1L) {                        # 1D
            if(accuracy) {                   # POA != POI
                RiceParam <- getRiceParam(xy, doRob=doRob)
                NU <- RiceParam$nu
                SD <- RiceParam$sigma["sigma"]
                pnorm(rNew, mean=NU, sd=SD) - pnorm(-rNew, mean=NU, sd=SD)
            } else {                         # POA = POI -> half normal distribution
                HalfNormParam <- getRayParam(xy, doRob=doRob)
                SD <- HalfNormParam$sigma["sigma"]
                pnorm(rNew, mean=0, sd=SD) - pnorm(-rNew, mean=0, sd=SD)
            }
        } else if(p == 2L) {                 # 2D
            if(accuracy) {                   # POA != POI -> Rice distribution
                RiceParam <- getRiceParam(xy, doRob=doRob)
                pRice(rNew, nu=RiceParam$nu, sigma=RiceParam$sigma["sigma"])
            } else {                         # POA = POI -> Rayleigh distribution
                RayParam <- getRayParam(xy, doRob=doRob)
                pRayleigh(rNew, scale=RayParam$sigma["sigma"])
            }
        } else if(p == 3L) {                 # 3D
            MaxParam <- getRayParam(xy, doRob=doRob)
            if(accuracy) {                   # offset circle probability
                ## circular covariance matrix with estimated M-B param sigma
                sigMat <- diag(rep(MaxParam$sigma["sigma"]^2, times=p))
                pmvnEll(r=rNew, sigma=sigMat, mu=numeric(p), x0=ctr, e=diag(p))
            } else {                         #  -> Maxwell-Boltzmann distribution
                pMaxwell(rNew, sigma=MaxParam$sigma["sigma"])
            }
        }
    }

    #####-----------------------------------------------------------------------
    ## only report the chosen estimates
    CEPinvMat <- cbind(CorrNormal=CorrNormal,
                       GrubbsPearson=GrubbsPearson,
                       GrubbsPatnaik=GrubbsPatnaik,
                       GrubbsLiu=GrubbsLiu,
                       Rayleigh=Rayleigh)

    rownames(CEPinvMat) <- paste0("R", r)

    return(CEPinvMat[ , type, drop=TRUE])
}
