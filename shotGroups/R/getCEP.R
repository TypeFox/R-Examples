## TODO:
## return: CEP is list with one component per CEPlevel
## check RMSE for 3D data
## check Valstar

getCEP <-
function(xy, CEPlevel=0.5, dstTarget=100, conversion="m2cm", accuracy=FALSE,
         type="CorrNormal", doRob=FALSE) {
    UseMethod("getCEP")
}

getCEP.data.frame <-
function(xy, CEPlevel=0.5, dstTarget=100, conversion="m2cm", accuracy=FALSE,
         type="CorrNormal", doRob=FALSE) {
    xy <- getXYmat(xy, xyTopLeft=FALSE, relPOA=FALSE)
    NextMethod("getCEP")
}

getCEP.default <-
function(xy, CEPlevel=0.5, dstTarget=100, conversion="m2cm", accuracy=FALSE,
         type="CorrNormal", doRob=FALSE) {
    xy <- as.matrix(xy)
    if(!is.numeric(xy))        { stop("xy must be numeric") }
    if(!is.numeric(CEPlevel))  { stop("CEPlevel must be numeric") }
    CEPlevel <- CEPlevel[CEPlevel > 0]
    if(length(CEPlevel) == 0L) { stop("CEPlevel must be > 0") }

    CEPlevel <- sort(unique(CEPlevel))

    type <- match.arg(type,
                      choices=c("CorrNormal", "GrubbsPearson", "GrubbsLiu",
                                "GrubbsPatnaik", "Rayleigh", "Krempasky",
                                "Ignani", "RMSE", "Ethridge", "RAND", "Valstar"),
                      several.ok=TRUE)

    ## check if CEPlevel is given in percent
    CEPlevelBefore <- CEPlevel
    CEPlevel <- vapply(CEPlevel, function(x) {
        if(x >= 1) { while(x >= 1) { x <- x / 100 } } else { x } }, numeric(1))

    if(any(CEPlevelBefore != CEPlevel)) {
        warning(c("CEPlevel must be in (0,1) and was set to ", CEPlevel))
    }

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
        ctr   <- rob$center                       # robust estimate: group center
        sigma <- rob$cov
    } else {
        ctr   <- colMeans(xy)
        sigma <- cov(xy)
    }

    ## error ellipse characteristics -> radii = sqrt of eigenvalues
    ## aspect ratio of ellipse = sqrt of condition index kappa
    aspRat <- sqrt(kappa(sigma, exact=TRUE))
    flat   <- 1 - (1/aspRat)             # flattening

    #####-----------------------------------------------------------------------
    ## CEP estimates
    CorrNormal <- if("CorrNormal" %in% type) {
        CEPCorrNormal(CEPlevel, ctr=ctr, sigma=sigma, accuracy=accuracy)
    } else { rep(NA_real_, length(CEPlevel)) }

    GrubbsPatnaik <- if("GrubbsPatnaik" %in% type) {
        CEPGrubbsPatnaik(CEPlevel, ctr=ctr, sigma=sigma, accuracy=accuracy)
    } else { rep(NA_real_, length(CEPlevel)) }

    GrubbsPearson <- if("GrubbsPearson" %in% type) {
        CEPGrubbsPearson(CEPlevel, ctr=ctr, sigma=sigma, accuracy=accuracy)
    } else { rep(NA_real_, length(CEPlevel)) }

    GrubbsLiu <- if("GrubbsLiu" %in% type) {
        CEPGrubbsLiu(CEPlevel, ctr=ctr, sigma=sigma, accuracy=accuracy)
    } else { rep(NA_real_, length(CEPlevel)) }

    Rayleigh <- if("Rayleigh" %in% type) {
        if((aspRat > 2)) {
            warning(c("Aspect ratio of error ellipse is ",
                     round(aspRat, digits=2) , " (> 2),\n",
                     "probably more than what Rayleigh CEP should be considered for"))
        }

        CEPRayleigh(CEPlevel, ctr=ctr, sigma=sigma, accuracy=accuracy, doRob=doRob, xy=xy)
    } else { rep(NA_real_, length(CEPlevel)) }

    Krempasky <- if("Krempasky" %in% type) {
        CEPKrempasky(CEPlevel, sigma=sigma, accuracy=accuracy)
    } else { rep(NA_real_, length(CEPlevel)) }

    Ignani <- if("Ignani" %in% type) {
        CEPIgnani(CEPlevel, sigma=sigma, accuracy=accuracy)
    } else { rep(NA_real_, length(CEPlevel)) }

    RMSE <- if("RMSE" %in% type) {
        CEPRMSE(CEPlevel, sigma=sigma, accuracy=accuracy, xy=xy)
    } else { rep(NA_real_, length(CEPlevel)) }

    Ethridge <- if("Ethridge" %in% type) {
        CEPEthridge(CEPlevel, accuracy=accuracy, xy=xy)
    } else { rep(NA_real_, length(CEPlevel)) }

    RAND <- if("RAND" %in% type) {
        if((aspRat > 3)) {
            warning(c("Aspect ratio of error ellipse is ",
                      round(aspRat, digits=2) , " (> 3),\n",
                      "probably more than what RAND CEP should be considered for"))
        }

        CEPRAND(CEPlevel, ctr=ctr, sigma=sigma, accuracy=accuracy)
    } else { rep(NA_real_, length(CEPlevel)) }

    Valstar <- if("Valstar" %in% type) {
        if((aspRat > 3)) {
            warning(c("Aspect ratio of error ellipse is ",
                      round(aspRat, digits=2) , " (> 3),\n",
                      "probably more than what Valstar CEP should be considered for"))
        }

        CEPValstar(CEPlevel, ctr=ctr, sigma=sigma, accuracy=accuracy)
    } else { rep(NA_real_, length(CEPlevel)) }

    #####-----------------------------------------------------------------------
    ## only report the chosen estimates
    CEPDF <- data.frame(CorrNormal=CorrNormal,
                        GrubbsPatnaik=GrubbsPatnaik,
                        GrubbsPearson=GrubbsPearson,
                        GrubbsLiu=GrubbsLiu,
                        Rayleigh=Rayleigh,
                        Krempasky=Krempasky,
                        Ignani=Ignani,
                        Ethridge=Ethridge,
                        RMSE=RMSE,
                        RAND=RAND,
                        Valstar=Valstar)[type]

    CEPL <- setNames(split(data.matrix(CEPDF), CEPlevel), paste0("CEP", CEPlevel))
    CEP  <- lapply(CEPL, function(x) {
        M <- as.matrix(makeMOA(x, dst=dstTarget, conversion=conversion))
        colnames(M) <- type
        M
    })

    return(c(CEP=list(CEP),
             ellShape=list(c(aspectRatio=aspRat, flattening=flat)),
             ctr=list(ctr)))
}
