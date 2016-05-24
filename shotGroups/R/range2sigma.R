## asumme Rayleigh case
## convert extreme spread / figure of merit / bounding box diagonal
## to Rayleigh sigma using lookup table from 1000000 runs for each
## combination of n*nGroups
## http://ballistipedia.com/index.php?title=Range_Statistics
range2sigma <-
    function(x, stat="ES", n=5, nGroups=1, CIlevel=0.95, collapse=TRUE,
             dstTarget=100, conversion="m2cm") {
    n       <- as.integer(n[1])
    nGroups <- as.integer(nGroups[1])
    stopifnot(all(x > 0),
              n       > 1L, n       <= max(shotGroups::DFdistr$n),
              nGroups > 0L, nGroups <= max(shotGroups::DFdistr$nGroups),
              CIlevel > 0)
    stat <- match.arg(toupper(stat), choices=c("ES", "FOM", "D"), several.ok=TRUE)

    argL <- recycle(x, stat)
    x    <- argL[[1]]
    stat <- argL[[2]]
    x    <- setNames(x, stat)

    ## check if CIlevel is given in percent
    if(CIlevel >= 1) {
        while(CIlevel >= 1) { CIlevel <- CIlevel / 100 }
        warning(c("CIlevel must be in (0,1) and was set to ", CIlevel))
    }

    alpha    <- 1 - CIlevel
    idxGroup <- which(shotGroups::DFdistr$nGroups == nGroups)
    haveN    <- unique(shotGroups::DFdistr$n[idxGroup])
    haveCI   <- c(0.5, 0.8, 0.9, 0.95, 0.99)

    ## Rayleigh sigma: use lookup table or monotone spline interpolation
    ## CI: use lookup table or Grubbs-Patnaik chi^2 approximation
    if(n %in% haveN) {
        ## can use lookup table for sigma estimate
        idx <- which((shotGroups::DFdistr$n == n) & (shotGroups::DFdistr$nGroups == nGroups))
        M   <- setNames(numeric(length(x)), stat)
        M[names(M) %in% "ES"]  <- shotGroups::DFdistr$ES_M[idx]
        M[names(M) %in% "FOM"] <- shotGroups::DFdistr$FoM_M[idx]
        M[names(M) %in% "D"]   <- shotGroups::DFdistr$D_M[idx]

        ## for Grubbs-Patnaik ES-CI approximation
        if(!(CIlevel %in% haveCI)) {
            ES_V    <- shotGroups::DFdistr$ES_V[idx]
            ESSQ_M  <- shotGroups::DFdistr$ESSQ_M[idx]
            ESSQ_V  <- shotGroups::DFdistr$ESSQ_V[idx]
            ES_SKEW <- shotGroups::DFdistr$ES_SKEW[idx]
            ES_KURT <- shotGroups::DFdistr$ES_KURT[idx]
        }
    } else {
        ## spline interpolation for sigma estimate - M is monotonically increasing
        warning("Rayleigh sigma estimate based on monotone spline interpolation for n")

        ## only interpolate for requested number of groups
        M <- setNames(numeric(length(x)), stat)
        M[names(M) %in% "ES"]  <- splinefun(haveN, shotGroups::DFdistr$ES_M[idxGroup],  method="monoH.FC")(n)
        M[names(M) %in% "FOM"] <- splinefun(haveN, shotGroups::DFdistr$FoM_M[idxGroup], method="monoH.FC")(n)
        M[names(M) %in% "D"]   <- splinefun(haveN, shotGroups::DFdistr$D_M[idxGroup],   method="monoH.FC")(n)

        ## spline interpolation for Grubbs-Patnaik ES-CI approximation
        ## M is monotonically increasing, but V, SKEW, KURT are not
        if(!(CIlevel %in% haveCI)) {
            ES_V    <- splinefun(haveN, shotGroups::DFdistr$ES_V[idxGroup],    method="fmm")(n)
            ESSQ_M  <- splinefun(haveN, shotGroups::DFdistr$ESSQ_M[idxGroup],  method="monoH.FC")(n)
            ESSQ_V  <- splinefun(haveN, shotGroups::DFdistr$ESSQ_V[idxGroup],  method="fmm")(n)
            #ES_SKEW <- splinefun(haveN, shotGroups::DFdistr$ES_SKEW[idxGroup], method="fmm")(n)
            #ES_KURT <- splinefun(haveN, shotGroups::DFdistr$ES_KURT[idxGroup], method="fmm")(n)
        }
    }

    ## Rayleigh sigma estimate
    sigma <- x/M

    ## need extreme spread for sigma CI
    xES  <- setNames(x[names(x) %in% "ES"],  NULL)
    xFoM <- setNames(x[names(x) %in% "FOM"], NULL)
    xD   <- setNames(x[names(x) %in% "D"],   NULL)
    if(CIlevel %in% haveCI) {
        CIlo <- sprintf("%04.1f", round((1-(alpha/2))*100, digits=1))
        CIhi <- sprintf("%04.1f", round(   (alpha/2) *100, digits=1))
        if(n %in% haveN) {
            ES_CIlo  <-  xES/shotGroups::DFdistr[[sub("\\.", "", paste0("ES_Q",  CIlo))]][idx]
            ES_CIhi  <-  xES/shotGroups::DFdistr[[sub("\\.", "", paste0("ES_Q",  CIhi))]][idx]
            FoM_CIlo <- xFoM/shotGroups::DFdistr[[sub("\\.", "", paste0("FoM_Q", CIlo))]][idx]
            FoM_CIhi <- xFoM/shotGroups::DFdistr[[sub("\\.", "", paste0("FoM_Q", CIhi))]][idx]
            D_CIlo   <-   xD/shotGroups::DFdistr[[sub("\\.", "", paste0("D_Q",   CIlo))]][idx]
            D_CIhi   <-   xD/shotGroups::DFdistr[[sub("\\.", "", paste0("D_Q",   CIhi))]][idx]
        } else {
            ## TODO this does not work
            ES_CIlo  <-  xES/splinefun(haveN, shotGroups::DFdistr[[sub("\\.", "", paste0("ES_Q",  CIlo))]][idxGroup],
                                       method="fmm")(n)
            ES_CIhi  <-  xES/splinefun(haveN, shotGroups::DFdistr[[sub("\\.", "", paste0("ES_Q",  CIhi))]][idxGroup],
                                       method="fmm")(n)
            FoM_CIlo <- xFoM/splinefun(haveN, shotGroups::DFdistr[[sub("\\.", "", paste0("FoM_Q", CIlo))]][idxGroup],
                                       method="fmm")(n)
            FoM_CIhi <- xFoM/splinefun(haveN, shotGroups::DFdistr[[sub("\\.", "", paste0("FoM_Q", CIhi))]][idxGroup],
                                       method="fmm")(n)
            D_CIlo   <-   xD/splinefun(haveN, shotGroups::DFdistr[[sub("\\.", "", paste0("D_Q",   CIlo))]][idxGroup],
                                       method="fmm")(n)
            D_CIhi   <-   xD/splinefun(haveN, shotGroups::DFdistr[[sub("\\.", "", paste0("D_Q",   CIhi))]][idxGroup],
                                       method="fmm")(n)
        }
    } else {
        if("ES" %in% stat) {
            warning("CI estimate based on Grubbs-Patnaik chi^2 approximation")
            ES_M <- setNames(M[names(M) %in% "ES"], NULL)
            ## Patnaik chi-approximation
            m <- ESSQ_M
            v <- ESSQ_V
            # m <- ES_M^2 + ES_V
            # v <- ES_KURT*ES_V^2 + 4*ES_SKEW*sqrt(ES_V)^3*ES_M + 4*ES_V*ES_M^2 - ES_V^2
            ES_CIlo  <- xES/qChisqGrubbs(1-(alpha/2), m=m, v=v, n=2*m^2/v)
            ES_CIhi  <- xES/qChisqGrubbs(   alpha/2,  m=m, v=v, n=2*m^2/v)
            FoM_CIlo <- NULL
            FoM_CIhi <- NULL
            D_CIlo   <- NULL
            D_CIhi   <- NULL
        } else {
            warning("CI from Grubbs-Patnaik approximation requires extreme spread")
            ES_CIlo  <- NULL
            ES_CIhi  <- NULL
            FoM_CIlo <- NULL
            FoM_CIhi <- NULL
            D_CIlo   <- NULL
            D_CIhi   <- NULL
        }
    }

    ## convert CIs to MOA
    sigma <- makeMOA(sigma, dst=dstTarget, conversion=conversion)
    sigmaESCI  <- lapply(seq_along(ES_CIlo),  function(i) {
        makeMOA(c("sigma ("=ES_CIlo[i], "sigma )" =ES_CIhi[i]),
                dst=dstTarget, conversion=conversion) })
    sigmaESCI  <- if(length(sigmaESCI) > 0L) {
        setNames(sigmaESCI,  paste0("ES",  round(xES, digits=2)))
    } else { NULL }

    sigmaFoMCI <- lapply(seq_along(FoM_CIlo), function(i) {
        makeMOA(c("sigma ("=FoM_CIlo[i], "sigma )"=FoM_CIhi[i]),
                dst=dstTarget, conversion=conversion) })
    sigmaFoMCI <- if(length(sigmaFoMCI) > 0L) {
        setNames(sigmaFoMCI, paste0("FoM", round(xFoM, digits=2)))
    } else { NULL }

    sigmaDCI   <- lapply(seq_along(D_CIlo),   function(i) {
        makeMOA(c("sigma ("=D_CIlo[i], "sigma )"  =D_CIhi[i]),
                dst=dstTarget, conversion=conversion) })
    sigmaDCI   <- if(length(sigmaDCI) > 0L) {
        setNames(sigmaDCI,   paste0("D",   round(xD,   digits=2)))
    } else { NULL }

    ## weed out non existing CIs
    sigmaCI <- Filter(Negate(is.null), list(ES=sigmaESCI, FoM=sigmaFoMCI, D=sigmaDCI))

    ## collapse sigmaCI list if required and possible
    if(collapse) {
        for(i in seq_along(sigmaCI)) {
            if(length(sigmaCI[[i]]) == 1L) { sigmaCI[[i]] <- sigmaCI[[c(i, 1)]] }
        }

        if(length(sigmaCI) == 1L) { sigmaCI <- sigmaCI[[1]] }
    }

    ## sigmaCI might be empty list when no ES is given
    Filter(function(l) { !is.null(l) && (length(l) > 0L) },
           list(sigma=sigma, sigmaCI=sigmaCI))
}
