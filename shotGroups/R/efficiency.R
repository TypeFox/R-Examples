## http://ballistipedia.com/index.php?title=Range_Statistics
efficiency <-
function(n, nGroups, CIlevel=0.95, CIwidth,
         stat=c("Rayleigh", "ES", "FoM", "D")) {
    stopifnot(is.numeric(n), is.numeric(CIlevel),
              all(n > 1L), all(n <= max(shotGroups::DFdistr$n)),
              all(CIlevel > 0))
    n       <- sort(unique(as.integer(n)))
    CIlevel <- round(CIlevel[1], digits=2)
    stat    <- match.arg(toupper(stat),
                         choices=c("RAYLEIGH", "ES", "FOM", "D"),
                         several.ok=FALSE)

    ## check if CIlevel is given in percent
    if(CIlevel >= 1) {
        while(CIlevel >= 1) { CIlevel <- CIlevel / 100 }
        warning(c("CIlevel must be in (0,1) and was set to ", CIlevel))
    }
    
    idx   <- which((shotGroups::DFdistr$n %in% n) & (shotGroups::DFdistr$nGroups == 1L))
    idx1  <- which(shotGroups::DFdistr$nGroups == 1L)
    n1    <- shotGroups::DFdistr[idx1, "n", drop=TRUE]
    alpha <- 1 - CIlevel
    z     <- qnorm(1-(alpha/2), mean=0, sd=1)

    ## can use lookup table for ES_CV and RS_CV or do spline interpolation
    if(all(n %in% n1)) {
        nActual <- shotGroups::DFdistr$n[idx]
        ES_CV   <- shotGroups::DFdistr$ES_CV[idx]
        FoM_CV  <- shotGroups::DFdistr$FoM_CV[idx]
        D_CV    <- shotGroups::DFdistr$D_CV[idx]
        RS_CV   <- shotGroups::DFdistr$RS_CV[idx]
    } else {
        ## spline interpolation for ES_CV/FoM_CV/D_CV/RS_CV
        nActual <- n
        ES_CV   <- splinefun(n1, shotGroups::DFdistr[idx1, "ES_CV"],  method="monoH.FC")(n)
        FoM_CV  <- splinefun(n1, shotGroups::DFdistr[idx1, "FoM_CV"], method="monoH.FC")(n)
        D_CV    <- splinefun(n1, shotGroups::DFdistr[idx1, "D_CV"],   method="monoH.FC")(n)
        RS_CV   <- splinefun(n1, shotGroups::DFdistr[idx1, "RS_CV"],  method="monoH.FC")(n)
    }

    m <- if(missing(nGroups) && !missing(CIwidth)) {
        ## nGroups is requested, CI width is given
        stopifnot(CIwidth > 0)

        ## check if CIwidth is given in percent
        if(CIwidth >= 1) {
            while(CIwidth >= 1) { CIwidth <- CIwidth / 100 }
            warning(c("CIwidth must be in (0,1) and was set to ", CIwidth))
        }

        E <- CIwidth/2
        nGroupsReq <- if(stat == "RAYLEIGH") {  # estimate from  Rayleigh sigma
            (z*RS_CV/E)^2
        } else if(stat == "ES") {               # estimate from extreme spread
            (z*ES_CV/E)^2
        } else if(stat == "FOM") {              # estimate from figure of merit
            (z*FoM_CV/E)^2
        } else if(stat == "D") {                # estimate from bounding box diagonal
            (z*D_CV/E)^2
        }

        nGroupsReqCeil <- ceiling(nGroupsReq)
        nShotsReq      <- nGroupsReq*nActual
        nShotsReqCeil  <- ceiling(nGroupsReq)*nActual

        data.frame(n=nActual,
                   nGroupsReq=nGroupsReq,
                   nGroupsReqCeil=nGroupsReqCeil,
                   nShotsReq=nShotsReq,
                   nShotsReqCeil=nShotsReqCeil,
                   CIlevel=CIlevel,
                   CIwidth=CIwidth)
    } else if(!missing(nGroups) && missing(CIwidth)) {
        ## nGroups is given, CI width is requested
        stopifnot(is.numeric(nGroups), nGroups > 0)
        nGroups <- as.integer(nGroups[1])

        E <- if(stat == "RAYLEIGH") {
            z*RS_CV/sqrt(nGroups)
        } else if(stat == "ES") {
            z*ES_CV/sqrt(nGroups)
        } else if(stat == "FOM") {
            z*FoM_CV/sqrt(nGroups)
        } else if(stat == "D") {
            z*D_CV/sqrt(nGroups)
        }

        data.frame(n=nActual,
                   nGroups=nGroups,
                   nShots=nGroups*nActual,
                   CIlevel=CIlevel,
                   CIwidth=2*E)

    } else { stop("One of nGroups or CIwidth must be supplied (but not both)") }

    return(m)
}
