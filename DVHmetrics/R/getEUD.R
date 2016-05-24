## we need ... because getMetric() will also pass parameters
## intended for other functions through ...
getEUD <-
function(x, EUDa, EUDfd=NULL, EUDab=NULL, ...) {
    UseMethod("getEUD")
}

getEUD.DVHs <-
function(x, EUDa, EUDfd=NULL, EUDab=NULL, ...) {
    if(length(EUDa) > 1) {
		warning(paste0("Will only use EUDa=", EUDa[1]))
		EUDa <- EUDa[1]
	}

    if(isTRUE(all.equal(EUDa, 0))) {
		warning("'EUDa' must not be zero")
		return(NA_real_)
	}

    ## special cases
	if(isTRUE(all.equal(EUDa, 1))) {
		return(getMetric(x, "DMEAN")$DMEAN)
	} else if(is.infinite(EUDa) && (EUDa > 0)) {  # +Inf
		return(getMetric(x, "DMAX")$DMAX)
	} else if(is.infinite(EUDa) && (EUDa < 0)) {  # -Inf
		return(getMetric(x, "DMIN")$DMIN)
	}

    ## get differential DVH
    xD <- convertDVH(x, toType="differential", toDoseUnit="asis", perDose=FALSE)

    ## convert dose to EQD2 if possible
    volume <- xD$dvhDiff[ , "volume"]
    dose   <- if(!is.null(EUDfd) && !is.null(EUDab)) {
        getEQD2(D=xD$dvhDiff[ , "dose"], fd=EUDfd, ab=EUDab)$EQD2
    } else {
        xD$dvhDiff[ , "dose"]
    }

    volDose <- volume*dose^EUDa / xD$structVol
    wtMean  <- sum(volDose[volume > 0], na.rm=TRUE)
    gEUD    <- wtMean^(1/EUDa)
	if(!is.finite(gEUD)) {
		warning("Numerical problems encountered, NA returned")
		gEUD <- NA_real_
	}

    data.frame(EUD=gEUD,
               patID=x$patID,
               structure=x$structure,
               stringsAsFactors=FALSE)
}

getEUD.DVHLst <-
function(x, EUDa, EUDfd=NULL, EUDab=NULL, ...) {
    EUDl  <- Map(getEUD, x,
                 EUDa=list(EUDa), EUDfd=list(EUDfd), EUDab=list(EUDab))
    EUDdf <- do.call("rbind", EUDl)
    rownames(EUDdf) <- NULL
    EUDdf
}

getEUD.DVHLstLst <-
function(x, EUDa, EUDfd=NULL, EUDab=NULL, ...) {
    EUDl  <- Map(getEUD, x,
                 EUDa=list(EUDa), EUDfd=list(EUDfd), EUDab=list(EUDab))
    EUDdf <- do.call("rbind", EUDl)
    rownames(EUDdf) <- NULL
    EUDdf
}
