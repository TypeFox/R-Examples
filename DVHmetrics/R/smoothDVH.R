## DVH smoothing
## ... in getKSmooth(), getSmoothSpl(), getInterpLin() is ignored
## and used to catch method="FMM" or method="monoH.FC" meant for getInterpSpl()
## when called from convertDVH()

## function for cubic local polynomial smoothing
## dose, dose rel, volume, volume rel, N DVH nodes
getKSmooth <- function(d, dR, v, vR, nodes=NULL, rangeD=NULL, ...) {
    nodes <- if(is.null(nodes)) {
        length(d)
    } else {
        nodes
    }
    
    ## smooth
    smDV <- if(!all(is.na(v))) {
        bwDV <- try(KernSmooth::dpill(d, v))
        try(KernSmooth::locpoly(d, v,  bandwidth=bwDV,  gridsize=nodes, degree=3))
    } else {
        NA_real_
    }

    smDVR <- if(!all(is.na(vR))) {
        bwDVR <- try(KernSmooth::dpill(d, vR))
        try(KernSmooth::locpoly(d, vR, bandwidth=bwDVR, gridsize=nodes, degree=3))
    } else {
        NA_real_
    }

    ## dose
    dose <- if(!inherits(smDV, "try-error") && !is.na(smDV)) {
        smDV$x
    } else if(!inherits(smDVR, "try-error") && !is.na(smDVR)) {
        smDVR$x
    } else {
        NA_real_
    }

    ## dose rel -> just use equally spaced grid points as is done in locpoly()
    doseRel <- if(!any(is.na(dR))) {
        rangeD <- range(dR)
        seq(rangeD[1], rangeD[2], length.out=nodes)
    } else {
        NA_real_
    }

    volume <- if(!inherits(smDV, "try-error") && !is.na(smDV)) {
        smDV$y
    } else {
        NA_real_
    }

    volumeRel <- if(!inherits(smDVR, "try-error") && !is.na(smDVR)) {
        smDVR$y
    } else {
        NA_real_
    }

    return(cbind(dose=dose, doseRel=doseRel, volume=volume, volumeRel=volumeRel))
}

## function for cubic local polynomial smoothing
## dose, dose rel, volume, volume rel, N DVH nodes
getSmoothSpl <- function(d, dR, v, vR, nodes=NULL, rangeD=NULL, ...) {
    nodes <- if(is.null(nodes)) {
        length(d)
    } else {
        nodes
    }

    if(is.null(rangeD)) { rangeD  <- range(d) }

    ## smooth
    smDV <- if(!all(is.na(v))) {
        try(smooth.spline(d, v))
    } else {
        NA_real_
    }

    smDVR <- if(!all(is.na(vR))) {
        try(smooth.spline(d, vR))
    } else {
        NA_real_
    }

    ## dose, dose rel -> just use equally spaced grid points
    dose    <- seq(rangeD[1],  rangeD[2], length.out=nodes)
    doseRel <- if(!any(is.na(dR))) {
        rangeDR <- range(dR)
        seq(rangeDR[1], rangeDR[2], length.out=nodes)
    } else {
        NA_real_
    }

    volume <- if(!inherits(smDV, "try-error") && !is.na(smDV)) {
        predict(smDV,  dose)$y
    } else {
        NA_real_
    }

    volumeRel <- if(!inherits(smDVR, "try-error") && !is.na(smDVR)) {
        predict(smDVR, dose)$y
    } else {
        NA_real_
    }

    return(cbind(dose=dose, doseRel=doseRel, volume=volume, volumeRel=volumeRel))
}

## function for cubic spline interpolation
## dose, dose rel, volume, volume rel, N DVH nodes
getInterpSpl <- function(d, dR, v, vR, nodes=NULL, rangeD=NULL,
                         method=c("fmm", "monoH.FC")) {
    method <- match.arg(method)

    nodes <- if(is.null(nodes)) {
        length(d)
    } else {
        nodes
    }

    if(is.null(rangeD)) { rangeD  <- range(d) }
    
    ## interpolation
    smDV <- if(!all(is.na(v))) {
        try(splinefun(d, v, method=method))
    } else {
        function(x) { return(NA_real_) }
    }

    smDVR <- if(!all(is.na(vR))) {
        try(splinefun(d, vR, method=method))
    } else {
        function(x) { return(NA_real_) }
    }
    
    ## dose, dose rel -> just use equally spaced grid points
    dose    <- seq(rangeD[1],  rangeD[2], length.out=nodes)
    doseRel <- if(!any(is.na(dR))) {
        rangeDR <- range(dR)
        seq(rangeDR[1], rangeDR[2], length.out=nodes)
    } else {
        NA_real_
    }

    volume <- if(!inherits(smDV, "try-error")) {
        smDV(dose)
    } else {
        NA_real_
    }

    volumeRel <- if(!inherits(smDVR, "try-error")) {
        smDVR(dose)
    } else {
        NA_real_
    }

    return(cbind(dose=dose, doseRel=doseRel, volume=volume, volumeRel=volumeRel))
}

## function for linear interpolation
## dose, dose rel, volume, volume rel, N DVH nodes
getInterpLin <- function(d, dR, v, vR, nodes=NULL, rangeD=NULL, ...) {
    nodes <- if(is.null(nodes)) {
        length(d)
    } else {
        nodes
    }

    if(is.null(rangeD)) {
        rangeD <- range(d)
    }
    
    ## interpolation
    smDV <- if(!all(is.na(v))) {
        try(approxfun(d, v,  method="linear", rule=2, ties=max))
    } else {
        function(x) { return(NA_real_) }
    }

    smDVR <- if(!all(is.na(vR))) {
        try(approxfun(d, vR, method="linear", rule=2, ties=max))
    } else {
        function(x) { return(NA_real_) }
    }
    
    ## dose, dose rel -> just use equally spaced grid points
    dose    <- seq(rangeD[1],  rangeD[2], length.out=nodes)
    doseRel <- if(!any(is.na(dR))) {
        rangeDR <- range(dR)
        seq(rangeDR[1], rangeDR[2], length.out=nodes)
    } else {
        NA_real_
    }

    volume <- if(!inherits(smDV, "try-error")) {
        smDV(dose)
    } else {
        NA_real_
    }

    volumeRel <- if(!inherits(smDVR, "try-error")) {
        smDVR(dose)
    } else {
        NA_real_
    }

    return(cbind(dose=dose, doseRel=doseRel, volume=volume, volumeRel=volumeRel))
}
