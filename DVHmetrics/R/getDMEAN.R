getDMEAN <-
function(x, interp=c("linear", "spline", "ksmooth", "smoothSpl"), nodes=5001L) {
    UseMethod("getDMEAN")
}

getDMEAN.DVHs <-
function(x, interp=c("linear", "spline", "ksmooth", "smoothSpl"), nodes=5001L) {
    interp <- match.arg(interp)
    nodes  <- max(nodes, nrow(x$dvh))

    ## median from cumulative DVH
    doseMed <- if(interp == "linear") {
        res <- if(!all(is.na(x$dvh[ , "volumeRel"]))) {
            try(approx(x$dvh[ , "volumeRel"], x$dvh[ , "dose"], 50, method="linear", rule=1, ties=max)$y)
        } else {
            NA_real_
        }

        if(!inherits(res, "try-error")) {
            res
        } else {
            NA_real_
        }
    } else if(interp == "spline") {
        sfun <- if(!all(is.na(x$dvh[ , "volumeRel"]))) {
            try(splinefun(x$dvh[ , "volumeRel"], x$dvh[ , "dose"], method="monoH.FC", ties=max))
        } else {
            NA_real_
        }

        if(!is.na(sfun) && !inherits(sfun, "try-error")) {
            sfun(50)
        } else {
            NA_real_
        }
    } else if(interp == "ksmooth") {
        sm <- if(!all(is.na(x$dvh[ , "volumeRel"]))) {
            bw <- try(KernSmooth::dpill(x$dvh[ , "volumeRel"], x$dvh[ , "dose"]))
            try(KernSmooth::locpoly(x$dvh[ , "volumeRel"], x$dvh[ , "dose"],
                                    bandwidth=bw, gridsize=nodes, degree=3))
        } else {
            NA_real_
        }

        if(!is.na(sm) && !inherits(sm, "try-error")) {
            idx <- which.min(abs(sm$x-50))
            sm$y[idx]
        } else {
            NA_real_
        }
    } else if(interp == "smoothSpl") {
        sm <- if(!all(is.na(x$dvh[ , "volumeRel"]))) {
            try(smooth.spline(x$dvh[ , "volumeRel"], x$dvh[ , "dose"]))
        } else {
            NA_real_
        }

        if(!is.na(sm) && !inherits(sm, "try-error")) {
            predict(sm, 50)$y
        } else {
            NA_real_
        }
    } else {
        NA_real_
    }

    ## convert to differential DVH - but not per unit dose
    xDiff <- if(interp == "linear") {
        convertDVH(x, toType="differential", interp=interp, nodes=nodes, perDose=FALSE)
    } else {
        warning("non-linear interpolation of differential DVH not recommended for calculation of DMEAN, DMIN, DMAX, ...")
        convertDVHsmooth(x, toType="differential", interp=interp, nodes=nodes, perDose=FALSE)
    }

    ## dose category mid-points
    doseMidPt <- xDiff$dvhDiff[ , "dose"]
    volRel    <- xDiff$dvhDiff[ , "volumeRel"]

    ## available volume
    volume <- if(!all(is.na(xDiff$dvhDiff[ , "volume"]))) {
        xDiff$dvhDiff[ , "volume"]
    } else {
        xDiff$dvhDiff[ , "volumeRel"]
    }

    doseMin <- min(xDiff$dvhDiff[volume > 0, "dose"])
    doseMax <- max(xDiff$dvhDiff[volume > 0, "dose"])
    doseAvg <- sum(doseMidPt*volRel/100)
    doseSD  <- sqrt(sum(doseMidPt^2*volRel/100) - doseAvg^2)

    ## for mode, abs or rel volume does not matter
    doseMode <- if(!all(is.na(xDiff$dvhDiff[ , "volume"]))) {
        xDiff$dvhDiff[which.max(xDiff$dvhDiff[ , "volume"]),    "dose"]
    } else if(!all(is.na(xDiff$dvhDiff[ , "volumeRel"]))) {
        xDiff$dvhDiff[which.max(xDiff$dvhDiff[ , "volumeRel"]), "dose"]
    } else {
        NA_real_
    }

    doseAvgTPS <- if(!is.null(x$doseAvg)) {
        x$doseAvg
    } else {
        NA_real_
    }

    doseMedTPS <- if(!is.null(x$doseMed)) {
        x$doseMed
    } else {
        NA_real_
    }

    doseMinTPS <- if(!is.null(x$doseMin)) {
        x$doseMin
    } else {
        NA_real_
    }

    doseMaxTPS <- if(!is.null(x$doseMax)) {
        x$doseMax
    } else {
        NA_real_
    }

    metrics <- data.frame(patID=x$patID, structure=x$structure,
                          doseMin=doseMin, doseMax=doseMax, doseAvg=doseAvg,
                          doseMed=doseMed, doseSD=doseSD, doseMode=doseMode,
                          doseAvgTPS=doseAvgTPS, doseMedTPS=doseMedTPS,
                          doseMinTPS=doseMinTPS, doseMaxTPS=doseMaxTPS)

    return(metrics)
}

getDMEAN.DVHLst <-
function(x, interp=c("linear", "spline", "ksmooth", "smoothSpl"), nodes=5001L) {
    interp <- match.arg(interp)
    ml <- Map(getDMEAN, x, interp=interp, nodes=nodes)
    df <- do.call("rbind", ml)
    rownames(df) <- NULL
    df
}

getDMEAN.DVHLstLst <-
function(x, interp=c("linear", "spline", "ksmooth", "smoothSpl"), nodes=5001L) {
    interp <- match.arg(interp)
    ml <- Map(getDMEAN, x, interp=interp, nodes=nodes)
    df <- do.call("rbind", ml)
    rownames(df) <- NULL
    df
}

## mean from differential DVH from
## integration of dose * (spline fit derivative)
#     xD <- convertDVH(x, toType="differential", smooth=interp, nodes=1001L)
#     ## linear
#     lfun <- approxfun(dose, 100-volume, method="linear", rule=2)
#     meanLin <- tryCatch(pracma::quadgk(function(y) {
#         y*numDeriv::grad(lfun, x=y)/100 }, 0, max(dose)),
#         error=function(e) return(NA_real_))
#
#     ## monotone Hermite spline
#     sfun <- try(splinefun(dose, 100-volume, method="monoH.FC"))
#     meanMHSpl <- if(!inherits(sfun, "try-error")) {
#         tryCatch(pracma::quadgk(function(y) {
#             y*sfun(y, deriv=1)/100 }, 0, max(dose)),
#             error=function(e) return(NA_real_))
#     } else {
#         NA_real_
#     }
#
#     ## locpoly
#     bw <- KernSmooth::dpill(dose, volume)
#     lp <- KernSmooth::locpoly(dose, 100-volume, drv=0, gridsize=10001L,
#                               bandwidth=bw, degree=3)
#     lpfun <- approxfun(lp$x, lp$y, method="linear", rule=2)
#     meanLP <- tryCatch(pracma::quadgk(function(y) {
#         y*numDeriv::grad(lpfun, x=y)/100 }, 0, max(dose)),
#         error=function(e) return(NA_real_))
