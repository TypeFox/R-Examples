## S3 generic function
checkConstraint <-
function(x, constr, byPat=TRUE, semSign=FALSE,
         sortBy=c("none", "observed", "compliance", "structure",
                  "constraint", "patID", "deltaV", "deltaD", "dstMin", "dstMinRel"),
         interp=c("linear", "spline", "smooth"), ...) {
    UseMethod("checkConstraint")
}

checkConstraint.DVHs <-
function(x, constr, byPat=TRUE, semSign=FALSE,
         sortBy=c("none", "observed", "compliance", "structure",
                  "constraint", "patID", "deltaV", "deltaD", "dstMin", "dstMinRel"),
         interp=c("linear", "spline", "smooth"), ...) {
    x <- if(byPat) {
        setNames(list(x), x$structure)
    } else {
        setNames(list(x), x$patID)
    }

    class(x) <- "DVHLst"
    attr(x, which="byPat") <- byPat

    #NextMethod("checkConstraint")
    checkConstraint.DVHLst(x, constr=constr, byPat=byPat, semSign=semSign,
                           sortBy=sortBy, interp=interp, ...)
}

## with byPat=TRUE:  x is a list of DVHs (1 per structure)
## with byPat=FALSE: x is a list of DVHs (1 per ID)
checkConstraint.DVHLst <-
function(x, constr, byPat=TRUE, semSign=FALSE,
         sortBy=c("none", "observed", "compliance", "structure",
                  "constraint", "patID", "deltaV", "deltaD", "dstMin", "dstMinRel"),
         interp=c("linear", "spline", "smooth"), ...) {
    interp <- match.arg(interp)

    ## make sure DVH list is organized as required for byPat
    if(is.null(attributes(x)$byPat) || attributes(x)$byPat != byPat) {
        stop(c("DVH list organization by-patient / by-structure",
               " either could not be determined or is different from byPat"))
    }

    ## restrict x and constr to the same IDs/structures
    xConstrSub  <- harmoConstrDVH(x, constr=constr, byPat=byPat)

    ## determine dose and volume units in x
    xDoseUnits <- vapply(xConstrSub$x, function(y) y$doseUnit,   character(1))
    xVolUnits  <- vapply(xConstrSub$x, function(y) y$volumeUnit, character(1))
    
    ## parse constraint and convert dose/volume units if necessary
    constrParse <- Map(parseConstraint, xConstrSub$constr,
                       doseUnit=xDoseUnits,
                        volUnit=xVolUnits)

    ## calculate difference between 1 observed metric value and 1 constraint
    cmpMetrics <- function(observed, smInv, DV, valCmp, cmp, valRef, dstInfo) {
        valDiff     <- observed - valCmp
        valRelPC    <- 100 * (valDiff/valCmp)
        valDiffInv  <- smInv - valRef
        valRelPCInv <- 100 * (valDiffInv/valRef)
        dstMinRel   <- 100 * dstInfo$dstMinRel

        deltaV   <- if(DV == "V") { valDiff  } else if(DV == "D") { valDiffInv  }
        deltaVpc <- if(DV == "V") { valRelPC } else if(DV == "D") { valRelPCInv }
        deltaD   <- if(DV == "D") { valDiff  } else if(DV == "V") { valDiffInv  }
        deltaDpc <- if(DV == "D") { valRelPC } else if(DV == "V") { valRelPCInv }

        list(observed=observed, deltaV=deltaV, deltaVpc=deltaVpc,
             deltaD=deltaD, deltaDpc=deltaDpc, valCmp=valCmp,
             cmp=cmp, dstMin=dstInfo$dstMin, dstMinRel=dstMinRel,
             ptMinD=dstInfo$ptMinD, ptMinV=dstInfo$ptMinV)
    }

    ## calculate metric and corresponding inverse metrics from given
    ## constraint set - vectorized in metrics
    stridMetrics <- function(dvh, cnstr, interp) {
        observed <- ifelse(cnstr$valid,
                           unlist(getMetric(dvh, metric=cnstr$metric, interp=interp, ...)),
                           NA_real_)

        
        ## inverse metrics
        smInv <- ifelse(cnstr$valid,
                        unlist(getMetric(dvh, metric=cnstr$metricInv, interp=interp, ...)),
                        NA_real_)

        names(observed) <- cnstr$constraint
        names(smInv)    <- cnstr$constraint

        ## distance constraint to closest point on DVH curve
        Dcoord  <- ifelse(cnstr$DV == "D", cnstr$valCmp, cnstr$valRef)
        Vcoord  <- ifelse(cnstr$DV == "V", cnstr$valCmp, cnstr$valRef)
        Dcoord  <- ifelse(cnstr$valid, Dcoord, NA_real_)
        Vcoord  <- ifelse(cnstr$valid, Vcoord, NA_real_)
        volRel  <- ((cnstr$DV == "D") & (cnstr$unitRef == "%")) |
                   ((cnstr$DV == "V") & (cnstr$unitCmp == "%"))
        doseRel <- ((cnstr$DV == "D") & (cnstr$unitCmp == "%")) |
                   ((cnstr$DV == "V") & (cnstr$unitRef == "%"))

        dstDVH <- dvhDistance(dvh,
                              DV=data.frame(D=Dcoord, V=Vcoord,
                                            volRel=volRel, doseRel=doseRel,
                                            unitRef=cnstr$unitRef, unitCmp=cnstr$unitCmp))

        Map(cmpMetrics, observed=observed, smInv=smInv,
            DV=cnstr$DV, valCmp=cnstr$valCmp, cmp=cnstr$cmp, valRef=cnstr$valRef,
            dstInfo=dstDVH)
    }

    ## get requested metrics for each structure/id
    metrics <- Map(stridMetrics, xConstrSub$x, constrParse, interp=interp)

    ## check constraints
    leqGeq <- function(y) {
        observed <- y$observed
        valCmp   <- y$valCmp
        fun      <- y$cmp
        do.call(fun, list(observed, valCmp))
    }

    cmpFun <- function(y) {
        Map(leqGeq, y)
    }

    compL <- Map(cmpFun, metrics)
    metrL <- melt(metrics, value.name="observed")

    ## transform into data frame
    compDF <- melt(compL, value.name="compliance")
    metrDF <- dcast(metrL, L1 + L2 ~ L3, value.var="observed")
    resDF  <- merge(compDF, metrDF)
    names(resDF)[names(resDF) == "L1"] <- if(byPat) {
        "structure"
    } else {
        "patID"
    }

    names(resDF)[names(resDF) == "L2"] <- "constraint"

    ## make sign of deltaD/deltaV semantically indicate compliance
    ## negative -> compliance, positive -> no compliance
    ## compliance may be NA -> catch first to leave sign unchanged
    complianceTF <- resDF$compliance
    complianceTF[is.na(complianceTF)] <- TRUE

    if(semSign) {
        resDF$deltaD    <- ifelse(complianceTF,
                                  -1*resDF$deltaD   *sign(resDF$deltaD),
                                     resDF$deltaD   *sign(resDF$deltaD))
        resDF$deltaDpc  <- ifelse(complianceTF,
                                  -1*resDF$deltaDpc *sign(resDF$deltaDpc),
                                     resDF$deltaDpc *sign(resDF$deltaDpc))
        resDF$deltaV    <- ifelse(complianceTF,
                                  -1*resDF$deltaV   *sign(resDF$deltaV),
                                     resDF$deltaV   *sign(resDF$deltaV))
        resDF$deltaVpc  <- ifelse(complianceTF,
                                  -1*resDF$deltaVpc *sign(resDF$deltaVpc),
                                     resDF$deltaVpc *sign(resDF$deltaVpc))
        resDF$dstMin    <- ifelse(complianceTF,
                                  -1*resDF$dstMin   *sign(resDF$dstMin),
                                     resDF$dstMin   *sign(resDF$dstMin))
        resDF$dstMinRel <- ifelse(complianceTF,
                                  -1*resDF$dstMinRel*sign(resDF$dstMinRel),
                                     resDF$dstMinRel*sign(resDF$dstMinRel))
    }

    ## reorder variables to return
    if(byPat) {
        patID     <- xConstrSub$x[[1]]$patID
        structure <- resDF$structure
    } else {
        patID     <- resDF$patID
        structure <- xConstrSub$x[[1]]$structure
    }

    finDF <- data.frame(patID, structure,
                        constraint=resDF$constraint,
                        observed=resDF$observed,
                        compliance=resDF$compliance,
                        deltaV=resDF$deltaV,
                        deltaVpc=resDF$deltaVpc,
                        deltaD=resDF$deltaD,
                        deltaDpc=resDF$deltaDpc,
                        dstMin=resDF$dstMin,
                        dstMinRel=resDF$dstMinRel,
                        ptMinD=resDF$ptMinD,
                        ptMinV=resDF$ptMinV,
                        stringsAsFactors=FALSE)

    ## sort
    if(!("none" %in% sortBy)) {
        sortBySub <- sortBy[sortBy %in% names(finDF)]
        idx       <- do.call("order", lapply(sortBySub, function(y) with(finDF, get(y))))
        finDF     <- finDF[idx, ]
    }

    return(finDF)
}

## x is a DVH list (1 component per id/structure) of lists (1 DVH for each id/structure)
checkConstraint.DVHLstLst <-
function(x, constr, byPat=TRUE, semSign=FALSE,
         sortBy=c("none", "observed", "compliance", "structure",
                  "constraint", "patID", "deltaV", "deltaD", "dstMin", "dstMinRel"),
         interp=c("linear", "spline", "smooth"), ...) {

    interp <- match.arg(interp)

    ## re-organize x into by-patient or by-structure form if necessary
    isByPat <- attributes(x)$byPat
    xRO <- if(is.null(isByPat) || (isByPat != byPat)) {
        reorgByPat(x, byPat=byPat)
    } else {
        x
    }

    ## make sure DVH and constraint have the same IDs / structures
    xConstrSub <- harmoConstrDVH(xRO, constr=constr, byPat=byPat)

    ## check the constraints
    dots <- list(...)
    args <- list(f=checkConstraint,
                 x=xConstrSub$x,
                 constr=xConstrSub$constr,
                 byPat=byPat,
                 semSign=semSign,
                 interp=interp)
    res <- do.call("Map", c(args, dots))

    ## transform into data frame
    resL    <- melt(res, id.vars=c("patID", "structure", "constraint", "compliance"))
    resL$L1 <- NULL   # redundant with patID (bPat=TRUE) or structure (byPat=FALSE)
    resDF   <- dcast(resL, patID + structure + constraint + compliance ~ variable,
                     value.var="value")

    ## reorder variables to return
    finDF <- data.frame(patID=resDF$patID,
                        structure=resDF$structure,
                        constraint=resDF$constraint,
                        observed=resDF$observed,
                        compliance=resDF$compliance,
                        deltaV=resDF$deltaV,
                        deltaVpc=resDF$deltaVpc,
                        deltaD=resDF$deltaD,
                        deltaDpc=resDF$deltaDpc,
                        dstMin=resDF$dstMin,
                        dstMinRel=resDF$dstMinRel,
                        ptMinD=resDF$ptMinD,
                        ptMinV=resDF$ptMinV,
                        stringsAsFactors=FALSE)

    ## sort
    if(!("none" %in% sortBy)) {
        sortBySub <- sortBy[sortBy %in% names(finDF)]
        idx <- do.call("order", lapply(sortBySub, function(y) with(finDF, get(y))))
        finDF <- finDF[idx, ]
    }

    return(finDF)
}
