## S3 generic function
showConstraint <-
function(x, constr, byPat=TRUE, rel=TRUE, guessX=TRUE, thresh=1) {
    UseMethod("showConstraint")
}

showConstraint.DVHs <-
function(x, constr, byPat=TRUE, rel=TRUE, guessX=TRUE, thresh=1) {
    x <- if(byPat) {
        setNames(list(x), x$structure)
    } else {
        setNames(list(x), x$patID)
    }

    class(x) <- "DVHLst"
    attr(x, which="byPat") <- byPat

    #NextMethod("showConstraint")
    showConstraint.DVHLst(x, constr=constr, byPat=byPat,
                           guessX=guessX, rel=rel, thresh=thresh)
}

## for byPat=TRUE
## x is a list of DVHs (1 per structure)
## for byPat=FALSE
## x is a list of DVHs (1 per patient)
showConstraint.DVHLst <-
function(x, constr, byPat=TRUE, rel=TRUE, guessX=TRUE, thresh=1) {
    ## make sure DVH list is organized as required for byPat
    if(is.null(attributes(x)$byPat) || attributes(x)$byPat != byPat) {
        stop(c("DVH list organization by-patient / by-structure",
               " either could not be determined or is different from byPat"))
    }

    ## restrict x and constr to the same IDs/structures
    xConstrSub <- harmoConstrDVH(x, constr=constr, byPat=byPat)
    
    ## determine dose and volume units in x
    xDoseUnits <- vapply(xConstrSub$x, function(y) y$doseUnit,   character(1))
    xVolUnits  <- vapply(xConstrSub$x, function(y) y$volumeUnit, character(1))
    
    ## parse constraint and convert dose/volume units if necessary
    constrParse <- Map(parseConstraint, xConstrSub$constr,
                       doseUnit=xDoseUnits,
                        volUnit=xVolUnits)

    ## data frame from parsed constraint list
    cDFL <- lapply(constrParse, function(y) data.frame(y, stringsAsFactors=FALSE))
    cDF  <- melt(cDFL, id.vars=c("constraint", "valid", "metric", "DV",
                                 "valRef", "unitRef", "cmp", "unitCmp", "metricInv"),
                 value.name="valCmp")
    cDF$variable <- NULL
    names(cDF)[names(cDF) == "L1"] <- if(byPat) {
        "structure"
    } else {
        "patID"
    }

    ## distance constraint to closest point on DVH curve
    getMinDstPt <- function(dvh, cnstr) {
        Dcoord  <- ifelse(cnstr$DV == "D", cnstr$valCmp, cnstr$valRef)
        Vcoord  <- ifelse(cnstr$DV == "V", cnstr$valCmp, cnstr$valRef)
        Dcoord  <- ifelse(cnstr$valid, Dcoord, NA_real_)
        Vcoord  <- ifelse(cnstr$valid, Vcoord, NA_real_)
        volRel  <- ((cnstr$DV == "D") & (cnstr$unitRef == "%")) |
                   ((cnstr$DV == "V") & (cnstr$unitCmp == "%"))
        doseRel <- ((cnstr$DV == "D") & (cnstr$unitCmp == "%")) |
                   ((cnstr$DV == "V") & (cnstr$unitRef == "%"))
        
        dstL <- dvhDistance(dvh, 
                            DV=data.frame(D=Dcoord, V=Vcoord,
                                          volRel=volRel, doseRel=doseRel,
                                          unitRef=cnstr$unitRef, unitCmp=cnstr$unitCmp))
        setNames(dstL, cnstr$constraint)
    }

    minDst   <- Map(getMinDstPt, xConstrSub$x, constrParse)
    minDstL  <- melt(minDst)
    minDstDF <- dcast(minDstL, L1 + L2 ~ L3, value.var="value")
    names(minDstDF)[names(minDstDF) == "L1"] <- if(byPat) {
        "structure"
    } else {
        "patID"
    }

    names(minDstDF)[names(minDstDF) == "L2"] <- "constraint"
    
    ## merge distance to DVH with remaining data
    allDF <- merge(cDF, minDstDF)

    ## from given info, calculate dose, volume, relative volume
    getDoseVol <- function(DV, valRef, unitRef, cmp, unitCmp, valCmp, DVH,
                           ptMinD, ptMinV, y) {
        if(!(DVH %in% names(y))) {
            warning(c("DVH for ", DVH, " not found"))
            return(NA)
        }

        if(is.na(ptMinD) || is.na(ptMinV) || is.na(valRef)) {
            doseAbs   <- NA_real_
            volAbs    <- NA_real_
            volRel    <- NA_real_
            ptMinDabs <- NA_real_
            ptMinVabs <- NA_real_
            ptMinVrel <- NA_real_
        } else if(DV == "V") {
            ## input  is absolute or relative dose
            ## output is absolute or relative volume
            if(unitRef == "%") {
                ## input is relative dose -> convert rel to abs
                ## -> we need doseRx
                stopifnot(    !is.null(y[[DVH]]$doseRx),
                          all(!is.na(  y[[DVH]]$doseRx)))
                doseAbs   <- (valRef/100) * y[[DVH]]$doseRx
                ptMinDabs <- (ptMinD/100) * y[[DVH]]$doseRx
            } else if(unitRef %in% c("GY", "CGY")) {
                ## input is absolute dose -> nothing to do
                doseAbs   <- valRef
                ptMinDabs <- ptMinD
            }

            if(unitCmp == "%") {
                ## output is relative volume -> convert rel to abs
                volAbs    <- (valCmp/100) * y[[DVH]]$structVol
                volRel    <- valCmp
                ptMinVabs <- (ptMinV/100) * y[[DVH]]$structVol
                ptMinVrel <- ptMinV
            } else if(unitCmp == "CC") {
                ## output is absolute volume -> convert abs to rel
                volAbs    <- valCmp
                volRel    <- 100*valCmp / y[[DVH]]$structVol[1]
                ptMinVabs <- ptMinV
                ptMinVrel <- 100*ptMinV / y[[DVH]]$structVol[1]
            }
        } else if(DV == "D") {
            ## input  is absolute or relative volume
            ## output is absolute or relative dose
            if(unitRef == "%") {
                ## input is relative volume -> convert rel to abs
                volAbs    <- (valRef/100) * y[[DVH]]$structVol
                volRel    <- valRef
                ptMinVabs <- (ptMinV/100) * y[[DVH]]$structVol
                ptMinVrel <- ptMinV
            } else if(unitRef == "CC") {
                ## input is absolute volume -> convert abs to rel
                volAbs    <- valRef
                volRel    <- 100*valRef / y[[DVH]]$structVol
                ptMinVabs <- ptMinV
                ptMinVrel <- 100*ptMinV / y[[DVH]]$structVol
            }

            if(unitCmp == "%") {
                ## output is relative dose -> convert rel to abs
                ## -> we need doseRx
                stopifnot(    !is.null(y[[DVH]]$doseRx),
                          all(!is.na(  y[[DVH]]$doseRx)))
                doseAbs   <- (valCmp/100) * y[[DVH]]$doseRx
                ptMinDabs <- (ptMinD/100) * y[[DVH]]$doseRx
            } else if(unitCmp %in% c("GY", "CGY")) {
                ## output is absolute dose -> nothing to do
                doseAbs   <- valCmp
                ptMinDabs <- ptMinD
            }
        }

        return(data.frame(doseAbs=doseAbs, volAbs=volAbs, volRel=volRel,
                          ptMinDabs=ptMinDabs,
                          ptMinVabs=ptMinVabs, ptMinVrel=ptMinVrel))
    }

    ## get dose and (abs/rel) volume coords
    doseVolL <- Map(getDoseVol,
                    DV=allDF$DV, valRef=allDF$valRef, unitRef=allDF$unitRef,
                    cmp=allDF$cmp, unitCmp=allDF$unitCmp, valCmp=allDF$valCmp,
                    DVH=if(byPat) { allDF$structure } else { allDF$patID },
                    ptMinD=allDF$ptMinD, ptMinV=allDF$ptMinV,
                    y=list(xConstrSub$x))

    doseVolDF <- do.call("rbind", doseVolL)

    ## which volume should be plotted?
    if(rel) {
        doseVolDF$volPlot    <- doseVolDF$volRel
        doseVolDF$ptMinVplot <- doseVolDF$ptMinVrel
    } else {
        doseVolDF$volPlot    <- doseVolDF$volAbs
        doseVolDF$ptMinVplot <- doseVolDF$ptMinVabs
    }

    ## add to existing data and remove invalid constraints
    diagDF <- cbind(allDF, doseVolDF)
    diagDF <- diagDF[diagDF$valid, ]

    ## create separate arrow in absolute size - point lower left and upper right
    grLen <- 0.7
    arrLL <- segmentsGrob(x0=grLen, y0=grLen, x1=0, y1=0, default.units="cm",
                arrow=arrow(angle=45, length=unit(grLen, "cm")),
                gp=gpar(lwd=3, lineend="butt", linejoin="mitre"))

    arrUR <- segmentsGrob(x0=-grLen, y0=-grLen, x1=0, y1=0, default.units="cm",
                arrow=arrow(angle=45, length=unit(grLen, "cm")),
                gp=gpar(lwd=3, lineend="butt", linejoin="mitre"))

    constraintArr <- function(dat) {
        getArr <- function(x, y, cmp) {
            if(grepl("^<.*", cmp)) {
                annotation_custom(grobTree(arrLL), xmin=x, ymin=y, xmax=x+1, ymax=y+1)
            } else if(grepl("^>.*", cmp)) {
                annotation_custom(grobTree(arrUR), xmin=x, ymin=y, xmax=x+1, ymax=y+1)
            }
        }

        annL <- Map(getArr, dat$doseAbs, dat$volPlot, dat$cmp)
    }

    ## do the actual plotting
    ## get DVH plot for relevant structures / IDs
    ## constraint may be more extreme than actual doses
    doseMax <- max(vapply(x, function(y) { max(y$dvh[ , "dose"]) }, numeric(1)))
    if(isTRUE(any(diagDF$doseAbs > doseMax))) {
        guessX <- max(diagDF$doseAbs)
    }

    diag <- showDVH(xConstrSub$x,
                    cumul=TRUE, byPat=byPat, rel=rel, guessX=guessX,
                    thresh=thresh, show=FALSE)

    ## add constraint arrow to plot
    diagC <- if(byPat) {
        diag + constraintArr(diagDF) +
            geom_point(data=diagDF,
                       aes_string(x="doseAbs", y="volPlot",
                                  color="structure", shape="constraint"),
                       size=5) +
            geom_point(data=diagDF,
                       aes_string(x="ptMinDabs", y="ptMinVplot",
                                  color="structure", shape="constraint"),
                       size=5)
    } else {
        diag + constraintArr(diagDF) +
            geom_point(data=diagDF,
                       aes_string(x="doseAbs",   y="volPlot",
                                  color="patID", shape="constraint"),
                       size=5) +
            geom_point(data=diagDF,
                       aes_string(x="ptMinDabs", y="ptMinVplot",
                                  color="patID", shape="constraint"),
                       size=5)
    }

    ## for more than 6 constraints, we need more shapes
    nConstr <- nlevels(as.factor(diagDF$constraint))
    if(nConstr > 6) {
        diagC <- diagC + scale_shape_manual(values=seq_len(nConstr))
    }

    print(diagC)

    return(invisible(diagC))
}

## x is a DVH list (1 component per id) of lists
showConstraint.DVHLstLst <-
function(x, constr, byPat=TRUE, rel=TRUE, guessX=TRUE, thresh=1) {
    ## re-organize x into by-patient or by-structure form if necessary
    isByPat <- attributes(x)$byPat
    xRO <- if(is.null(isByPat) || (isByPat != byPat)) {
        reorgByPat(x, byPat=byPat)
    } else {
        x
    }

    ## make sure DVH and constraint have the same IDs / structures
    xConstrSub <- harmoConstrDVH(xRO, constr=constr, byPat=byPat)

    diagL <- Map(showConstraint,
                 x=xConstrSub$x, constr=xConstrSub$constr,
                 byPat=byPat, rel=rel, guessX=guessX, thresh=thresh)

    return(invisible(diagL))
}
