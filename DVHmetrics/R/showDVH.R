## combine DVH data from a DVHLst object
## interp is FALSE or gives number of nodes to interpolate
combineDVHs <-
function(x, cumul=TRUE, interp=FALSE, rangeD=NULL) {
    ## differential DVH or interpolation - convert first
    x <- if(cumul && !interp) {
        ## cumulative without interpolation -> nothing to do
        x
    } else if(cumul && interp) {
        ## cumulative with linear interpolation
        convertDVH.DVHLst(x, toType="asis", interp="linear",
                          nodes=interp, rangeD=rangeD, perDose=TRUE)
    } else if(!cumul && !interp) {
        ## differential without interpolation
        convertDVH.DVHLst(x, toType="differential", perDose=TRUE)
    } else if(!cumul && interp) {
        ## differential with linear interpolation
        convertDVH.DVHLst(x, toType="differential", interp="linear",
                          nodes=interp, rangeD=rangeD, perDose=TRUE)
    }

    ## extract actual DVH from each DVHs object
    dvhDFL <- if(cumul) {
        ## cumulative DVH
        lapply(x, function(y) {
            data.frame(y$dvh, patID=y$patID, structure=y$structure,
                       stringsAsFactors=FALSE)
        })
    } else {
        ## differential DVH
        lapply(x, function(y) {
            data.frame(y$dvhDiff, patID=y$patID, structure=y$structure,
                       stringsAsFactors=FALSE)
        })
    }
    
    do.call("rbind", dvhDFL)
}

## S3 generic method
showDVH <-
function(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL, rel=TRUE,
         guessX=TRUE, thresh=1, addMSD=FALSE, show=TRUE, fixed=TRUE) {
    UseMethod("showDVH")
}

## plots 1 DVH file for 1 id and 1 structure
showDVH.DVHs <-
function(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL, rel=TRUE,
         guessX=TRUE, thresh=1, addMSD=FALSE, show=TRUE, fixed=TRUE) {
    x <- if(byPat) {
        setNames(list(x), x$structure)
    } else {
        setNames(list(x), x$patID)
    }

    class(x) <- "DVHLst"
    attr(x, which="byPat") <- byPat

    showDVH.DVHLst(x, cumul=cumul, byPat=byPat, patID=patID, structure=structure,
                   rel=rel, guessX=guessX, thresh=thresh, addMSD=addMSD, show=show)
}

## plots 1 list of DVH objects
## for byPat=TRUE:  1 patient   -> multiple structures
## for byPat=FALSE: 1 structure -> multiple patients
showDVH.DVHLst <-
function(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL, rel=TRUE,
         guessX=TRUE, thresh=1, addMSD=FALSE, show=TRUE, fixed=TRUE) {

    ## make sure DVH list is organized as required for byPat
    if(is.null(attributes(x)$byPat) || attributes(x)$byPat != byPat) {
        stop(c("DVH list organization by-patient / by-structure ",
               "either could not be determined or is different from byPat"))
    }

    guessX <- as.numeric(guessX)

    ## if patIDs are selected, filter them here -> strips DVHLst class
    if(!is.null(patID)) {
        p <- trimWS(patID, side="both")
        x <- Filter(function(y) { 
            if(fixed) {
                any(y$patID %in% p)
            } else {
                any(grepl(paste(p, collapse="|"), y$patID))
            }
        }, x)

        if(length(x) < 1L) { stop("No selected patient found") }
    }

    ## if structures are selected, filter them here
    if(!is.null(structure)) {
        s <- trimWS(structure, side="both")
        x <- Filter(function(y) {
            if(fixed) {
                any(y$structure %in% s)
            } else {
                any(grepl(paste(s, collapse="|"), y$structure))
            }
        }, x)

        if(length(x) < 1L) { stop("No selected structure found") }
    }

    ## combine all data frames in the DVHs
    ## find number of dose values in DVHs
    if(addMSD) {
        rangeD  <- c(0, max(vapply(x, function(z) {    max(z$dvh[ , "dose"]) }, numeric(1))))
        doseLen <-      max(vapply(x, function(z) { length(z$dvh[ , "dose"]) }, numeric(1)))
        
        ## coarser dose grid for M+SD but with at least 100 nodes
        doseLen <- max(100, ceiling(doseLen/3))
    } else {
        rangeD  <- NULL
        doseLen <- FALSE
    }

    dvhDF <- combineDVHs(x, cumul=cumul, interp=doseLen, rangeD=rangeD)

    ## choose upper x-axis limit (dose) - add 10%
    ## TODO: up to first dose for which all structures reach threshold
    xMax <- if((guessX == 1L) && !all(is.na(x[[1]]$dvh[ , "volumeRel"]))) {
        volGEQ <- lapply(x, function(y) {
            y$dvh[ , "dose"][y$dvh[ , "volumeRel"] >= thresh] })
        1.1*max(unlist(volGEQ), na.rm=TRUE)
    } else {
        1.1*max(c(guessX, dvhDF$dose))
    }
    
    ## set plot volume to absolute or relative
    ## and choose upper y-axis limit (volume) - add 3%
    if(rel) {
        ## relative volume - check if available
        if(all(is.na(dvhDF$volumeRel))) {
            warning("All relative volumes are missing, will try to show absolute volume")
            yMax <- 1.03*max(dvhDF$volume, na.rm=TRUE)
            rel  <- FALSE
            dvhDF$volPlot <- dvhDF$volume
        } else {
            yMax <- min(c(103, 1.03*max(dvhDF$volumeRel, na.rm=TRUE)))
            dvhDF$volPlot <- dvhDF$volumeRel
        }
    } else {
        ## absolute volume - check if available
        if(all(is.na(dvhDF$volume))) {
            warning("All absolute volumes are missing, will try to show relative volume")
            yMax <- min(c(103, 1.03*max(dvhDF$volumeRel, na.rm=TRUE)))
            rel  <- TRUE
            dvhDF$volPlot <- dvhDF$volumeRel
        } else {
            yMax <- 1.03*max(dvhDF$volume, na.rm=TRUE)
            dvhDF$volPlot <- dvhDF$volume
        }
    }

    ## check if absolute dose is available
    isDoseRel <- if(all(is.na(dvhDF$dose))) {
        warning("All absolute doses are missing, will try to show relative dose")
        dvhDF$dose <- dvhDF$doseRel
        TRUE
    } else {
        FALSE
    }

    ## check if point-wise volume mean and sd should be plotted
    if(addMSD) {
        ## generate point-wise mean/sd for dose -> aggregate over dose
        if(byPat) {
            dfM  <- aggregate(volPlot ~ patID + dose,     data=dvhDF, FUN=mean, na.rm=TRUE)
            dfSD <- aggregate(volPlot ~ patID + dose,     data=dvhDF, FUN=sd,   na.rm=TRUE)
        } else {
            dfM  <- aggregate(volPlot ~ structure + dose, data=dvhDF, FUN=mean, na.rm=TRUE)
            dfSD <- aggregate(volPlot ~ structure + dose, data=dvhDF, FUN=sd,   na.rm=TRUE)
        }

        ## rename columns in aggregated data frames
        ## dfM stays as is because "volPlot" is used in top layer of the plot
        namesDFSD <- names(dfSD)
        namesDFSD[namesDFSD == "volPlot"] <- "volSD"
        names(dfSD) <- namesDFSD

        dfMSD <- Reduce(merge, list(dfM, dfSD))
        dfMSD$volLo1SD <- dfMSD$volPlot -   dfMSD$volSD
        dfMSD$volHi1SD <- dfMSD$volPlot +   dfMSD$volSD
        dfMSD$volLo2SD <- dfMSD$volPlot - 2*dfMSD$volSD
        dfMSD$volHi2SD <- dfMSD$volPlot + 2*dfMSD$volSD
    }

    ## title string
    strTitle <- if(byPat) {
        paste0("patient ",   x[[1]]$patID)
    } else {
        paste0("structure ", x[[1]]$structure)
    }

    doseUnit <- if(isDoseRel) {
        "%"
    } else {
        x[[1]]$doseUnit
    }

    volUnit <- if(rel) {
        if(cumul) {
            "%"
        } else {
            paste0("%/", doseUnit)
        }
    } else {
        if(cumul) {
            x[[1]]$volumeUnit
        } else {
            paste0(x[[1]]$volumeUnit, "/", doseUnit)
        }
    }

    ## do we have many structure/ID categories? -> more legend columns
    nCateg <- if(byPat) {
        length(unique(dvhDF$structure))
    } else {
        length(unique(dvhDF$patID))
    }

    ## 15 entries per legend column
    nLegendCols <- (nCateg %/% 15) + 1            # number of categories

    ## do the actual plotting
    diag <- ggplot(dvhDF, aes_string(x="dose", y="volPlot"))

    ## rescale x-axis
    diag <- if(is.finite(xMax)) {
        diag + coord_cartesian(xlim=c(0, xMax))
    } else {
        diag
    }
    
    ## rescale y-axis
    diag <- if(is.finite(yMax)) {
        diag + coord_cartesian(ylim=c(0, yMax))
    } else {
        diag
    }

    ## point-wise 1-SD, 2-SD shaded areas?
    diag <- if(addMSD) {
        diag +
        geom_ribbon(data=dfMSD,
                    aes_string(x="dose", ymin="volLo1SD", ymax="volHi1SD"),
                    alpha=0.25, linetype="blank") +
        geom_ribbon(data=dfMSD,
                    aes_string(x="dose", ymin="volLo2SD", ymax="volHi2SD"),
                    alpha=0.2, linetype="blank")
    } else {
        diag
    }
    
    ## actual DVHs
    diag <- if(byPat) {
        diag + geom_line(data=dvhDF,
                         aes_string(x="dose", y="volPlot", color="structure"),
                         size=1.2)
    } else {
        diag + geom_line(data=dvhDF,
                         aes_string(x="dose", y="volPlot", color="patID"),
                         size=1.2)
    }

    ## add point-wise mean DVH?
    diag <- if(addMSD) {
        diag + geom_line(data=dfMSD,
                         aes_string(x="dose", y="volPlot"),
                         color="black", size=1.2)
    } else {
        diag
    }

    ## final diagram
    diag <- diag +
        ggtitle(strTitle) +
        theme_bw() +
        expand_limits(y=0) +                      # make sure 0 is included
        scale_y_continuous(expand=c(0, 0.6)) +    # no space below y=0
        guides(color=guide_legend(ncol=nLegendCols)) +
        xlab(paste0("Dose [",   doseUnit, "]")) +
        ylab(paste0("Volume [", volUnit,  "]"))

    if(show) {
        print(diag)
    }

    return(invisible(diag))
}

## x is a DVH list (1 per id or 1 per structure) of lists
## plots many DVH files
## either for many patients   -> multiple structures per DVH
## or     for many structures -> multiple patients   per DVH
showDVH.DVHLstLst <-
function(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL, rel=TRUE,
         guessX=TRUE, thresh=1, addMSD=FALSE, show=TRUE, fixed=TRUE) {

    ## re-organize x into by-patient or by-structure form if necessary
    isByPat <- attributes(x)$byPat
    if(is.null(isByPat) || (isByPat != byPat)) {
        x <- reorgByPat(x, byPat=byPat)
    }

    ## if byPat=TRUE and patIDs are selected, filter them here
    if(byPat && !is.null(patID)) {
        p <- trimWS(patID, side="both")
        x <- Filter(function(y) {
            if(fixed) {
                any(y[[1]]$patID %in% p)
            } else {
                any(grepl(paste(p, collapse="|"), y[[1]]$patID))
            }
        }, x)

        if(length(x) < 1L) { stop("No selected patient found") }
        patID <- NULL
    }

    ## if byPat=FALSE and structures are selected, filter them here
    if(!byPat && !is.null(structure)) {
        s <- trimWS(structure, side="both")
        x <- Filter(function(y) {
            if(fixed) {
                any(y[[1]]$structure %in% s)
            } else {
                any(grepl(paste(s, collapse="|"), y[[1]]$structure))
            }
        }, x)

        if(length(x) < 1L) { stop("No selected structure found") }
        structure <- NULL
    }

    diagL <- Map(showDVH, x, cumul=cumul, byPat=byPat, rel=rel,
                 patID=list(patID), structure=list(structure),
                 guessX=guessX, thresh=thresh, addMSD=addMSD, show=show, fixed=fixed)

    return(invisible(diagL))
}
