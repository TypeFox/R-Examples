## S3 generic method
getMeanDVH <-
function(x, fun=list(mean=mean, median=median, sd=sd),
         cumul=TRUE, thin=1, byPat=TRUE, patID=NULL, structure=NULL,
         fixed=TRUE) {
    UseMethod("getMeanDVH")
}

## for completeness sake - "mean" of just 1 DVH
getMeanDVH.DVHs <-
function(x, fun=list(mean=mean, median=median, sd=sd),
         cumul=TRUE, thin=1, byPat=TRUE, patID=NULL, structure=NULL,
         fixed=TRUE) {

    x <- if(byPat) {
        setNames(list(x), x$structure)
    } else {
        setNames(list(x), x$patID)
    }

    class(x) <- "DVHLst"
    attr(x, which="byPat") <- byPat

    getMeanDVH.DVHLst(x, fun=fun, cumul=cumul, thin, byPat=byPat,
                      patID=patID, structure=structure, fixed=fixed)
}

getMeanDVH.DVHLst <-
function(x, fun=list(mean=mean, median=median, sd=sd),
         cumul=TRUE, thin=1, byPat=TRUE, patID=NULL, structure=NULL,
         fixed=TRUE) {

    ## make sure DVH list is organized as required for byPat
    if(is.null(attributes(x)$byPat) || attributes(x)$byPat != byPat) {
        stop(c("DVH list organization by-patient / by-structure ",
               "either could not be determined or is different from byPat"))
    }

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

    ## get dose range and number of dose nodes
    rangeD <- c(0, max(vapply(x, function(z) {    max(z$dvh[ , "dose"]) }, numeric(1))))
    nodes  <-      max(vapply(x, function(z) { length(z$dvh[ , "dose"]) }, numeric(1)))
        
    ## coarser dose grid for M+SD but with at least 100 nodes
    nodes <- max(100, ceiling(nodes/thin))

    ## average cumulative or differential DVH?
    x <- if(cumul) {
        ## cumulative with linear interpolation
        convertDVH.DVHLst(x, toType="asis", interp="linear",
                          nodes=nodes, rangeD=rangeD, perDose=TRUE)
    } else {
        ## differential with linear interpolation
        convertDVH.DVHLst(x, toType="differential", interp="linear",
                          nodes=nodes, rangeD=rangeD, perDose=TRUE)
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
    
    ## combine list to data frame
    dvhDF <- do.call("rbind", dvhDFL)

    ## generate point-wise mean/sd for dose -> aggregate over dose
    getAggr <- function(fun, symbol) {
        dat <- if(byPat) {
            aggregate(cbind(volume, volumeRel) ~ patID + dose,
                      data=dvhDF, FUN=fun, na.action=na.pass)
        } else {
            aggregate(cbind(volume, volumeRel) ~ structure + dose,
                      data=dvhDF, FUN=fun, na.action=na.pass)
        }

        ## rename columns in aggregated data frames
        namesDat <- names(dat)
        namesDat[namesDat == "volume"]    <- paste0("volume",    symbol)
        namesDat[namesDat == "volumeRel"] <- paste0("volumeRel", symbol)
        names(dat) <- namesDat
        dat
    }

    ## get all point-wise estimates
    dfL <- Map(getAggr, fun, toupper(names(fun)))

    ## combine point-wise estimates
    dfMSD <- Reduce(merge, dfL)
    rownames(dfMSD) <- NULL
    dfMSD
}

## x is a DVH list (1 per id or 1 per structure) of lists
## plots many DVH files
## either for many patients   -> multiple structures per DVH
## or     for many structures -> multiple patients   per DVH
getMeanDVH.DVHLstLst <-
function(x, fun=list(mean=mean, median=median, sd=sd),
         cumul=TRUE, thin=1, byPat=TRUE, patID=NULL, structure=NULL,
         fixed=TRUE) {

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

    resDFL <- Map(getMeanDVH, x, fun=list(fun), cumul=cumul, thin,
                  byPat=byPat, patID=list(patID), structure=list(structure),
                  fixed=fixed)

    resDF <- do.call("rbind", resDFL)
    rownames(resDF) <- NULL
    resDF
}

