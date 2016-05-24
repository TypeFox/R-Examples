## parse character vector from Cadplan DVH file
parseCadplan <- function(x, planInfo=FALSE, courseAsID=FALSE) {
    ## function to extract one information element from a number of lines
    getElem <- function(pattern, ll, trim=TRUE, iCase=FALSE, collWS=TRUE) {
        line <- ll[grep(pattern, ll)]
        elem <- sub("^.+?:[[:blank:]]*([[:alnum:][:punct:][:blank:]]+$)", "\\1",
                    line, ignore.case=iCase, perl=TRUE)
        elem <- if(trim) {
            trimWS(elem, side="both")
        } else {
            elem
        }

        if(collWS) {
            collWS(elem)
        } else {
            elem
        }
    }

    getPlan <- function(pattern, ll, trim=TRUE, iCase=FALSE, collWS=TRUE) {
        line <- ll[grep(pattern, ll)]
        elem <- sub("^[[:alpha:]]+[[:blank:]]+([[:alnum:]]+[[:alnum:][:blank:]]*)", "\\1",
                    line, ignore.case=iCase, perl=TRUE)
        elem <- if(trim) {
            trimWS(elem, side="both")
        } else {
            elem
        }

        if(collWS) {
            collWS(elem)
        } else {
            elem
        }
    }

    getDose <- function(pattern, ll, doseRx, percent=TRUE) {
        line <- ll[grep(pattern, ll)]
        elem <- sub("^.+?:[[:blank:]]+([[:alnum:][:punct:]]+[[:blank:]]*$)", "\\1", line, perl=TRUE)
        num  <- trimWS(elem)
        if(percent && any(grepl("%", line))) {
            if(!missing(doseRx)) {
                doseRx * as.numeric(num)/100
            } else {
                NA_real_
            }
        } else {
            as.numeric(num)
        }
    }

    getDoseUnit <- function(ll) {
        line <- ll[grep("^Prescr\\. dose.+:", ll)]
        elem <- sub("^.+\\((GY|CGY)\\)[[:blank:]]*:.+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimWS(elem))
    }

    getVolUnit <- function(ll) {
        line <- ll[grep("^Volume.+:", ll)]
        elem <- sub("^.+\\((CC)\\)[[:blank:]]*:.+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimWS(elem))
    }

    getDVHtype <- function(ll) {
        line <- ll[grep("^(Cumulative|Differential) (Dose Volume.+|DVH)", ll)]
        elem <- sub("^(Cumulative|Differential).+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        tolower(trimWS(elem))
    }

    ## split file into list of structure sections
    sStart <- grep("^Histogram[[:blank:]]+: [[:alnum:][:punct:]]+", x)  # start of sections
    sLen   <- diff(c(sStart, length(x)+1))        # length of sections
    if((length(sLen) < 1L) || all(sLen < 1L)) {
        stop("No structures found")
    }

    structList <- Map(function(start, len) x[start:(start+len-1)], sStart, sLen)

    ## extract file header and header info
    header  <- x[seq_len(sStart[1]-1)]                       # header
    patName <- getElem("Patient Name[[:blank:]]*:", header)  # patient name
    patID   <- getElem("^Patient ID[[:blank:]]*:",  header)  # patient id
    plan    <- getPlan("^PLAN", header, iCase=TRUE, collWS=TRUE)  # treatment plan
    DVHdate <- getElem("^Date[[:blank:]]+:", header)         # export date
    DVHtype <- getDVHtype(header)

    ## prescribed dose may be encoded in plan
    if(tolower(planInfo) == "doserx") {
        doseRxUnit <- toupper(sub("^.+_([.[:digit:]]+)(GY|CGY)_.*", "\\2",
                                  plan, perl=TRUE, ignore.case=TRUE))
    
        if(!grepl("^(GY|CGY)$", doseRxUnit)) {
            warning("Could not determine dose Rx unit")
            doseRxUnit <- NA_character_
        }
        
        drx <- sub("^[[:alnum:]]+_([.[:digit:]]+)(GY|CGY)_[[:alnum:]]*", "\\1",
                   plan, perl=TRUE, ignore.case=TRUE)
        doseRx <- as.numeric(drx)
    } else {
        doseRx     <- NA_real_
        doseRxUnit <- NA_character_
    }

    ## extract DVH from one structure section and store in a list
    ## with DVH itself as a matrix
    getDVH <- function(strct, info) {
        plan <- info$plan

        ## extract structure, prescribed dose, volume, completed to 100% dose Rx,
        ## dose min, max, mean, median, mode, and sd
        structure <- getElem("^Histogram.*:", strct)

        isoDoseRx0 <- getElem("^% for dose[[:blank:]]*:", strct)
        ## check if sum plan
        isoDoseRx <- if((length(isoDoseRx0) > 0) && (isoDoseRx0 != "not defined")) {
            as.numeric(isoDoseRx0)
        } else {                                        # sum plan -> use plan info?
            if(tolower(planInfo) == "doserx") {
                warning("Iso-dose-Rx is assumed to be 100")
                100
            } else {
                warning("No info on % for dose")
                NA_real_
            }
        }

        doseUnit <- getDoseUnit(strct)
        if(!grepl("^(GY|CGY)$", doseUnit)) {
            doseUnit <- NA_character_
            warning("Could not determine dose measurement unit")
        }

        doseRx0 <- getElem("^Prescr\\. dose.*:", strct)
        ## check if sum plan
        doseRx <- if((length(doseRx0) > 0) && (doseRx0 != "not defined")) {
            doseRxUnit <- doseUnit
            getDose("^Prescr\\. dose.*:", strct)
        } else {                                        # sum plan
            ## doseRx may be encoded in plan name
            if(tolower(planInfo) == "doserx") {
                doseRxUnit <- info$doseRxUnit
                info$doseRx
            } else {
                warning("No info on prescribed dose")
                doseRxUnit <- NA_character_
                NA_real_
            }
        }

        ## check if we have dose Rx
        ## if so, does it have the same unit as doseUnit -> convert
        if(!is.na(doseUnit) && !is.na(doseRxUnit)) {
            if((doseUnit == "GY") && (doseRxUnit == "CGY")) {
                doseRx <- doseRx/100
            } else if((doseUnit == "CGY") && (doseRxUnit == "GY")) {
                doseRx <- doseRx*100
            }
        }

        volumeUnit <- getVolUnit(strct)
        volumeUnit <- if(grepl("^CC$", volumeUnit)) {
            "CC"
        } else if(grepl("^%", volumeUnit)) {
            "PERCENT"
        } else {
            warning("Could not determine volume measurement unit")
            NA_character_
        }

        structVol <- as.numeric(getElem("^Volume.*:", strct))
        doseMin   <- getDose("^Dose minimum.*:", strct, doseRx)
        doseMax   <- getDose("^Dose maximum.*:", strct, doseRx)
        doseAvg   <- getDose("^Dose mean.*:",    strct, doseRx)
        doseMed   <- getDose("^Dose median.*:",  strct, doseRx)
        doseMod   <- getDose("^Dose modal.*:",   strct, doseRx)
        doseSD    <- getDose("^Standard dev.*:", strct, doseRx, percent=FALSE)

        ## find DVH
        ## DVH column headers
        colHead  <- grep("DOSE[[:blank:]]*\\((%|GY|CGY)\\).+VOLUME", strct,
                         ignore.case=TRUE, perl=TRUE)
        dvhStart <- colHead+1            # first numeric line of DVH
        dvhLen   <- length(strct) - dvhStart + 1
        if((length(dvhLen) < 1L) || dvhLen < 1L) {
            stop("No DVH data found")
        }

        ## column headers
        vars1 <- unlist(strsplit(strct[colHead],
                        split="\\([[:alpha:]%]+\\)", fixed=FALSE, perl=TRUE))
        ## remove leading and trailing white space
        vars2 <- tolower(trimWS(vars1))

        ## make sure we recognize all columns in the DVH
        patDose    <- "^dose"
        patDoseRel <- "^relative dose"
        patVol     <- "^volume.+cm3"
        hits <- sum(c(grepl(patDose, vars2), grepl(patDoseRel, vars2), grepl(patVol, vars2)))
        if(length(vars2) != hits) {
            stop(c("Could not identify all DVH columns"),
        		 paste(vars2, collapse=", "))
        }

        ## replace column headers
        vars3 <- vars2
        vars3[grep(patDose,    vars2)] <- "dose"
        vars3[grep(patDoseRel, vars2)] <- "doseRel"
        vars3[grep(patVol,     vars2)] <- "volume"

        ## extract DVH as a matrix and store preceding information
        ## read line length(strct) for cases where file does not end with a
        ## blank line -> this will then be last DVH line, otherwise blank
        ## check if dvh is all blank -> no data
        if(all(!nzchar(strct[dvhStart:length(strct)]))) {
            return(NULL)
        }

        dvh <- data.matrix(read.table(text=strct[dvhStart:length(strct)],
                                      header=FALSE, stringsAsFactors=FALSE,
                                      colClasses=rep("numeric", length(vars3)),
                                      comment.char="", nrows=dvhLen))

        ## rename
        colnames(dvh) <- vars3

        ## catch special case: structVol is 0.0 due to limited precision
        if("volume" %in% vars3) {
            structVol <- if(info$DVHtype == "cumulative") {
                max(c(structVol, dvh[ , "volume"]))
            } else {
                ## reconstruct volumes -> volume is per gray -> mult with bin width
                volBin <- dvh[ , "volume"]*diff(c(-dvh[1, "dose"], dvh[ , "dose"]))
                max(c(structVol, sum(volBin)))
            }
        }

        ## add information we don't have yet: relative/absolute volume
        if((       "volumeRel" %in% vars3) && !("volume"    %in% vars3)) {
            dvh <- cbind(dvh, volume=structVol*(dvh[ , "volumeRel"]/100))
        } else if(("volume"    %in% vars3) && !("volumeRel" %in% vars3)) {
            dvh <- cbind(dvh, volumeRel=100*(dvh[ , "volume"]/structVol))
        }

        ## add information we don't have yet: relative/absolute dose
        ## considering isoDoseRx
        if((    "doseRel" %in% vars3) && !("dose"    %in% vars3)) {
            dvh <- cbind(dvh, dose=dvh[ , "doseRel"]*doseRx / isoDoseRx)
            # (doseRx/(isoDoseRx/100))*(dvh$doseRel/100)
        } else if(("dose" %in% vars3) && !("doseRel" %in% vars3)) {
            dvh <- cbind(dvh, doseRel=dvh[ , "dose"]*isoDoseRx / doseRx)
            # 100*(dvh$dose/(doseRx/(isoDoseRx/100)))
        }

        ## check if dose is increasing
        stopifnot(isIncreasing(dvh))

        DVH <- list(dvh=dvh,
                    patID=info$patID,
                    patName=info$patName,
                    date=info$date,
                    DVHtype=info$DVHtype,
                    plan=info$plan,
                    structure=structure,
                    structVol=structVol,
                    doseUnit=doseUnit,
                    volumeUnit=volumeUnit,
                    doseMin=doseMin,
                    doseMax=doseMax,
                    doseRx=doseRx,
                    doseRxUnit=doseRxUnit,
                    isoDoseRx=isoDoseRx,
                    doseAvg=doseAvg,
                    doseMed=doseMed,
                    doseSD=doseSD)

        ## convert differential DVH (not per unit dose!) to cumulative
        ## and add differential DVH separately
        if(info$DVHtype == "differential") {
            DVH$dvh     <- convertDVH(dvh, toType="cumulative",
                                      toDoseUnit="asis", perDose=FALSE)
            DVH$dvhDiff <- dvh
        }

        ## set class
        class(DVH) <- "DVHs"
        return(DVH)
    }

    ## list of DVH data frames with component name = structure
    info <- list(patID=patID, patName=patName, plan=plan, date=DVHdate,
                 DVHtype=DVHtype, doseRx=doseRx, doseRxUnit=doseRxUnit)
    dvhL <- lapply(structList, getDVH, info)
    dvhL <- Filter(Negate(is.null), dvhL)
    names(dvhL) <- sapply(dvhL, function(y) y$structure)
    if(length(unique(names(dvhL))) < length(dvhL)) {
        warning("Some structures have the same name - this can lead to problems")
    }

    class(dvhL) <- "DVHLst"
    attr(dvhL, which="byPat") <- TRUE

    return(dvhL)
}
