#####---------------------------------------------------------------------------
## combine info from one Pinnacle patient directory = DVHLst
mergePinnaclePat  <- function(x, planInfo=FALSE, courseAsID=FALSE) {
    planInfo  <- as.character(planInfo)
    
    ## we need these files
    infoFiles <- c("DoseInfo.csv", "PatInfo.csv", "PlanInfo.csv", "Data/Info.csv")
    reqFiles  <- paste0(x, "/", infoFiles)
    fRead <- if(!all(file_test(op="-f", reqFiles))) {
        warning(paste("Required files not found in directory", basename(x)))
        NULL
    } else {
        Map(function(y) { read.csv(y, header=TRUE, stringsAsFactors=FALSE, comment.char="") }, reqFiles)
    }
    
    if(!is.null(fRead)) {
        names(fRead) <- basename(names(fRead))
        doseElem <- grep("PrescriptionDose", names(fRead$DoseInfo.csv), ignore.case=TRUE)
        doseUnit <- sub("^[[:alpha:]]+[.](cGy|Gy)$", "\\1", names(fRead$DoseInfo.csv)[doseElem], ignore.case=TRUE)
        
        doseRxIdx   <- grep("Dosis.+(cGy|Gy)",   names(fRead$DoseInfo.csv), ignore.case=TRUE)
        fracDoseIdx <- grep("PrescriptionDose",  names(fRead$DoseInfo.csv), ignore.case=TRUE)
        
        info <- list(patID=removeWS(fRead$PatInfo.csv$MedicalRecordNumber),
                     patName=collWS(trimWS(paste(fRead$PatInfo.csv$FirstName, fRead$PatInfo.csv$LastName))),
                     plan=collWS(trimWS(fRead$PlanInfo.csv$PlanName)),
                     fractionDose=fRead$DoseInfo[[fracDoseIdx]],
                     fractionN=fRead$DoseInfo$NumberOfFractions,
                     doseRx=fRead$DoseInfo[[doseRxIdx]],
                     doseUnit=toupper(doseUnit))
        
        ## read structures
        structs <- split(fRead$Info.csv, f=fRead$Info.csv$RegionOfInterestName)
        DVHraw  <- Map(readLines, paste0(x, "/Data/", fRead$Info.csv$Filename, ".csv"))
        dvhL    <- Map(parsePinnacleDVH, DVHraw, structs, list(info))
        dvhL    <- Filter(Negate(is.null), dvhL)
        names(dvhL) <- sapply(dvhL, function(y) y$structure)
        if(length(unique(names(dvhL))) < length(dvhL)) {
            warning("Some structures have the same name - this can lead to problems")
        }

        class(dvhL) <- "DVHLst"
        attr(dvhL, which="byPat") <- TRUE
        dvhL
    } else {
        NULL
    }

    return(dvhL)
}

#####---------------------------------------------------------------------------
## parse character vector from Pinnacle DVH file
parsePinnacleDVH <- function(x, structInfo, info) {
    ## function to extract one information element from a number of lines
    ## make sure only first : is matched -> not greedy
    getElem <- function(pattern, ll, trim=TRUE, iCase=FALSE, collWS=TRUE) {
        line <- ll[grep(pattern, ll)]
        elem <- sub("^.+?=[[:blank:]]+([[:alnum:][:punct:][:blank:]]+);$", "\\1",
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

    dvhCols   <- as.numeric(getElem("NumberOfDimensions", x))
    dvhLenAdv <- as.numeric(getElem("NumberOfPoints", x))

    ## DVH column headers
    dvhStart  <- grep("Points.+{", x, ignore.case=TRUE, perl=TRUE) + 1
    dvhStop   <- grep("};", x, ignore.case=TRUE, perl=TRUE) - 1
    dvhLen    <- dvhStop - dvhStart + 1
    if((length(dvhStart) < 1L) || length(dvhStop) < 1L) {
        stop("No DVH data found")
    }
    
    if(dvhLenAdv != dvhLen) {
        warning("Advertised DVH length differs from actual length")
    }

    ## extract DVH as a matrix and store preceding information
    ## check if dvh is all blank -> no data
    if(all(!nzchar(x[dvhStart:dvhStop]))) {
        return(NULL)
    }
    
    ## strip trailing ,
    x[dvhStart:dvhStop] <- gsub("(^.+),$", "\\1", x[dvhStart:dvhStop])
    
    dvh <- data.matrix(read.csv(text=x[dvhStart:dvhStop],
                                header=FALSE, stringsAsFactors=FALSE,
                                colClasses=rep("numeric", dvhCols),
                                comment.char="", nrows=dvhLen))

    ## set column names
    ## TODO: this should be exported
    vars3 <- c("dose", "volume")
    colnames(dvh) <- vars3

    ## TODO: this should be exported
    warning("Iso-dose-Rx is assumed to be 100")
    info$isoDoseRx <- 100

    ## structure volume
    volumeIdx <- grep("^Volume", names(structInfo), ignore.case=TRUE)
    info$structVol <- structInfo[[volumeIdx]]

    ## add information we don't have yet: relative/absolute volume
    if((       "volumeRel" %in% vars3) && !("volume"    %in% vars3)) {
        dvh <- cbind(dvh, volume=info$structVol*(dvh[ , "volumeRel"]/100))
    } else if(("volume"    %in% vars3) && !("volumeRel" %in% vars3)) {
        dvh <- cbind(dvh, volumeRel=100*(dvh[ , "volume"]/info$structVol))
    }

    ## add information we don't have yet: relative/absolute dose
    ## considering isoDoseRx
    if((    "doseRel" %in% vars3) && !("dose"    %in% vars3)) {
        dvh <- cbind(dvh, dose=dvh[ , "doseRel"]*info$doseRx / info$isoDoseRx)
        # (info$doseRx/(info$isoDoseRx/100))*(dvh$doseRel/100)
    } else if(("dose" %in% vars3) && !("doseRel" %in% vars3)) {
        dvh <- cbind(dvh, doseRel=dvh[ , "dose"]*info$isoDoseRx / info$doseRx)
        # 100*(dvh$dose/(info$doseRx/(info$isoDoseRx/100)))
    }

    ## TODO: DVHdate
    volumeElem <- grep("Volume", names(structInfo), ignore.case=TRUE)
    volumeUnit <- sub("^[[:alpha:]]+[.](.+)$", "\\1", names(structInfo)[volumeElem], ignore.case=TRUE)
    if(toupper(volumeUnit) == "CCM") {
        volumeUnit <- "CC"
    }

    ## check if dose is increasing
    stopifnot(isIncreasing(dvh))

    ## differential or cumulative DVH
    DVHtype <- dvhType(dvh)

    ## TODO: identify list elements for structVol, dose*, independent of unit
    doseMinIdx <- grep("^DoseMin",  names(structInfo), ignore.case=TRUE)
    doseMaxIdx <- grep("^DoseMax",  names(structInfo), ignore.case=TRUE)
    doseAvgIdx <- grep("^DoseMean", names(structInfo), ignore.case=TRUE)
    DVH <- list(dvh=dvh,
                patName=info$patName,
                patID=info$patID,
                date=info$date,
                DVHtype=DVHtype,
                plan=info$plan,
                fractionDose=info$fractionDose,
                fractionN=info$fractionN,
                quadrant=info$quadrant,
                structure=structInfo$RegionOfInterestName,
                structVol=structInfo[[volumeIdx]],
                doseUnit=info$doseUnit,
                volumeUnit=toupper(volumeUnit),
                doseRx=info$doseRx,
                isoDoseRx=info$isoDoseRx,
                doseMin=structInfo[[doseMinIdx]],
                doseMax=structInfo[[doseMaxIdx]],
                doseAvg=structInfo[[doseAvgIdx]])

    ## convert differential DVH to cumulative
    ## and add differential DVH separately
    if(DVHtype == "differential") {
        DVH$dvh <- convertDVH(dvh, toType="cumulative",
                              toDoseUnit="asis", perDose=FALSE)
        DVH$dvhDiff <- dvh
    }
    
    ## set class
    class(DVH) <- "DVHs"
    return(DVH)
}
