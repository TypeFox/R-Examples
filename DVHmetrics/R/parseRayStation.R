#####---------------------------------------------------------------------------
## parse character vector from RayStation DVH file
parseRayStation <- function(x, planInfo=FALSE, courseAsID=FALSE) {
    planInfo <- as.character(planInfo)

    ## function to extract one information element from a number of lines
    ## make sure only first : is matched -> not greedy
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

    getDoseUnit <- function(ll) {
        line <- ll[grep("^#Dose unit:", ll)]
        elem <- sub("^.+:[[:blank:]]+(GY|CGY)$", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimWS(elem))
    }

    sStart <- grep("^#RoiName:[[:alnum:][:punct:]]", x) # start of sections
    sLen   <- diff(c(sStart, length(x)+1))              # length of sections
    if((length(sLen) < 1L) || all(sLen < 1L)) {
        stop("No structures found")
    }

    structList <- Map(function(start, len) x[start:(start+len-1)], sStart, sLen)

    ## extract file header and header info
    header  <- x[seq_len(sStart[1]-1)]            # header
    patName <- getElem("#PatientName:", header)   # patient name
    patID   <- getElem("^#PatientId:",  header)   # patient id
    plan    <- getElem("^#Dosename:",   header)   # treatment plan
    DVHdate <- NA_character_

    doseRx <- if(tolower(planInfo) == "doserx") {
        doseRxUnit <- toupper(sub("^.+[[:blank:]][.[:digit:]]+(c?Gy).*$", "\\1",
                                  plan, perl=TRUE, ignore.case=TRUE))

        if(!grepl("^(GY|CGY)$", doseRxUnit)) {
            warning("Could not determine dose Rx unit")
            doseRxUnit <- NA_character_
        }

        drx <- sub("^.+[[:blank:]]([.[:digit:]]+)c?Gy.*$", "\\1",
                   plan, perl=TRUE, ignore.case=TRUE)
        as.numeric(drx)
    } else {
        doseRxUnit <- NA_character_
        warning("No info on prescribed dose")
        NA_real_
    }

    isoDoseRx <- if(tolower(planInfo) == "doserx") {
        warning("Iso-dose-Rx is assumed to be 100")
        100
    } else {
        warning("No info on % for dose")
        NA_real_
    }

    ## extract DVH from one structure section and store in a list
    ## with DVH itself as a matrix
    getDVH <- function(strct, info) {
        ## extract information from info list
        doseRx     <- info$doseRx
        doseRxUnit <- info$doseRxUnit
        isoDoseRx  <- info$isoDoseRx

        ## extract structure, volume, dose min, max, mean, median and sd
        structure  <- getElem("^#RoiName:", strct)
        structVol  <- NA_real_
        doseMin    <- NA_real_
        doseMax    <- NA_real_
        doseAvg    <- NA_real_
        doseMed    <- NA_real_
        doseMode   <- NA_real_
        doseSD     <- NA_real_
        volumeUnit <- "CC"

        ## set dose unit
        doseUnit <- getDoseUnit(strct)
        if(!grepl("^(GY|CGY)$", doseUnit)) {
            warning("Could not determine dose measurement unit")
            doseUnit <- NA_character_
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

        ## find DVH
        ## DVH column headers
        colHead  <- grep("^#Dose unit:", strct, ignore.case=TRUE, perl=TRUE)
        dvhStart <- colHead+1                 # first numeric line of DVH
        dvhLen   <- length(strct) - dvhStart + 1
        if((length(dvhLen) < 1L) || dvhLen < 1L) {
            stop("No DVH data found")
        }

        ## extract DVH as a matrix and store preceding information
        ## read line length(strct) for cases where file does not end with a
        ## blank line -> this will then be last DVH line, otherwise blank
        ## check if dvh is all blank -> no data
        if(all(!nzchar(strct[dvhStart:length(strct)]))) {
            return(NULL)
        }
        
        dvh <- data.matrix(read.table(text=strct[dvhStart:length(strct)],
                                      header=FALSE, stringsAsFactors=FALSE,
                                      colClasses=rep("numeric", 2),
                                      comment.char="", nrows=dvhLen))

        ## rename columns
        colnames(dvh) <- c("dose", "volumeRel")

        ## add information we don't have yet: absolute volume
        dvh <- cbind(dvh, volume=structVol*(dvh[ , "volumeRel"]/100))

        ## add information we don't have yet: relative dose
        ## considering isoDoseRx
        dvh <- cbind(dvh, doseRel=dvh[ , "dose"]*isoDoseRx / doseRx)
        
        ## check if dose is increasing
        stopifnot(isIncreasing(dvh))

        ## differential or cumulative DVH
        DVHtype <- dvhType(dvh)

        DVH <- list(dvh=dvh,
                    patName=info$patName,
                    patID=info$patID,
                    date=info$date,
                    DVHtype=DVHtype,
                    plan=info$plan,
                    course=info$course,
                    quadrant=info$quadrant,
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
                    doseMode=doseMode,
                    doseSD=doseSD)

        ## convert differential DVH to cumulative
        ## and add differential DVH separately
        if(DVHtype == "differential") {
            warning("I assume differential DVH is per unit dose\nbut I have no information on this")
            DVH$dvh <- convertDVH(dvh, toType="cumulative",
                                  toDoseUnit="asis", perDose=TRUE)
            DVH$dvhDiff <- dvh
        }

        ## set class
        class(DVH) <- "DVHs"
        return(DVH)
    }

    ## list of DVH data frames with component name = structure
    info <- list(patID=patID, patName=patName, date=DVHdate, plan=plan,
                 doseRx=doseRx, doseRxUnit=doseRxUnit, isoDoseRx=isoDoseRx)
    dvhL <- lapply(structList, getDVH, info=info)
    dvhL <- Filter(Negate(is.null), dvhL)
    names(dvhL) <- sapply(dvhL, function(y) y$structure)
    if(length(unique(names(dvhL))) < length(dvhL)) {
        warning("Some structures have the same name - this can lead to problems")
    }

    class(dvhL) <- "DVHLst"
    attr(dvhL, which="byPat") <- TRUE

    return(dvhL)
}
