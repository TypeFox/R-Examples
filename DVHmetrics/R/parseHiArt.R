#####---------------------------------------------------------------------------
## parse character vector from Tomo HiArt DVH file
parseHiArt <- function(x, planInfo=FALSE, courseAsID=FALSE) {
    planInfo <- as.character(planInfo)

    ## find columns for structure, dose, volume
    vars <- as.matrix(read.csv(text=x[1], header=FALSE,
                               stringsAsFactors=FALSE, comment.char="")[1, ])

    structIdx <- seq(1, length(vars), by=3)
    doseIdx   <- seq(2, length(vars), by=3)
    volumeIdx <- seq(3, length(vars), by=3)

    ## Tomo HiArt has no patient name or id in file
    patName   <- gsub("[^a-z0-9]", "\\1", tempfile(pattern="", tmpdir=""))
    patID     <- gsub("[^a-z0-9]", "\\1", tempfile(pattern="", tmpdir="")) # as.character(trunc(runif(1, 0, 10000001)))
    plan      <- NA_character_
    doseRx    <- NA_real_
    isoDoseRx <- NA_real_
    DVHdate   <- NA_character_
    
    ## get dose and volume units
    varDose  <- vars[grepl("Dose",   vars)][1]
    varVol   <- vars[grepl("Volume", vars)][1]
    if(grepl("Relative", varDose, ignore.case=TRUE)) {
        ## TODO: need example file for this
        warning("HiArt files with relative dose are not implemented")
        isDoseRel <- TRUE
        doseUnit  <- "PERCENT"
    } else {
        isDoseRel <- FALSE
        doseUnit  <- toupper(sub("Dose \\((.+)\\)", "\\1", varDose))
    }

    if(grepl("Relative", varVol, ignore.case=TRUE)) {
        isVolRel   <- TRUE
        volumeUnit <- "PERCENT"
    } else {
        warning("HiArt files with absolute volume are not implemented")
        ## TODO: need example file for this
        isVolRel   <- FALSE
        volumeUnit <- NA_character_
    }

    ## read all data
    datAll <- data.matrix(read.csv(text=x[-1], header=FALSE,
                                   stringsAsFactors=FALSE, comment.char=""))

    ## extract DVH from one structure section and store in a list
    ## with DVH itself as a matrix
    getDVH <- function(strIdx, dIdx, vIdx, info) {
        structure <- sub("(.+)\\(STANDARD\\)$", "\\1", vars[strIdx])

        ## extract DVH as a matrix and set variable names
        dvh <- datAll[ , c(dIdx, vIdx)]
        haveVars <- if(isVolRel) {
            if(isDoseRel) {
                c("doseRel", "volumeRel")
            } else {
                c("dose",    "volumeRel")
            }
        } else {
            if(isDoseRel) {
                c("doseRel", "volume")
            } else {
                c("dose",    "volume")
            }
        }
        
        colnames(dvh) <- haveVars

        ## add information we don't have yet
        ## relative/absolute volume/dose
        if(!("volume" %in% haveVars)) {
            dvh <- cbind(dvh, volume=NA_real_)
        }

        if(!("volumeRel" %in% haveVars)) {
            dvh <- cbind(dvh, volumeRel=NA_real_)
        }

        if(!("dose" %in% haveVars)) {
            dvh <- cbind(dvh, dose=NA_real_)
        }

        if(!("doseRel" %in% haveVars)) {
            dvh <- cbind(dvh, doseRel=NA_real_)
        }

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
                    structure=structure,
                    structVol=NA_real_,
                    doseUnit=info$doseUnit,
                    volumeUnit=info$volumeUnit,
                    doseRx=doseRx,
                    isoDoseRx=isoDoseRx,
                    doseMin=NA_real_,
                    doseMax=NA_real_,
                    doseAvg=NA_real_,
                    doseMed=NA_real_,
                    doseMode=NA_real_,
                    doseSD=NA_real_)

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
    info <- list(patID=patID, patName=patName, date=DVHdate,
                 plan=plan, doseRx=doseRx, isoDoseRx=isoDoseRx,
                 doseUnit=doseUnit, volumeUnit=volumeUnit)
    dvhL <- Map(getDVH, structIdx, doseIdx, volumeIdx, info=list(info))
    dvhL <- Filter(Negate(is.null), dvhL)
    names(dvhL) <- sapply(dvhL, function(y) y$structure)
    if(length(unique(names(dvhL))) < length(dvhL)) {
        warning("Some structures have the same name - this can lead to problems")
    }

    class(dvhL) <- "DVHLst"
    attr(dvhL, which="byPat") <- TRUE

    return(dvhL)
}
