## parse vector of constraint strings into components
## optionally convert dose unit and volume unit
parseConstraint <-
function(x, doseUnit=NULL, volUnit=NULL) {
    if(is.list(x)) { x <- unlist(x) }

    ## remove whitespace and convert to upper case
    x <- toupper(removeWS(x))

    ## special metrics that are recognized when prefixed with D, e.g., DMEAN, DEUD
    specMetr <- getSpecialMetrics()

    ## split constraints into metric and comparison
    cnstrSpl <- strsplit(x, ">|<|<=|>=")
    metric   <- vapply(cnstrSpl, function(x) head(x, n=1), character(1)) # metric part
    comp     <- vapply(cnstrSpl, function(x) tail(x, n=1), character(1)) # comparison part
    cmp      <- sub("^.+(>|<|<=|>=).+", "\\1", x) # the comparison

    ## metric details
    pm <- parseMetric(metric, doseUnit=doseUnit, volUnit=volUnit)
    
    ## identify cases of special dose metric DMEAN etc.
    specDose <- pm$valRef %in% specMetr

    ## comparison details
    pattern    <- "^([.[:digit:]]+)([%]|GY|CGY|CC)$"
    valCmp     <- sub(pattern, "\\1", comp)
    valCmpNum  <- as.numeric(valCmp)
    unitCmp    <- sub(pattern, "\\2", comp)

    ## convert absolute dose units if requested
    if(!is.null(doseUnit)) {
        doseUnit <- toupper(removeWS(doseUnit))

        ## output is dose -> don't convert relative dose
        ## consider DMEAN case
        idxD      <- (pm$DV == "D") & (unitCmp != "%") & !specDose
        valCmpNum <- ifelse(idxD,
                            valCmpNum*suppressWarnings(getConvFac(paste0(unitCmp, "2", doseUnit))),
                            valCmpNum)
        
        valCmp  <- ifelse(idxD, as.character(valCmpNum), valCmp)
        unitCmp <- ifelse(idxD, doseUnit, unitCmp)
    }
    
    ## convert absolute volume units if requested
    if(!is.null(volUnit)) {
        volUnit <- toupper(removeWS(volUnit))

        ## output is volume -> don't convert relative volume
        idxV      <- (pm$DV == "V") & (unitCmp != "%")
        valCmpNum <- ifelse(idxV,
                            valCmpNum*suppressWarnings(getConvFac(paste0(unitCmp, "2", volUnit))),
                            valCmpNum)

        valCmp  <- ifelse(idxV, as.character(valCmpNum), valCmp)
        unitCmp <- ifelse(idxV, volUnit, unitCmp)
    }

    ## clean constraint string
    DV         <- pm$DV
    unitRef    <- pm$unitRef
    metric     <- ifelse(is.na(pm$unitDV), paste0(pm$metric, "_", unitCmp), pm$metric)
    constraint <- paste0(metric, " ", cmp, " ", valCmp, unitCmp)

    ## check comparison validity
    validPattern <- grepl(pattern, comp)
    validUnitCmp <- ((pm$DV == "D") & (unitCmp %in% c("%", "GY", "CGY"))) |
                    ((pm$DV == "V") & (unitCmp %in% c("%", "CC")))
    validCmp     <- cmp %in% c(">", "<", "<=", ">=")
    validUnitDV  <- ifelse(is.na(pm$unitDV), TRUE, pm$unitDV == unitCmp)
    valid        <- pm$valid & validPattern & validUnitCmp & validCmp & validUnitDV

    if(!all(valid)) {
        warning(c("Constraint ", paste(constraint[!valid], collapse=", "),
                  " is invalid"))
    }

    ## inverse metric
    DVinv      <- ifelse(pm$DV == "V", "D", "V")
    valRefInv  <- valCmp
    unitRefInv <- unitCmp

    ## special cases
    metricInv <- ifelse((pm$DV == "D") &
                        (pm$valRef %in% specMetr),
                        NA_character_,
                        paste0(DVinv, valRefInv, unitRefInv, "_", pm$unitRef))

    return(list(constraint=constraint, valid=valid, metric=metric,
                DV=pm$DV, valRef=pm$valRefNum, unitRef=pm$unitRef,
                cmp=cmp, valCmp=valCmpNum, unitCmp=unitCmp, metricInv=metricInv))
}

## expand globs in constraint definition
## constr is a constraint data frame
## we also need IDs and structures to match against
## dvhID     is a list - 1 component for each structure
## dvhStruct is a list - 1 component for each dvhID
expandConstraint <-
function(constr, byPat=TRUE, dvhID, dvhStruct) {
    constrDF <- as.data.frame(constr, stringsAsFactors=FALSE)

    ## if constr has 1 column -> this is the constraint
    if(ncol(constrDF) == 1L) {
        names(constrDF) <- "constraint"
    }

    ## if constr has >= 2 columns -> we need variable names
    dfNames <- names(constrDF)

    ## if we don't have IDs -> match to all
    if(!("patID" %in% dfNames)) {
        constrDF$patID <- "*"
    }
    
    ## if we don't have structures -> match to all
    if(!("structure" %in% dfNames)) {
        constrDF$structure <- "*"
    }

    retDF <- if(byPat) {
        if(missing(dvhID)) { dvhID <- "" }
        if(length(dvhID) > length(dvhStruct)) {
            warning("More IDs than structure vectors -> IDs get truncated")
            dvhID <- dvhID[seq_along(dvhStruct)]
        }
    
        if(length(dvhID) < length(dvhStruct)) {
            warning("Fewer IDs than structure vectors -> IDs get recycled")
            dvhID <- rep(dvhID, length(dvhStruct))[seq_along(dvhStruct)]
        }

        ## 1st expand IDs
        ## match IDs in constrDF to those available in dvhID
        ## replace patID pattern with actual patIDs
        IDpat  <- glob2rx(constrDF$patID)
        IDL    <- lapply(IDpat, function(x) dvhID[grepl(x, dvhID, ignore.case=TRUE)])
        IDlens <- vapply(IDL, length, 1)
        expID  <- constrDF[rep(seq_len(nrow(constrDF)), IDlens), ]
        expID$patID <- unlist(IDL)
        
        ## 2nd within IDs: expand structures
        ## split according to ID
        expIDL <- split(expID, expID$patID)
    
        ## for one ID: expand structures
        expandStructs <- function(x, struct) {
            ## match structures in x to those available in struct
            ## replace structure pattern with actual structures
            strPat  <- glob2rx(x$structure)
            strL    <- lapply(strPat, function(y) struct[grepl(y, struct, ignore.case=TRUE)])
            strLens <- vapply(strL, length, 1)
            expStr  <- x[rep(seq_len(nrow(x)), strLens), ]
            expStr$structure <- unlist(strL)
            expStr
        }
    
        ## only for matched patient IDs
        dvhStructSub <- dvhStruct[unique(expID$patID)]
        expStrL <- Map(expandStructs, expIDL, dvhStructSub)
        do.call("rbind", expStrL)
    } else {
        ## byPat=FALSE
        if(missing(dvhStruct)) { dvhStruct <- "" }
        if(length(dvhStruct) > length(dvhID)) {
            warning("More structures than ID vectors -> structures get truncated")
            dvhStruct <- dvhStruct[seq_along(dvhID)]
        }
    
        if(length(dvhStruct) < length(dvhID)) {
            warning("Fewer structures than ID vectors -> structures get recycled")
            dvhID <- rep(dvhStruct, length(dvhID))[seq_along(dvhID)]
        }

        ## 1st expand structures
        ## match structures in constrDF to those available in dvhStruct
        ## replace structure pattern with actual structures
        strPat  <- glob2rx(constrDF$structure)
        strL    <- lapply(strPat, function(x) dvhStruct[grepl(x, dvhStruct, ignore.case=TRUE)])
        strLens <- vapply(strL, length, 1)
        expStr  <- constrDF[rep(seq_len(nrow(constrDF)), strLens), ]
        expStr$structure <- unlist(strL)
        
        ## 2nd within structures: expand IDs
        ## split according to structure
        expStrL <- split(expStr, expStr$structure)
    
        ## for one structure: expand IDs
        expandIDs <- function(x, ID) {
            ## match IDs in x to those available in ID
            ## replace ID pattern with actual IDs
            IDpat  <- glob2rx(x$patID)
            IDL    <- lapply(IDpat, function(y) ID[grepl(y, ID, ignore.case=TRUE)])
            IDlens <- vapply(IDL, length, 1)
            expID  <- x[rep(seq_len(nrow(x)), IDlens), ]
            expID$patID <- unlist(IDL)
            expID
        }
    
        ## only for matched structures
        dvhIDsub <- dvhID[unique(expStr$structure)]
        expIDL   <- Map(expandIDs, expStrL, dvhIDsub)
        do.call("rbind", expIDL)
    }
    
    return(unique(retDF))
}

## convert constraint data frame for 1 id/structure to hierarchical list
constrDF2L1 <-
function(x, byPat=TRUE, expand=FALSE, ...) {
    xExp <- if(expand) {
        ## first expand globs
        ## -> ... should then be IDs or structures to match against
        expandConstraint(x, byPat=byPat, ...)
    } else {
        x
    }

    constr <- if(byPat) {
        lapply(split(xExp, xExp$structure), function(y) { 
            lapply(split(y, y$constraint),  function(z) z$constraint) })
    } else {
        lapply(split(xExp, xExp$patID),     function(y) {
            lapply(split(y, y$constraint),  function(z) z$constraint) })
    }
    
    attr(constr, which="byPat") <- byPat
    constr
}

## convert constraint data frame for many ids to hierarchical list
constrDF2L <-
function(x, byPat=TRUE, expand=FALSE, ...) {
    xExp <- if(expand) {
        ## first expand globs
        ## -> ... should then be IDs and structures to match against
        expandConstraint(x, byPat=byPat, ...)
    } else {
        x
    }

    ## split by patient or by structure
    xSpl <- if(byPat) {
        lapply(split(xExp, xExp$patID),     function(y) split(y, y$structure))
    } else {
        lapply(split(xExp, xExp$structure), function(y) split(y, y$patID))
    }

    ## within each patient/structure -> split by constraint
    constr <- lapply(xSpl, function(y) lapply(y, function(z) {
        lapply(split(z, z$constraint), function(zz) zz$constraint) }))

    attr(constr, which="byPat") <- byPat
    constr
}
