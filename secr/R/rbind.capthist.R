## 2011-09-13 moved from methods.R
## rbind.capthist returns single-session input object unchanged
## rbind.capthist accepts single session objects in pool list
## pooled sessions named correctly
## verify option
## rbind.capthist merges polygons
## 2015-01-11 tweak rownames in rbind.capthist (default component names)
## 2015-01-23 replace long object names (>100ch)

flatten <- function(x) {
## 7/6/2010, 12/9/2011
## x is expected to be a list whose components are non-list and list objects
## create a new list from a concatenation of the non-list objects and the first-order
## components of the lists
    if (!is.list(x))
        stop ("can only flatten a list")
    temp <- lapply(x,
                  function(y)
                      if (is.list(y))
                          return(y)
                      else {
                          return(list(y))
                      }
                 )
    unlist(temp, recursive=F)
}

MS.capthist <- function (...) {
    # make a list of capthist objects, one for each session
    # modified 7/6/2010 for more general input:
    #    multiple single-session objects (as in old version)
    #    list of single-session objects
    #    combination of existing MS objects
    #    combination of existing MS objects and single-session objects?
    # modified 2/7/2010 so always named
    # modified 12/9/2011 so creates names as needed

    dots <- match.call(expand.dots = FALSE)$...
    oldsess <- unlist(sapply(list(...), session))
    MS <- flatten(list(...))
    class(MS) <- c('list', 'capthist')
    if (is.null(names(MS))) {
        names(MS) <- rep("",length(MS))
    }
    if (any(names(MS)=="")) {
        names(MS)[names(MS)==""] <- oldsess[names(MS)==""]
    }
    if (any(duplicated(names(MS)))) {
        warning ("session names replaced to avoid duplication")
        names(MS) <- 1:length(MS)
    }
    session(MS) <- names(MS)
    MS
}
###############################################################################

rbind.capthist <- function (..., renumber = TRUE, pool = NULL, verify = TRUE)
## NOT S3 method (for now) because then naming lost...

{
    dots <- match.call(expand.dots = FALSE)$...
    allargs <- list(...)

    ##############################################################
    ## 2014-11-23
    ## don't understand the purpose of this line, and it breaks secrdesign
    ## temporarily(?) suppress
    ## names(allargs) <- lapply(dots, as.character)

    ## 2015-01-11
    ## Aha! It provides input object names as base for rownames later in:
    ## source <- rep(names(object), sapply(object, nrow))

    inputnames <- lapply(dots, as.character)
    ## if (any(is.na(inputnames) | duplicated(inputnames)))
    if (any(is.na(inputnames) | duplicated(inputnames) | (nchar(inputnames)>100)))
        inputnames <- as.character(1:length(allargs))
    names(allargs) <- inputnames
    ##############################################################

    if (length(dots)==1) object <- allargs[[1]]
    else object <- allargs
    newMCP <- TRUE ## option 2013-11-20

    ## Catch singleton - added 2011-09-12
    if ((length(dots) == 1) & !ms(object) )
        return(object)       ## unchanged!
    if ((length(dots) == 1) & ms(object) & (length(object) == 1) )
        return(object[[1]])  ## unchanged!

    ## Case 1 DEPRECATED 2011-09-12
    ## several lists or a combination of list & elementary capthist objects
    ## concatenate lists, including elementary objects (ignore 'pool')
    ## objects may differ in traps etc.

    if ((length(dots)>1) & any(sapply(allargs, is.list)))
        stop ("invalid input to rbind.capthist; ",
              "use MS.capthist to concatenate sessions")

    ## Case 2
    ## a single MS capthist (i.e. a list)
    ## rbind components as identified in 'pool'
    ## recursive call for each component of 'pool'

    if((length(dots)==1) & (is.list(object)) & !is.null(pool)) {
        if (!is.list(pool)) {
            pool <- list(combined=1:length(object))
            warning ("list not specified, pooling all components")
        }
        else if (any (sapply(unlist(pool),
                             function(x) length(object[[x]])==0)))
                ## prempted by 'subscript out of bounds'
                stop ("invalid pooling indices")

        ## 2011-09-12
        getpooled <- function (x) {
            temphist <- object[x]
            class(temphist) <- c('list', 'capthist')
             ## recursive call
            rbind.capthist(temphist, renumber = renumber, pool=NULL, verify = FALSE)
        }

        temp <- lapply (pool, getpooled)
        if (length(temp)==1) {
            temp <- temp[[1]]
            class(temp) <- 'capthist'
        }
        else {
            class (temp) <- c('list', 'capthist')
            if (is.null(names(pool)) | any(names(pool) == ""))
                names(temp) <- sapply(temp,session)
            else {
                session(temp) <- names(pool)
            }
        }
        ## do it once
        if (verify) {
            verify(temp)
        }
        return(temp)
    }
    else {

    ## Case 3
    ## 1 to several several elementary capthist objects
    ## conventional rbind, given compatible traps, covariates, noccasions
    ## optional renumbering

        check <- function (x) {
            if (!is(x,'capthist'))
                stop ("all arguments must be 'capthist' objects")
            if (is.null(covariates(x)) != is.null(covariates(object[[1]]) ))
                stop ("covariates must be provided for all or none")
            if (is.null(Tu(x)) != is.null(Tu(object[[1]])))
                stop ("unmarked sightings Tu must be provided for all or none")
            if (is.null(Tm(x)) != is.null(Tm(object[[1]])))
                stop ("nonID sightings Tu must be provided for all or none")
            if (any(dim(x)[-1] != dim(object[[1]])[-1]))
                stop ("varying numbers of occasions and/or detectors ",
                      "in rbind.capthist", call. = FALSE)
            notPoolPoly <- !(detector(traps(object[[1]])) %in% c('polygon','polygonX', 'telemetry'))
            if (!identical(traps(x), traps(object[[1]])) & notPoolPoly)
                stop ("cannot pool capthist with different",
                      " detector arrays in rbind.capthist", call. = FALSE)
        }

        ## 2011-09-12 reinstated this check
        sapply (object, check)

        ## form new object
        temp <- abind(..., along = 1)
        class(temp) <- c('capthist')
        trps <- traps(object[[1]])
        mergepoly <- detector(trps) %in% c('polygon','telemetry')
        if (mergepoly) {
        
            srl <- lapply(traps(object), function(x) Polygon(as.matrix(x)))
            tmp <- Polygons(srl,1)
            if (!requireNamespace ('maptools', quietly = TRUE))
                stop("maptools required")
            tmp2 <- maptools::unionSpatialPolygons(SpatialPolygons(list(tmp)), 1)
            ## tmp2 <- unionSpatialPolygons(SpatialPolygons(list(tmp)), 1)
            trps <- as.data.frame(getcoord(tmp2)[[1]])
            rownames(trps) <- 1:nrow(trps)
            class(trps) <- c('traps', 'data.frame')
            detector(trps) <- detector( traps(object[[1]]))
            polyID(trps) <- rep(1,nrow(trps))
            ## note any covariates have been abandoned
        }
        traps(temp) <- trps

        ## 2011-09-13 common covariates
        tempcov <- covariates(object)
        covnamelist <- lapply (tempcov, names)
        covnames <- Reduce(intersect, covnamelist)

        if (length(covnames) > 0) {
            tempcov <- lapply(tempcov, function(x) x[,covnames, drop = FALSE])
            tempcov <- do.call (rbind, tempcov)
            covariates(temp) <- tempcov
        }
        else
            covariates(temp) <- NULL

        ##################################################
        ## sightings
        ## either all-scalar or all-matrix
        if (!is.null(Tu(object[[1]])))
            Tu(temp) <- sum(Tu(object))
        if (!is.null(Tm(object[[1]])))
            Tm(temp) <- sum(Tm(object))
        
        ##################################################
        ## polygon or transect coordinates xy
        tempxy <-  lapply(object, xy)
        xy(temp) <- do.call(rbind, tempxy)

        ##################################################
        ## signal
        tempsig <- lapply(object, signalframe)
        signalframe(temp) <- do.call(rbind, tempsig)
#        if (!is.null(signalframe))
#            stop("rbind.capthist not yet updated for signalframe structure")
#            signal(temp) <- do.call(c, tempsig)
#            cutvals <- sapply(object, function(x) attr(x,'cutval'))
#            attr(temp, 'cutval') <- max(cutvals)
#            temp <- subset(temp, cutval = max(cutvals))
#        }
#        else
#            if (!all(sapply(tempsig, is.null)))
#                stop ("signal attribute missing in one or more sessions")

        ##################################################
        ## messy problem of correct order of detections
        if (!is.null(xy(temp)) | !is.null(signalframe(temp))) {
            occ <- unlist(lapply(object, occasion))
            ID  <- lapply(object, animalID, names = FALSE)
            maxID <- suppressWarnings(sapply(ID, max))
            nID <- sapply(ID, length)
            ID <- unlist(ID)
            uniqueID <- ID + rep(c(0, cumsum(maxID[-length(maxID)])), nID)
            trp <- unlist(lapply(object, trap))
            neworder <- order (occ, uniqueID, trp)
            if (!is.null(xy(temp)))
                xy(temp) <- xy(temp)[neworder,,drop=F]
            if (!is.null(signalframe(temp)))
                signalframe(temp) <- signalframe(temp)[neworder,,drop=F]
        }
        ##################################################

        ## name new sessions
        session (temp) <- paste(names(object), collapse='+')

        ## optionally construct unique row names
        if (renumber) {
            ID <- unlist(sapply(object, rownames))
            source <- rep(names(object), sapply(object, nrow))
            rownames(temp) <- paste(source, ID, sep='.')
        }

        ## optionally verify
        if (verify) {
            verify(temp)
        }
        temp
    }
}
###############################################################################

merge.capthist <- function (..., renumber = TRUE, remove.dupl.traps = FALSE,
                            concurrent = TRUE, verify = TRUE, tol = 1)
{
    dots <- match.call(expand.dots = FALSE)$...
    allargs <- list(...)
    names(allargs) <- lapply(dots, as.character)
    if (length(dots)==1) object <- allargs[[1]]
    else object <- allargs

    ## Catch singleton - added 2011-09-12
    if ((length(dots) == 1) & !ms(object) )
        return(object)       ## unchanged!
    if ((length(dots) == 1) & ms(object) & (length(object) == 1) )
        return(object[[1]])  ## unchanged!

    ## several lists or a combination of list & elementary capthist objects
    if ((length(dots)>1) & any(sapply(allargs, is.list)))
        stop ("invalid input to merge.capthist; ",
              "use MS.capthist to concatenate sessions")

    ## 1 to several several elementary capthist objects
    ## conventional rbind, given compatible traps, covariates, noccasions
    ## optional renumbering

    check <- function (x) {
      if (!is(x,'capthist'))
          stop ("all arguments must be 'capthist' objects")
    }
    sapply (object, check)

    ## resolve traps
    temptrp <- traps(object)
    newtraps <- do.call(rbind, temptrp)
    trp <- sapply(temptrp, nrow)
    if (remove.dupl.traps) {
        temp <- as.matrix(dist(newtraps))
        diag(temp) <- 1e10
        temp[upper.tri(temp)] <- 1e10
        drop <- apply(temp, 1, function(x) any(x<0.001))
        temp <- temp[,!drop]


############# use matrix to find matches with distance < tol
        stop()

        trpindex <- 1
    }
    else {
        ntrp <- sum(trp)
        trpindex <- split(1:ntrp, rep(1:length(trp), trp))
    }

    ## resolve occasions
    occ <- sapply (object, ncol)
    if (concurrent) {
        nocc <- max(occ)
        occindex <- lapply(occ, function(S) 1:S)
    }
    else {  ## consecutive
        nocc <- sum(occ)
        occindex <- split(1:nocc, rep(1:length(occ), occ))
    }

    ## resolve animals

    ## form new object
    temp <- abind(..., along=1)
    class(temp) <- c('capthist')
    traps(temp) <- traps(object[[1]])

    ## 2011-09-13 common covariates
    tempcov <- covariates(object)
    covnamelist <- lapply (tempcov, names)
    covnames <- Reduce(intersect, covnamelist)
    if (length(covnames) > 0) {
        tempcov <- lapply(tempcov, function(x) x[,covnames, drop = FALSE])
        tempcov <- do.call (rbind, tempcov)
        covariates(temp) <- tempcov
    }
    else
        covariates(temp) <- NULL

    ## polygon or transect coordinates xy
    tempxy <-  lapply(object, xy)
    xy(temp) <- do.call(rbind, tempxy)

    ## signal
    tempsig <- lapply(object, signal)
    if (!any(sapply(tempsig, is.null))) {
        signal(temp) <- do.call(c, tempsig)
        cutvals <- sapply(object, function(x) attr(x,'cutval'))
        attr(temp, 'cutval') <- max(cutvals)
        temp <- subset(temp, cutval = max(cutvals))
    }
    else
        if (!all(sapply(tempsig, is.null)))
            stop ("signal attribute missing in one or more sessions")

    ##################################################
    ## messy problem of correct order of detections
    if (!is.null(xy(temp)) | !is.null(signal(temp))) {
        occ <- unlist(lapply(object, occasion))
        ID  <- lapply(object, animalID, names = FALSE)
        maxID <- sapply(ID, max)
        nID <- sapply(ID, length)
        ID <- unlist(ID)
        uniqueID <- ID + rep(c(0, cumsum(maxID[-length(maxID)])), nID)
        trp <- unlist(lapply(object, trap))
        neworder <- order (occ, uniqueID, trp)
        if (!is.null(xy(temp)))
            xy(temp) <- xy(temp)[neworder,,drop=F]
        if (!is.null(signal(temp)))
            signal(temp) <- signal(temp)[neworder]
    }
    ##################################################

    ## name new sessions
    session (temp) <- paste(names(object), collapse='+')

    ## optionally construct unique row names
    if (renumber) {
        ID <- unlist(sapply(object, rownames))
        source <- rep(names(object), sapply(object, nrow))
        rownames(temp) <- paste(source, ID, sep='.')
    }

    ## optionally verify
    if (verify) {
        verify(temp)
    }
    temp
}
###############################################################################

