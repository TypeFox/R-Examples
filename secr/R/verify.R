############################################################################################
## package 'secr'
## verify.R
## 2009 09 18, 2009 09 19, 2009 09 20, 2009 10 02, 2009 11 05, 2009 11 13
## 2010 05 02 removed erroneous ref to 'areabinary' detector
## 2011 02 15 validpoly
## 2011 03 19 adjustments to allow unmarked; need more complete check of unmarked
## 2012 02 02 signal check applied at level of whole sound
## 2012-10-22 xylist checked
## 2013-05-09 tweak to avoid error in checkcovariatelevels when no covariates
## 2015-01-06 stop if mixture of NULL and non-NULL covariates
## 2015-10-03 resight data
## 2015-10-12 verify.traps error messages
## 2015-11-02 xyinpoly moved to utility.R
############################################################################################

verify <- function (object, report, ...) UseMethod("verify")

verify.default <- function (object, report, ...) {
  cat ('no verify method for objects of class', class(object), '\n')
}
############################################################################################

## from is.integer help page
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

overlapcells <- function (xy) {
    vertexinside <- function (a,b) {
        OK <- FALSE
        for (k in 1:4) {
            temp <- .C('inside',  PACKAGE = 'secr',
                as.double (a[k,]),
                as.integer (0),
                as.integer (3),
                as.integer (5),
                as.double (b),
                result = integer(1))
            if (any(as.logical(temp$result))) OK <- TRUE
            temp <- .C('inside',  PACKAGE = 'secr',
                as.double (b[k,]),
                as.integer (0),
                as.integer (3),
                as.integer (5),
                as.double (a),
                result = integer(1))
            if (any(as.logical(temp$result))) OK <- TRUE
        }
        OK
    }
    spacex <- attr(xy, 'spacex')
    spacey <- attr(xy, 'spacey')
    fuzz <- 1e-10
    spx2 <- spacex/2 - fuzz
    spy2 <- spacey/2 - fuzz
    xy <- as.matrix(xy)
    nr <- nrow(xy)
    if (nr<2)
        FALSE
    else {
        pixel <- matrix(ncol = 2, c(-spx2,-spx2,spx2,spx2,-spx2,-spy2,spy2,spy2,-spy2,-spy2))
        overlap <- matrix(FALSE, nrow = nr, ncol = nr)
        for (i in 1:(nr-1))
            for (j in (i+1):nr)
                {
                    verti <- t(apply(pixel, 1, function (x) xy[i,] + x))
                    vertj <- t(apply(pixel, 1, function (x) xy[j,] + x))
                    if ((length(verti)>0) && (length(vertj)>0))
                    overlap[i,j] <- vertexinside(verti,vertj)
                }
        any (overlap, na.rm=T)
    }
}
############################################################################################

overlappoly <- function (xy, polyID) {
    vertexinside <- function (a,b) {
        OK <- FALSE
        n.a <- nrow(a)
        n.b <- nrow(b)
        a <- as.matrix(a)
        b <- as.matrix(b)
        for (k in 1:n.a) {
            temp <- .C('inside',  PACKAGE = 'secr',
                as.double (a[k,]),
                as.integer (0),
                as.integer (n.b-1),
                as.integer (n.b),
                as.double (b),
                result = integer(1))
            if (any(as.logical(temp$result))) OK <- TRUE
        }
        for (k in 1:n.b) {
            temp <- .C('inside',  PACKAGE = 'secr',
                as.double (b[k,]),
                as.integer (0),
                as.integer (n.a-1),
                as.integer (n.a),
                as.double (a),
                result = integer(1))
            if (any(as.logical(temp$result))) OK <- TRUE
        }
        OK
    }
    lxy <- split (xy, levels(polyID))
    nr <- length(lxy)
    if (nr<2)
        FALSE
    else {
        overlap <- matrix(FALSE, nrow = nr, ncol = nr)
        for (i in 1:(nr-1))
            for (j in (i+1):nr)
                {
                    overlap[i,j] <- vertexinside(lxy[[i]], lxy[[j]])
                }
        any (overlap, na.rm=T)
    }
}
############################################################################################

validpoly <- function (xy, polyID, nx = 500) {
    ## check intersections of perimeters with vertical lines
    OKpoly <- function (xy) {
        cross <- rep(0, nx)
        x <- seq(min(xy[,1]), max(xy[,1]), length=nx)
        for (i in 1:nx) {
            tempx <- xy[,1] - x[i]
            cross[i] <- sum (sign(tempx[-1]) != sign(tempx[-length(tempx)]))
        }
        cross <= 2
    }
    lxy <- split (xy, levels(polyID))
    temp <- lapply(lxy, OKpoly)
    all(unlist(temp))
}
############################################################################################

xyontransect <- function (xy, trps, tol=0.01) {
    ptontransect <- function (i,k) {
        ## is point i on transect k?
        transectxy <- as.matrix(lxy[[k]])
        nr <- nrow(transectxy)
        temp <- .C('ontransect',  PACKAGE = 'secr',
            as.double (xy[i,]),
            as.integer (0),
            as.integer (nr-1),
            as.integer (nr),
            as.double (transectxy),
            as.double (tol),
            result = integer(1))
        as.logical(temp$result)
    }
    lxy <- split (trps, levels(transectID(trps)))
    firsttransect <- function (i) {
        for (k in 1:length(lxy))
            if (ptontransect(i,k)) return(k)
        0
    }
    sapply(1:nrow(xy), firsttransect)
}
############################################################################################

checkcovariatelevels <- function (cov) {
    ## cov is a list of dataframes of covariates
    ## tweak 2013-05-09 to avoid error when no covariates
    xfactor <- function(x) {
        if ((nrow(x)>0) & (ncol(x)>0))
            names(x)[sapply(x,is.factor)]
        else
            NULL
    }
    factornames <- sapply(cov, xfactor)
    factornames <- unique(unlist(factornames))
    if (length(factornames) > 0) {
        if (any(sapply(cov, function(x) !all(factornames %in% names(x))))) {
            return(FALSE)
        }
        else {
            checklevels <- function(variable) {
                baselevels <- levels(cov[[1]][,variable])
                lev <- lapply(cov, function(x) levels(x[,variable]))
                all(sapply(lev[-1], function(x) identical(x,baselevels)))
            }
            return(all(sapply(factornames, checklevels)))
        }
    }
    else
        return(TRUE)
}
############################################################################################

verify.traps <- function (object, report = 2, ...) {

## Check internal consistency of 'traps' object
##
## -- Number of rows in dataframe of detector covariates differs expected
## -- Number of detectors in usage matrix differs from expected
## -- Occasions with no used detectors

    if (!inherits(object, 'traps')) {
         stop ("object must be of class 'traps'")
    }

    if (inherits(object, 'list')) {
        temp <- lapply (object, verify, report = min(report,1))
        anyerrors <- any(sapply(temp, function(x) x$errors))

        ## check covariate factor levels conform across sessions
        
        if (!all(sapply(covariates(object), is.null))) {
            if (any(sapply(covariates(object), is.null)))
                stop ("mixture of NULL and non-NULL trap covariates in different sessions")
            trapcovariatelevelsOK <- checkcovariatelevels(covariates(object))
            if (!trapcovariatelevelsOK & report>0) {
                warning ('Levels of factor trap covariate(s) differ between sessions')
            }
        }

        if ((report == 2) && !anyerrors)
            cat('No errors found :-)\n')
        invisible(list(errors = anyerrors, bysession = temp))
    }
    else {

        single <- detector(object) %in% c('single')
        area <- FALSE
        poly <- detector(object) %in% c('polygon','polygonX')
        telem <- detector(object) == 'telemetry'

        usagedetectorsOK <- TRUE
        usagenonzeroOK <- TRUE
        markoccOK <- TRUE
        markoccOK2 <- TRUE
        areaOK <- TRUE
        polyIDOK <- TRUE
        polyconvexOK <- TRUE
        trapcovariatesOK <- TRUE

        if (!is.null(covariates(object)))
            if ((ncol(covariates(object)) == 0 ) |
                (nrow(covariates(object)) == 0 )) covariates(object) <- NULL

        ## 1
        trapNAOK <- !any(is.na(object))

        ## 2
        if (!is.null(covariates(object)))
            trapcovariatesOK <- nrow(covariates(object)) == ndetector(object)

        ## 'usage' of traps
        if (!is.null(usage(object))) {
            ## 3
            usagedetectorsOK <- nrow(usage(object)) == ndetector(object)

            ## 4
            usagecount <- apply(usage(object),2,sum)
            usagenonzeroOK <- !any(usagecount == 0)
        }
        else usagecount <- rep(NA, ncol(object))

        ## 5
        if (sighting(object)) {
            if (!is.null(usage(object))) 
                markoccOK <- length(markocc(object)) == ncol(usage(object))
            markoccOK2 <- all(markocc(object) %in% c(-1,0,1))
        }
        
        ## 6
        if (area) {
            areaOK <- !is.na(searcharea(object))
            areaOK <- areaOK & !overlapcells(object)
        }
        else
        if (poly | telem) {
            areaOK <- !overlappoly (object, levels(polyID(object)))
        }

        ## 7
        if (poly | telem) {
            polyIDOK <- (length(polyID(object)) == nrow(object)) &
                is.factor(polyID(object))
        }

        ## 8
        if (poly) {
            polyconvexOK <- validpoly (object, polyID(object))
        }

        errors <- !all(c(trapNAOK, trapcovariatesOK, usagedetectorsOK, usagenonzeroOK,
                         areaOK, polyIDOK, polyconvexOK))

        if (report > 0) {
            if (errors) {
                if (!trapNAOK) {
                    cat ('Missing detector coordinates not allowed\n')
                }
                if (!trapcovariatesOK) {
                    cat ('Wrong number of rows in dataframe of detector covariates\n')
                    cat ('traps :', ndetector(object), 'detectors\n')
                    cat ('covariates :', nrow(covariates(object)), 'detectors\n')
                }
                if (!usagedetectorsOK) {
                    cat ('Conflicting number of detectors in usage matrix\n')
                    cat ('traps :', ndetector(object), 'detectors\n')
                    cat ('usage(traps) :', nrow(usage(object)), 'detectors\n')
                }
                if (!usagenonzeroOK) {
                    cat ("Occasions when no detectors 'used'\n")
                    cat ((1:length(usagecount))[usagecount==0], '\n')
                }
                if (!markoccOK) {
                    cat ("Conflicting number of occasions in markocc and usage \n")
                    cat ('markocc : ', length(markocc(object)), 'occasions \n')
                    cat ('usage : ', ncol(usage(object)), 'occasions \n')
                }
                if (!markoccOK2) {
                    cat ("Values in markocc should be -1,0 or 1 \n")
                    cat ('markocc(object) \n')
                }
                if (!areaOK) {
                    cat ("Search areas overlap, or no search area specified \n")
                }
                if (!polyIDOK) {
                    cat ("Invalid polyID \n")
                }
                if (!polyconvexOK) {
                    cat ("The boundary of at least one polygon is concave east-west \n")
                }
            }
        }

        if ((report == 2) && !errors) cat('No errors found :-)\n')

        out <- list(errors = errors,
            trapNAOK = trapNAOK,
            trapcovariatesOK = trapcovariatesOK,
            usagedetectorsOK = usagedetectorsOK,
            usagenonzeroOK = usagenonzeroOK,
            markoccOK == markoccOK,
            areaOK = areaOK,
            polyIDOK = polyIDOK,
            polyconvexOK = polyconvexOK,
            usagecount = usagecount
        )

        invisible(out)

    }
}
############################################################################################

verify.capthist <- function (object, report = 2, tol = 0.01, ...) {

## Check internal consistency of 'capthist' object
##
## -- 'traps' component present
## -- verify(traps)
## -- No live releases
## -- Live detection(s) after reported dead
## -- More than one capture in single-catch trap(s)

## -- Number of rows in 'traps' object not compatible with reported detections
## -- Number of rows in dataframe of individual covariates differs from capthist
## -- Number of occasions in usage matrix differs from capthist
## -- Detections at unused detectors
    
## -- resighting attributes Tu, Tm compatible if present
## -- no resightings on marking occasions, or new animals on resighting occasions
    
    if (!inherits(object, 'capthist'))
        stop ("object must be of class 'capthist'")
    if (inherits(object, 'list')) {
        temp <- lapply (object, verify, report = min(report, 1))
        anyerrors <- any(sapply(temp, function(x) x$errors))

        ## check covariate factor levels conform across sessions
        if (!all(sapply(covariates(object), is.null))) {
            if (any(sapply(covariates(object), is.null)))
                stop ("mixture of NULL and non-NULL individual covariates in different sessions")            
            covariatelevelsOK <- checkcovariatelevels(covariates(object))
            if (!covariatelevelsOK & report>0) {
                warning ('Levels of factor covariate(s) differ between sessions')
            }
        }
        if ((report == 2) && !anyerrors)
            cat('No errors found :-)\n')
        invisible(list(errors = anyerrors, bysession = temp))
    }
    else {
        ## preliminaries
        dim3 <- length(dim(object)) == 3
        count <- detector(traps(object)) %in% .localstuff$countdetectors
        area <- FALSE
        binary <- detector(traps(object)) %in% c('proximity')
        single <- detector(traps(object)) %in% c('single')
        signal <- detector(traps(object)) %in% c('signal','signalnoise')
        poly <- detector(traps(object)) %in% c('polygon', 'polygonX')
        telem <- detector(traps(object)) %in% c('telemetry')
        transect <- detector(traps(object)) %in% c('transect', 'transectX')
        unmarked <- detector(traps(object)) %in% c('unmarked')
        markocc <- markocc(traps(object))
        sighting <- sighting(traps(object))
        allsighting <- if (sighting) !any(markocc>0) else FALSE

        NAOK <- TRUE
        deadOK <- TRUE
        usageOK <- TRUE
        usageoccasionsOK <- TRUE
        usagedetectorsOK <- TRUE
        usagenonzeroOK <- TRUE
        detectornumberOK <- TRUE
        detectorconflcts <- NULL
        singleOK <- TRUE
        binaryOK <- TRUE
        countOK <- TRUE
        cutvalOK <- TRUE
        signalOK <- TRUE
        xyOK <- TRUE
        xyinpolyOK <- TRUE
        xyontransectOK <- TRUE
        IDOK <- TRUE
        rownamesOK <- TRUE
        xylistOK <- TRUE
        sightingsOK <- TRUE
        sightingusageOK <- TRUE        
        MOK <- TRUE
        ROK <- TRUE
        
        if (!is.null(covariates(object)))
            if ((ncol(covariates(object)) == 0 ) |
                (nrow(covariates(object)) == 0 ))
                covariates(object) <- NULL

        ## 1
        trapspresentOK <- !is.null(traps(object))

        ## standalone check of detectors
        ## this is done one session at a time
        ## so does not check between session agreement of covariates
        if (trapspresentOK)
            trapcheck <- verify(traps(object), report = 0)  ## delay reporting
        else
            trapcheck <- list(errors=TRUE)

        ## 2
        trapsOK <- !trapcheck$errors

        ## 3
        if (length(object)==0) {
            detectionsOK <- unmarked   ## and presence? 2011-10-01
        }
        else  {

            detectionsOK <- sum(object[object>0]) > 0

            ## 4
            NAOK <- !any(is.na(object))

            ## 5
            if (signal) {
                ## must have cutval; deads not allowed
                if (length(attr(object,'cutval')) != 1) cutvalOK <- FALSE
                else {
                    maxbyanimal <- tapply(signal(object), animalID(object), max)
                    maxbyanimal <- maxbyanimal[!is.na(maxbyanimal)]  ## because NA OK 2012-02-10
                    if (any(maxbyanimal < attr(object,'cutval')))
                        cutvalOK <- FALSE
                }
                if (length(signal(object)) != sum(abs(object)))
                    signalOK <- FALSE
            }
            ## 6
            else {
                fn <- function(x) {
                    if (dim3) x <- apply(x,1,min)
                    (min(x)<0) && (tail(x[x!=0],1)>0)
                }
                undead <- apply(object, 1, fn)
                deadOK <- !any(undead)
                if (!deadOK) {
                    if (dim3)
                        reincarnated <- object[undead,,, drop=F]
                    else
                        reincarnated <- object[undead,, drop=F]
                }
            }

            ## 7
            if (single) {
                fn <- function (x) duplicated(abs(x)[x!=0])
                multiple <- apply(object, 2, fn)
                singleOK <- !any(unlist(multiple))
            }

            ## 8
            if (binary) {
                ## must be binary
                multiples <- sum(abs(object)>1)
                binaryOK <- multiples == 0
            }

            ## 9
## blocked 2010-12-01 - no problem with 'dead' count
##            if (count) {
##                countOK <- all (object>=0)
##            }
        }

        ## 10
        if (nrow(object) > 0) {
            if (poly | transect | telem) {
                detectornumberOK <- ifelse (dim3,
                    length(levels(polyID(traps(object)))) == dim(object)[3],
                    max(abs(object)) <= ndetector(traps(object)))
            }
            else
                detectornumberOK <- ifelse (dim3,
                  dim(object)[3] == nrow(traps(object)),
                  max(abs(object)) <= nrow(traps(object)))
        }

        ## 11
        covariatesOK <- ifelse(is.null(covariates(object)),
            TRUE,
            nrow(covariates(object)) == nrow(object))

        ## is 'usage' of traps consistent with reported detections?
        if (!is.null(usage(traps(object)))) {
            conflcts <- 0

            ## 12
            usageoccasionsOK <- ncol(usage(traps(object))) == ncol(object)

            if (detectionsOK) {
                # 2012-12-17
                # notused <- !usage(traps(object))   ## traps x occasions
                notused <- usage(traps(object)) == 0 ## traps x occasions
                if (dim3) {
                    if (usagedetectorsOK && usageoccasionsOK) {
                        tempobj <- aperm(object, c(2,3,1))   ## occasion, traps, animal sKn
                        # 2012-12-17
                        # tempuse <- array(t(usage(traps(object))), dim=dim(tempobj))
                        tempuse <- array(t(usage(traps(object))>0), dim=dim(tempobj)) # repl to fill
                        conflcts <- (abs(tempobj)>0) && (tempuse==0)
                        tempobjmat <- array(tempobj[,,1], dim= dim(tempobj)[1:2])
                        occasion <- rep(row(tempobjmat), dim(tempobj)[3])
                        detector <- rep(col(tempobjmat), dim(tempobj)[3])
                        ID <- rep(rownames(object), rep(prod(dim(tempobj)[1:2]), nrow(object)))
                        detectorconflcts <- as.data.frame(cbind(ID,detector,occasion)[conflcts,])
                    }
                }
                else {
                    if (usagedetectorsOK && usageoccasionsOK) {
                        occasion <- occasion(object)
                        ID <- animalID(object, names = FALSE)
                        detector <- trap(object, names = FALSE)
                        conflcts <- notused[cbind(detector, occasion)] > 0
                        detectorconflcts <- as.data.frame(cbind(ID,detector,occasion)[conflcts,])
                    }
                }
            }

            ## 13
            usageOK <- sum(conflcts)==0

        }

        ## 14
        if (poly | telem) {
            xy <- xy(object)
            if (detector(traps(object)) %in% c('polygon','telemetry'))
                xyOK <- nrow(xy) == sum(abs(object))
            else
                xyOK <- nrow(xy) == sum(abs(object)>0)
            inpoly <- xyinpoly(xy(object), traps(object))

            inpoly <- inpoly == trap(object, names = F)
            xyinpolyOK <- all(inpoly)
        }
        if (transect) {
            xy <- xy(object)
            ID <- as.numeric(animalID(object))   ## does this allow for alpha names?
            if (detector(traps(object))=='transect')
                xyOK <- nrow(xy) == sum(abs(object))
            else
                xyOK <- nrow(xy) == sum(abs(object)>0)
            ontransect <- xyontransect(xy(object), traps(object), tol = tol)
            ontransect <- ontransect == trap(object, names = F)
            xyontransectOK <- all(ontransect)
        }
        if (telem) {
            xyOK <- nrow(xy) == sum(abs(object))
        }

        ## 15
        rownamesOK <- !any(duplicated(rownames(object)))

        ## 16
        ## 2012-10-22
        if (nrow(object)>0) {
            zeros <- apply(abs(object)>0,1,sum)==0
            if (!allsighting) {
                xyl <- telemetryxy(object)
                if (!is.null(xyl) | any(zeros)) {
                    xylistOK <- all(names(xyl) %in% row.names(object))
                    if (!all(row.names(object)[zeros] %in% names(xyl)))
                        xylistOK <- FALSE
                }
            }
        }
        
        ## 17
        ## 2015-10-03
        ## -- resighting attributes Tu, Tm compatible if present
        ## -- no resightings on marking occasions, or new animals on resighting occasions
        
        if (sighting) {

            Tu <- Tu(object)
            Tm <- Tm(object)
            nocc <- ncol(object)
            K <- ndetector(traps(object))
            r <- numeric(nocc)
            usge <- usage(traps(object))
            if (is.null(usge)) usge <- 1
            if (length(markocc) != nocc) sightingsOK <- FALSE
            ## allow scalar summed sighting counts 2015-10-31
            if (!is.null(Tu)) {
                if (length(Tu) > 1) {  
                    if (any((Tu>0) & (usge==0))) sightingusageOK <- FALSE
                    if (ncol(Tu) != nocc) sightingsOK <- FALSE
                    if (nrow(Tu) != K) sightingsOK <- FALSE
                    r <- r + apply(Tu,2,sum)
                }
                if (any(Tu<0)) sightingsOK <- FALSE
                if (!all(is.wholenumber(Tu))) sightingsOK <- FALSE
            }
            if (!is.null(Tm)) {
                if (length(Tm) > 1) {
                    if (any((Tu>0) & (usge==0))) sightingusageOK <- FALSE
                    if (nrow(Tm) != K) sightingsOK <- FALSE
                    if (ncol(Tm) != nocc) sightingsOK <- FALSE
                    r <- r + apply(Tm,2,sum)
                }
                if (any(Tm<0)) sightingsOK <- FALSE
                if (!all(is.wholenumber(Tm))) sightingsOK <- FALSE
            }
            if (sightingsOK) {
                if (length(Tu)>1) {
                    u <- unlist(counts(object, 'u'))[1:nocc]
                    if ( any((u > 0) & (markocc<1) & !allsighting) ) MOK <- FALSE
                }
                if (length(Tm)>1)
                    if ( any((r > 0) & (markocc>0))) ROK <- FALSE
            }
        }

        errors <- !all(c(trapspresentOK, trapsOK, detectionsOK, NAOK,
            deadOK, singleOK, binaryOK, countOK, cutvalOK, signalOK,
            detectornumberOK, covariatesOK, usageoccasionsOK, usageOK,
            xyOK, xyinpolyOK, xyontransectOK, IDOK, rownamesOK, xylistOK,
            sightingsOK, sightingusageOK, MOK, ROK))

        if (report > 0) {
            if (errors) {
                cat ('Session', session(object), '\n')

                if (!trapspresentOK) {
                    cat ('No valid detectors\n')
                }
                if (!trapsOK) {
                    cat ('Errors in traps\n')
                    if (!trapcheck$trapNAOK) {
                        cat ('Missing detector coordinates not allowed\n')
                    }
                    if (!trapcheck$trapcovariatesOK) {
                        cat ('Wrong number of rows in dataframe of detector covariates\n')
                        cat ('traps(capthist) :', ndetector(traps(object)), 'detectors\n')
                        cat ('covariates(traps(capthist)) :', nrow(covariates(traps(object))), 'detectors\n')
                    }
                    if (!trapcheck$usagedetectorsOK) {
                        cat ('Conflicting number of detectors in usage matrix\n')
                        cat ('traps(capthist) :', ndetector(traps(object)), 'detectors\n')
                        cat ('usage(traps(capthist)) :', nrow(usage(traps(object))), 'detectors\n')
                    }
                    if (!trapcheck$usagenonzeroOK) {
                        cat ("Occasions when no detectors 'used'\n")
                        cat ((1:length(trapcheck$usagecount))[trapcheck$usagecount==0], '\n')
                    }
                }

                if (!detectionsOK) {
                    cat ('No live releases\n')
                }

                if (!NAOK) {
                    cat ('Missing values not allowed in capthist\n')
                }

                if (!deadOK) {
                    cat ('Recorded alive after dead\n')
                    print(reincarnated)
                }

                if (!singleOK) {
                    cat ('More than one capture in single-catch trap(s)\n')
                }

                if (!binaryOK) {
                    cat ('More than one detection per detector per occasion at proximity detector(s)\n')
                }

                if (!countOK) {
                    cat ('Count(s) less than zero\n')
                }

                if (!cutvalOK) {
                    cat ('Signal less than cutval or invalid cutval\n')
                }

                if (!signalOK) {
                    cat ('Signal attribute does not match detections\n')
                }

                if (!detectornumberOK) {
                    cat ('traps object incompatible with reported detections\n')
                    cat ('traps(capthist) :', ndetector(traps(object)), 'detectors\n')
                    if (dim3)
                        cat ('capthist :', dim(object)[3], 'detectors\n')
                    else
                        cat ('capthist :', max(abs(object)), 'max(detector)\n')
                }

                if (!covariatesOK) {
                    cat ('Wrong number of rows in dataframe of individual covariates\n')
                    cat ('capthist :', nrow(object), 'individuals\n')
                    cat ('covariates(capthist) :', nrow(covariates(object)), 'individuals\n')
                }
                if (!usageoccasionsOK) {
                    cat ('Conflicting number of occasions in usage matrix\n')
                    cat ('capthist :', ncol(object), 'occasions\n')
                    cat ('usage(traps(capthist)) :', ncol(usage(traps(object))), 'occasions\n')
                }
                if (!usageOK) {
                    cat ("Detections at 'unused' detectors\n")
                    print(detectorconflcts)
                }
                if (!xyOK) {
                    cat ("Polygon or telemetry detector xy coordinates of detections",
                         " do not match counts\n")
                }
                if (!xyinpolyOK) {
                    cat ("XY coordinates not in polygon\n")
                    print (xy(object)[!inpoly,])
                }
                if (!xyontransectOK) {
                    cat ("XY coordinates not on transect\n")
                    print (xy(object)[!ontransect,])
                }
                if (!IDOK) {
                    cat ("Polygon detector mismatch between ID attribute and counts\n")
                }
                if (!rownamesOK) {
                    cat ("Duplicated row names (animal ID)\n")
                }
                if (!xylistOK) {
                    cat ("Telemetry data (xylist) do not match capture histories\n")
                }
                if (!sightingsOK) {
                    cat("Incompatible dimensions of sighting attributes markocc, Tu or Tm\n")
                }
                if (!sightingusageOK) {
                    cat("Sightings at unused detectors\n")
                    Tu <- Tu(object)
                    Tm <- Tm(object)
                    bad <- (Tu>0) & (usage(traps(object))==0) & (length(Tu)>1)
                    if (sum(bad)>0) {
                        cat("Tu\n")
                        print(cbind(Detector = row(bad)[bad], Occasion = col(bad)[bad]))
                    }
                    bad <- (Tm>0) & (usage(traps(object))==0) & (length(Tm)>1)
                    if (sum(bad)>0) {
                        cat("Tm\n")
                        print(cbind(Detector = row(bad)[bad], Occasion = col(bad)[bad]))
                    }
                }
                if (!MOK) {
                    cat("New individual(s) on sighting-only occasion\n")
                    occ <- split(occasion(object), animalID(object, names=TRUE))
                    firstocc <- sapply(occ, '[', 1)
                    print(firstocc[firstocc %in% (1:ncol(object))[markocc(traps(object))<1]])
                }
                if (!ROK) {
                    cat("Sighting(s) on marking-only occasion\n")
                }
            }

            if ((report == 2) && !errors) cat('No errors found :-)\n')

        }

        out <- list(errors = errors, trapcheck = trapcheck)
        if (!is.null(detectorconflcts)) out$detections.at.unused.detectors <- detectorconflcts
        invisible(out)
    }
}
############################################################################################

verify.mask <- function (object, report = 2, ...) {

## Check internal consistency of 'mask' object
##
## valid x and y coordinates
## nrow(covariates) = nrow(object)
## ...also look at attributes?

    if (!inherits(object, 'mask'))
        stop ("object must be of class 'mask'")

    if (inherits(object, 'list')) {
        temp <- lapply (object, verify, report = min(report, 1))
        anyerrors <- any(sapply(temp, function(x) x$errors))

        ## check covariate factor levels conform across sessions
        if (!all(sapply(covariates(object), is.null))) {
            covariatelevelsOK <- checkcovariatelevels(covariates(object))
            if (!covariatelevelsOK & report>0) {
                warning ('Levels of factor mask covariate(s) differ between sessions')
            }
        }
        if ((report == 2) && !anyerrors)
            cat('No errors found :-)\n')
        invisible(list(errors = anyerrors, bysession = temp))
    }
    else {

        ## 1
        xyOK <- !(is.null(object$x) | is.null(object$y) | any(is.na(object)))
        xyOK <- xyOK && is.numeric(unlist(object))

        ## 2

        if (!is.null(covariates(object)))
            covariatesOK <- ifelse (nrow(covariates(object))>0,
            nrow(object) == nrow(covariates(object)), TRUE)
        else
            covariatesOK <- TRUE

        errors <- !all(c(xyOK, covariatesOK))

        if (report > 0) {
            if (errors) {
                ## cat ('Session', session(object), '\n')

                if (!xyOK) {
                    cat ('Invalid x or y coordinates in mask\n')
                }

                if (!covariatesOK) {
                    cat ('Number of rows in covariates(mask) differs from expected\n')
                }
            }

            if ((report == 2) && !errors) cat('No errors found :-)\n')
        }

        out <- list(errors = errors)
        invisible(out)
    }
}
############################################################################################

