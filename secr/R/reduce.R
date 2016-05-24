###########################################################################################
## package 'secr'
## reduce.R
## 2010-12-01 check for overlapping columns
## 2011-02-08 rewritten to include polygonX and transectX detectors, and simplified
## 2011-03-18 output to unmarked
## 2011-03-21 'by' argument
## 2012-12-17 non-binary usage
## 2012-12-21 re-write reduce.capthist to include spatial lumping
## 2012-12-21 amalgamated 'reduce' methods
## 2013-11-20 telemetry allowed
############################################################################################

#----------------------------------------------------------------------------------------------------
# Dimensions   2       2      3	         3      2        2          3       3+      3+
#----------------------------------------------------------------------------------------------------
#              single  multi  proximity  count  polygonX transectX  signal  polygon transect
# single	&#	&	*	   *      NA       NA         NA      NA      NA
# multi		&#	&	*	   *      NA       NA         NA      NA      NA
# proximity	&#	&	*	   *      NA       NA         NA      NA      NA
# count		&#@	&@	@	   *      NA       NA         NA      NA      NA
# polygonX	&#	&	*	   *      &$       NA         NA      NA      NA
# transectX	&#	&	*	   *      NA       &$         NA      NA      NA
# signal        &#	&	*	   @      @        @          $       NA      NA
# polygon     	&#@~	&@~	@~	   *~     @        NA         NA       *      NA
# transect     	&#@~	&@~	@~	   *~     NA       @          NA      NA       *

#----------------------------------------------------------------------------------------------------
#  * no loss of data
#  # must choose among animals (more than one animal)
#  & must choose among traps (animal caught more than once)
#  @ reduce to binary 
#  $ apply rule for combining signals or locations (first, last, random, min, max, mean)
#  ~ form new point detectors from mean of vertices (assumes symmetry)
#  NA not feasible
############################################################################################


reduce     <- function (object, ...) UseMethod("reduce")

reduce.default <- function (object, columns, ...) {
  object <- as.matrix(object)
  if (any(is.na(object)))
      warning ("NAs in input converted to zero")
  firsttrap <- function (y) y[abs(y)>0][1]    # first non-zero
  fnmulti   <- function (occ) apply (object[,occ,drop=F], 1, firsttrap)
  nrow <- nrow(object)
  nnew <- length(columns)
  temp <- sapply (columns, fnmulti)
  temp[is.na(temp)] <- 0
  temp
}

reduce.traps <- function (object, newtraps = NULL, span = NULL, rename = FALSE, ...) {
    if (ms(object)) {
        lapply(object, reduce, newtraps = newtraps, span = span, rename = rename, ...)
    }
    else {
        if (!inherits(object, 'traps'))
            stop ("requires traps object")
        ## allow for vector input, or distance threshold
        if (!is.null(span))
            newtraps <- cutree (hclust(dist(object), ...), h = span)
        if (missing(newtraps))
            newtraps <- as.list(1:ndetector(object))  ## no change
        if (!is.list(newtraps)) {
            if (is.null(names(newtraps)))
                names(newtraps) <- 1:length(newtraps)
            newtraps <- split(names(newtraps), newtraps)
        }
        if (!detector(object) %in% .localstuff$pointdetectors)
            stop ("reduce.traps is only for point detectors")
        if (any(duplicated(unlist(newtraps))))
            stop("traps should not appear in more than one group")

        if (any(sapply(newtraps, is.character))) {
            ## convert to numeric trap indices
            newtraps <- lapply(newtraps, function(x) match(x, rownames(object)))
        }

        nnew <- length(newtraps)
        if (rename)
            newtrapnames <- 1:nnew
        else if (any(sapply(newtraps, is.character))) {
            newtrapnames <- sapply(newtraps, function(x) paste(x, collapse='+'))
        }
        else {
            namefn <- function(x) paste(rownames(object)[x], collapse='+')
            newtrapnames <- sapply(newtraps, namefn)
        }

        g <- rep(1:nnew, sapply(newtraps,length))
        splitfactor <- numeric(ndetector(object))   ## initially zero
        splitfactor[unlist(newtraps)] <- g             ## levels are indices of groups 1, 2,...
        splitfactor <- factor(splitfactor)

        grouped <- split.data.frame(object, splitfactor) ## otherwise calls split.traps
        ## or could use for usage and cov below and save repeat splitting...
        if (any (splitfactor == 0)) grouped <- grouped[-1]  ## drop unwanted traps
        names(grouped) <- newtrapnames
        grouped <- lapply(grouped, function(df) apply(df,2,mean))
        newxy <- do.call(rbind, grouped)
        newxy <- as.data.frame(newxy)
        class (newxy)   <- c('traps', 'data.frame')
        detector(newxy) <- detector(object)
        markocc(newxy) <- markocc(object)

        if (!is.null(usage(object))) {
            usagelist <- split(as.data.frame(usage(object)), splitfactor)
            if (any (splitfactor == 0)) usagelist <- usagelist[-1]  ## drop unwanted traps
            temp <- lapply(usagelist, function(x) apply(x,2,sum))
            temp <- do.call(rbind, temp)
            temp <- as.matrix(temp)
            dimnames(temp) <- list(newtrapnames, 1:ncol(temp))
            usage(newxy) <- temp

            daily <- as.data.frame(temp)
            names(daily) <- paste('occ', names(daily), sep='')
        }

        if (!is.null(covariates(object))) {
            covlist <- split(covariates(object), splitfactor)
            if (any (splitfactor == 0)) covlist <- covlist[-1]  ## drop unwanted traps
            varying <- sapply(covlist, function(x) apply(x,2,function(y) length(unique(y))>1))
            if (any(varying))
                warning("covariates vary within groups; using only first")
            temp <- lapply(covlist, head, 1)
            temp <- do.call(rbind, temp)
            covariates(newxy) <- temp
            if (!is.null(usage(object)))
                covariates(newxy) <- cbind(covariates(newxy), daily)
            else
                covariates(newxy)$combined <- sapply(newtraps, length)
        }
        else {  ## new covariate dataframe
            if (!is.null(usage(object)))
                covariates(newxy) <- daily
            else
                covariates(newxy) <- data.frame(combined = sapply(newtraps, length))
        }

        sp <- spacing(newxy, recalculate = TRUE)
        if (!is.null(sp)) spacing(newxy) <- sp
        temp <- as.numeric(levels(splitfactor)[splitfactor])
        temp[temp==0] <- NA
        attr(newxy, 'newtrap') <- temp
        newxy
    }
}

poly2point <- function (object, detector = 'count') {
    if (!detector(object) %in% c('polygon','polygonX'))
        stop ("requires 'polygon' input")
    if (detector %in% .localstuff$polydetectors)
        stop ("requires non-polygon, non-transect output")
    temp <- split(object, polyID(object))
    temp <- lapply(temp, function(df) apply(df,2,mean))
    temp1 <- t(abind(temp, along=2))
    dimnames(temp1) <- list(levels(polyID(object)), c('x','y'))
    temp <- data.frame(temp1, row.names=NULL)
    class (temp)   <- c('traps', 'data.frame')
    detector(temp) <- detector
    usage(temp)    <- usage(object)
    covariates(temp) <- covariates(object)
    attr(temp,'spacex') <- 100 * (searcharea(object)/nrow(temp))^0.5
    attr(temp,'spacey') <- attr(temp,'spacex')
    temp
}
#-----------------------------------------------------------------------------

transect2point <- function (object, detector = 'count') {
    if (!detector(object) %in% c('transect','transectX'))
        stop ("requires 'transect' input")
    if (detector %in% c('transect', 'transectX'))
        stop ("requires non-transect output")
    temp <- split(object, transectID(object))
    temp <- lapply(temp, function(df) apply(df,2,mean))
    temp1 <- t(abind(temp, along=2))
    dimnames(temp1) <- list(levels(transectID(object)), c('x','y'))
    temp <- data.frame(temp1, row.names=NULL)
    class (temp)   <- c('traps', 'data.frame')
    detector(temp) <- detector
    usage(temp)    <- usage(object)
    covariates(temp) <- covariates(object)
    attr(temp,'spacex') <- mean(transectlength(object))/2   ## arbitrary
    attr(temp,'spacey') <- attr(temp,'spacex')
    temp
}
#-----------------------------------------------------------------------------

## function to make list in which each component is a
## subset of occasions (for use in reduce.capthist)
## MGE 2011-03-10

split.by <- function (x, by) {
    if ((length(x) == 1) & (x[1] > 1))
        x <- 1:x
    if (by < 1)
        stop ("invalid 'by' argument")
    index <- 1:length(x)
    gp <- trunc((index-1)/by) + 1
    split (index, gp)
}
#----------------------------------------------------------------------------------------------

reduce.capthist <- function (object, newtraps = NULL, span = NULL,
    rename = FALSE, newoccasions = NULL, by = 1, outputdetector =
    detector(traps(object)), select = 'last', dropunused = TRUE,
    verify = TRUE, sessions = NULL, ...) {

    # newoccasions - list, each component gives occasions to include in new capthist
    # newtraps     - list, each component gives traps to include in new capthist

    #--------------------------------
    seltrap <- function (y) {
        y <- t(y)                  # allow for occasion x trap matrix in proximity, count data
        y <- y[abs(y)>0]
        if (length(y)<1) y <- 0
        if (length(y) == 1) y
        else switch (select,
            first = head(y,1),    # first non-null
            last = tail(y,1),     # last non-null
            random = sample (size=1, y)    # random non-null, weighted by frequency in sample
        )
    }
    #--------------------------------
    selused <- function (y) {
        y <- y[abs(y)>0]
        if (length(y)<1) y <- 0
        if (length(y) == 1) y
        else switch (select,
            first = head(y,1),    # first non-null
            last = tail(y,1),     # last non-null
            random = sample (size=1, y)    # random non-null, weighted by frequency in sample
        )
    }
    #----------------------------------------------------------------------------
    # functions applied to collapse a set of occasions 'occ' to a single occasion
    # result is a vector (single, multi detectors)
    fnused <- function (occ, fn) {
        if (length(occ)>0) {
            temp <- usage(trps)[,occ,drop=F]
            apply (temp, 1, fn)
        }
        else NULL
    }
    #----------------------------------------------------------------------------
    collapse <- function (df) {
        ## reduce data frame to a single row
        if (nrow(df)>1) {
            df$alive <- rep(all(df$alive), nrow(df))
            allowedCriteria <- c('first','last','random')
            if (!(select %in% allowedCriteria))
                stop ("selection criterion for signal should be one of ",
                      paste(sapply(allowedCriteria, dQuote),collapse=','))
            index <- switch (select, first = 1, last = nrow(df),
                random = sample.int (nrow(df),1) )
            df <- df[index,,drop=FALSE]
        }
        df
    }
    #----------------------------------------------------------------------------

    # main line
    if (ms(object)) {
        if (is.null(sessions)) sessions <- 1:length(object)
        temp <- lapply (object[sessions], reduce,
            newoccasions = newoccasions,
            by = by,
            newtraps = newtraps,
            span = span,
            outputdetector = outputdetector,
            select = select,
            rename = rename,
            dropunused = dropunused,
            verify = verify,
            ...)
        class(temp) <- c('list', 'capthist')
        if (length(temp) == 1) temp <- temp[[1]]
        return(temp)
    }
    else {
        if (tolower(by) == 'all')
            by <- ncol(object)
        polygons <- c('polygon','polygonX')
        transects <- c('transect','transectX')
        inputdetector <- detector(traps(object))
        ntrap <- ndetector(traps(object))  ## npoly if 'polygon' or 'transect'

        nrw <- nrow(object)
        if (is.null(newoccasions)) {
            newoccasions <- split.by (1:ncol(object), by)
            if ((ncol(object) %% by) > 0)
                warning ("number of occasions is not a multiple of 'by'")
        }

        if (is.null(outputdetector)) outputdetector <- inputdetector
        if (!(outputdetector %in% .localstuff$validdetectors))
            stop ("'outputdetector' should be one of ",
                  paste(sapply(.localstuff$validdetectors, dQuote),collapse=','))
        if ((!(inputdetector %in% c('signal','signalnoise'))) && (outputdetector == 'signal'))
                stop ("cannot convert non-signal data to signal data")
        if ((!(inputdetector %in% c('signalnoise'))) && (outputdetector == 'signalnoise'))
                stop ("cannot convert non-signalnoise data to signalnoise data")
        if ((!(inputdetector %in% polygons)) && (outputdetector %in% polygons))
                stop ("cannot convert non-polygon data to 'polygon' data")
        if ((!(inputdetector %in% transects)) && (outputdetector %in% transects))
                stop ("cannot convert non-transect data to 'transect' data")

        ######################################
        ## add usage to count pooled occasions
        ## 2014-10-11
        if (outputdetector %in% .localstuff$countdetectors) {
            if (is.null(usage(traps(object)))) {
                usage(traps(object)) <- matrix(1, nrow = ndetector(traps(object)),
                                               ncol = ncol(object))
            }
        }

        ######################################
        ## check newoccasions
        for (i in length(newoccasions):1) {
            occ <- newoccasions[[i]]
            occ <- occ[occ %in% (1:ncol(object))]  ## discard nonexistent occ
            if (length(occ)==0)
                newoccasions[[i]] <- NULL
            else
                newoccasions[[i]] <- occ
        }
        cumocc <- numeric(0)
        for (i in length(newoccasions):1) {
            if (any (newoccasions[[i]] %in% cumocc))
                warning ("new occasions overlap")
            cumocc <- c(cumocc, newoccasions[[i]])
        }
        nnew <- length(newoccasions)
        newcols <- rep(1:nnew, sapply(newoccasions,length))
        newcols <- factor(newcols)

        ####################################
        ## check and build newtraps

        trps <- traps(object)
        reducetraps <- !is.null(newtraps) | !is.null(span)
        if (reducetraps) {
            trps <- reduce(trps, newtraps = newtraps, span =
                span, rename = rename, ...)
            newtrapID <- attr(trps, 'newtrap')
            ntrap <- ndetector(trps)
        }

        ####################################
        ## build dataframe of observations

        df <- data.frame(
            trap = trap(object, names = F),
            occ = occasion(object),
            ID = animalID(object, names = F),
            alive = alive(object))
        if (reducetraps)
            df$trap <- newtrapID[df$trap]
        if (outputdetector %in% c(polygons, transects, 'telemetry')) {   ## telemetry 2013-11-20
            df$x <- xy(object)[,1]
            df$y <- xy(object)[,2]
        }
        if (!is.null(attr(object,'signalframe'))) {
            df <- cbind(df, attr(object,'signalframe'))
        }
        df$newocc <- newcols[match(df$occ, unlist(newoccasions))]
        if (dropunused) {
            df$newocc <- factor(df$newocc)
            nnew <- length(levels(df$newocc))
        }
        df <- df[!is.na(df$newocc),]                   ## drop null obs
        df$newID <- factor(df$ID)                      ## assign newID
        if (outputdetector %in% .localstuff$exclusivedetectors) {
            ID.occ <- interaction(df$ID, df$newocc)
            dflist <- split(df, ID.occ)
            dflist <- lapply(dflist, collapse)
            df <- do.call(rbind, dflist)
        }
        
        if (outputdetector %in% c('single')) {
            occ.trap <- interaction(df$newocc,df$trap)
            dflist <- split(df, occ.trap)
            dflist <- lapply(dflist, collapse)
            df <- do.call(rbind, dflist)
        }
        df$newID <- factor(df$ID)                     ## re-assign newID

        ####################################
        ## build new object
        validrows <- (1:nrow(object)) %in% df$ID   ## or newID??? 2012-12-12

        if (outputdetector %in% .localstuff$exclusivedetectors) {
            alivesign <- df$alive*2 - 1
            tempnew <- matrix(0, nrow = sum(validrows), ncol = nnew)
            tempnew[cbind(df$newID, df$newocc)] <- df$trap * alivesign
        }
        else {
            df$trap <- factor(df$trap, levels=1:ntrap)
            tempnew <- table(df$newID, df$newocc, df$trap)
            alivesign <- tapply(df$alive, list(df$newID,df$newocc,df$trap),all)
            alivesign[is.na(alivesign)] <- TRUE
            alivesign <- alivesign * 2 - 1
            if (! (outputdetector %in% .localstuff$countdetectors)
                && (length(tempnew)>0)) {
                ## convert 'proximity' and 'signal' to binary
                tempnew[tempnew>0] <- 1
            }
            tempnew <- tempnew * alivesign
        }

        ################################
        ## general attributes
        class(tempnew) <- 'capthist'
        session(tempnew) <- session(object)
        attr(tempnew, 'n.mash') <- attr(object, 'n.mash')
        attr(tempnew, 'centres') <- attr(object, 'centres')
        if ((inputdetector %in% polygons) && !(outputdetector %in% polygons))
            traps(tempnew) <- poly2point(trps)
        else
            if ((inputdetector %in% transects) && !(outputdetector %in% transects))
                traps(tempnew) <- transect2point(trps)
        else
            traps(tempnew) <- trps
        detector(traps(tempnew)) <- outputdetector

        ################################
        ## covariates and ancillary data

        if (!is.null(covariates(object)))
             covariates(tempnew) <- covariates(object)[validrows,,drop=F]
     
        detectorder <- order(df$trap, df$newocc,df$ID)  ## CHECK!
        if (outputdetector %in% c(polygons, transects, 'telemetry'))  ## telemetry added 2013-11-20
            xy(tempnew) <- df[detectorder,c('x','y'),drop=FALSE]
        if (outputdetector %in% c('signal','signalnoise')) {
            sigcolumns <- names(attr(object,'signalframe'))
            attr(tempnew,'signalframe') <- df[detectorder, sigcolumns, drop=FALSE]
            attr(tempnew, 'cutval') <- attr(object, 'cutval')
        }
        ##################################
        ## fix usage for reduced occasions
        if (nrow(tempnew) > 0)
            dimnames(tempnew)[[1]] <- 1:nrow(tempnew)  ## temporary, for animalID in subset
        if (!is.null(usage(traps(tempnew)))) {
            usagematrix <- unlist(sapply (newoccasions, fnused, sum))
            usagematrix <- matrix(usagematrix, nrow = ndetector(traps(tempnew)))
            usage(traps(tempnew)) <- usagematrix

            if (dropunused) {
                OK <- apply(usagematrix, 1, sum) > 0
                tempnew <- subset(tempnew, traps = OK)
            }
        }
        tempnew[is.na(tempnew)] <- 0
        ################################
        ## dimnames
        if (nrow(tempnew) > 0) {
            indices <- (1:length(validrows))[validrows]
            rowname <- rownames(object)[indices]
        }
        else
            rowname <- NULL
        if (length(dim(tempnew)) == 3)
            dimnames(tempnew) <- list(rowname,1:nnew,NULL)   # renew numbering
        else
            dimnames(tempnew) <- list(rowname,1:nnew)

        if (verify) verify(tempnew, report=1)
        tempnew
    }
}
############################################################################################
