
# ------------------------ VALIDATION FUNCTIONS --------------------------

validunmarkedFrame <- function(object) {
    errors <- character(0)
    M <- nrow(object@y)
    J <- ncol(object@y)
    if(!is.null(object@siteCovs))
        if(nrow(object@siteCovs) != M)
            errors <- c(errors,
               "siteCovData does not have same size number of sites as y.")
    if(!is.null(obsCovs(object)) & !is.null(obsNum(object)))
        if(nrow(object@obsCovs) != M*obsNum(object))
            errors <- c(errors, "obsCovData does not have M*obsNum rows.")
    if(length(errors) == 0)
        TRUE
    else
        errors
}

# --------------------------- DATA CLASSES -------------------------------

# Class to hold data for analyses in unmarked.
setClass("unmarkedFrame",
    representation(y = "matrix",
        obsCovs = "optionalDataFrame",
        siteCovs = "optionalDataFrame",
        mapInfo = "optionalMapInfo",
        obsToY = "optionalMatrix"),
    validity = validunmarkedFrame)

## a class for multi-season data

setClass("unmarkedMultFrame",
    representation(numPrimary = "numeric",
        #data frame in site-major, year-minor order describing siteCovs
        yearlySiteCovs = "optionalDataFrame"),
    contains="unmarkedFrame")



## a class for distance sampling data
setClass("unmarkedFrameDS",
    representation(
        dist.breaks = "numeric",
        tlength = "numeric",
        survey = "character",
        unitsIn = "character"),
    contains = "unmarkedFrame",
    validity = function(object) {
        errors <- character(0)
        J <- numY(object)
        db <- object@dist.breaks
        if(J != length(db) - 1)
            errors <- c(errors, "ncol(y) must equal length(dist.breaks)-1")
        if(db[1] != 0)
            errors <- c(errors, "dist.breaks[1] must equal 0")
        if(!is.null(obsCovs(object)))
            "obsCovs cannot be used with distsamp"
        if(length(errors) == 0) TRUE
        else errors
        })


setClass("unmarkedFrameOccu",
		contains = "unmarkedFrame")

setClass("unmarkedFrameOccuFP",
         representation(
           type = "numeric"),
         contains = "unmarkedFrame")


setClass("unmarkedFramePCount",
		contains = "unmarkedFrame")


setClass("unmarkedFrameMPois",
		representation(
			samplingMethod = "character",
			piFun = "character"),
		contains = "unmarkedFrame")


setClass("unmarkedFrameG3",
         contains = "unmarkedMultFrame")


setClass("unmarkedFramePCO",
         representation(primaryPeriod = "matrix"),
         contains = "unmarkedMultFrame")


setClass("unmarkedFrameGMM",
    representation(
        piFun = "character",
        samplingMethod = "character"),
    contains = "unmarkedFrameG3")

setClass("unmarkedFrameGDS",
    representation(
        dist.breaks = "numeric",
        tlength = "numeric",
        survey = "character",
        unitsIn = "character"),
    contains = "unmarkedFrameG3")

setClass("unmarkedFrameGPC",
    contains = "unmarkedFrameG3")



# ------------------------------- CONSTRUCTORS ---------------------------


# Constructor for unmarkedFrames.
unmarkedFrame <- function(y, siteCovs = NULL, obsCovs = NULL, mapInfo,
                          obsToY) {
    if(!missing(obsToY))
        obsNum <- nrow(obsToY)
    if(class(obsCovs) == "list") {
        obsVars <- names(obsCovs)
        for(i in seq(length(obsVars))) {
            if(!(class(obsCovs[[i]]) %in% c("matrix", "data.frame")))
                stop("At least one element of obsCovs is not a matrix or data frame.")
            if(ncol(obsCovs[[i]]) != obsNum | nrow(obsCovs[[i]]) != nrow(y))
                stop("At least one matrix in obsCovs has incorrect number of dimensions.")
            }
        if(is.null(obsNum)) obsNum <- ncol(obsCovs[[1]]) #??
        obsCovs <- data.frame(lapply(obsCovs, function(x) as.vector(t(x))))
        }
    if(("data.frame" %in% class(y)) | ("cast_matrix" %in% class(y)))
        y <- as.matrix(y)
    if(missing(obsToY)) obsToY <- NULL
    if(missing(mapInfo)) mapInfo <- NULL

    umf <- new("unmarkedFrame", y = y, obsCovs = obsCovs, siteCovs = siteCovs,
        mapInfo = mapInfo, obsToY = obsToY)
    return(umf)
}


unmarkedFrameDS <- function(y, siteCovs = NULL, dist.breaks, tlength,
                            survey, unitsIn, mapInfo = NULL)
{
    if(missing(survey))
        stop("survey argument must be specified")
    if(missing(tlength) & survey == "point")
        tlength <- rep(NA_real_, nrow(y))
    if((survey=="line") & (length(tlength) != nrow(y)))
        stop("tlength should be a vector with length(tlength)==nrow(y)")
    umfds <- new("unmarkedFrameDS", y = y, obsCovs = NULL,
                 siteCovs = siteCovs, dist.breaks = dist.breaks,
                 tlength = tlength, survey = survey, unitsIn = unitsIn,
                 obsToY = matrix(1, 1, ncol(y)))
    return(umfds)
}



unmarkedFrameOccu <- function(y, siteCovs = NULL, obsCovs = NULL, mapInfo)
{
    J <- ncol(y)
    umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J),
                         mapInfo = mapInfo)
    umf <- as(umf, "unmarkedFrameOccu")
    umf
}

unmarkedFrameOccuFP <- function(y, siteCovs = NULL, obsCovs = NULL, type, mapInfo)
{
  J <- ncol(y)
  umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J),
                       mapInfo = mapInfo)
  umf <- as(umf, "unmarkedFrameOccuFP")
  umf@type <- type
  umf
}




unmarkedFramePCount <- function(y, siteCovs = NULL, obsCovs = NULL, mapInfo)
{
    J <- ncol(y)
    umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J),
        mapInfo = mapInfo)
    umf <- as(umf, "unmarkedFramePCount")
    umf
}



unmarkedFrameMPois <- function(y, siteCovs = NULL, obsCovs = NULL, type,
    obsToY, mapInfo, piFun)
{
    if(!missing(type)) {
        switch(type,
            removal = {
                obsToY <- matrix(1, ncol(y), ncol(y))
                obsToY[col(obsToY) < row(obsToY)] <- 0
                #obsToY <- diag(ncol(y))
                #obsToY[upper.tri(obsToY)] <- 1
                piFun <- "removalPiFun"
                },
            double = {
                #obsToY <- matrix(c(1, 0, 0, 1, 1, 1), 2, 3)
                obsToY <- matrix(1, 2, 3)
                piFun <- "doublePiFun"
                })
    } else {
        if(missing(obsToY))
            stop("obsToY is required for multinomial-Poisson data with no specified type.")
        type <- "userDefined"
        }
    umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = obsToY,
        mapInfo = mapInfo)
    umf <- as(umf, "unmarkedFrameMPois")
    umf@piFun <- piFun
    umf@samplingMethod <- type
    umf
}



# This function constructs an unmarkedMultFrame object.
unmarkedMultFrame <- function(y, siteCovs = NULL, obsCovs = NULL,
                              numPrimary, yearlySiteCovs = NULL)
{
    J <- ncol(y)
	  umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J))
    umf <- as(umf, "unmarkedMultFrame")
    umf@numPrimary <- numPrimary

    if(class(yearlySiteCovs) == "list") {
        yearlySiteVars <- names(yearlySiteCovs)
        for(i in seq(length(yearlySiteVars))) {
            if(!(class(yearlySiteCovs[[i]]) %in% c("matrix", "data.frame")))
                stop("At least one element of yearlySiteCovs is not a matrix or data frame.")
            if(ncol(yearlySiteCovs[[i]]) != numPrimary |
                nrow(yearlySiteCovs[[i]]) != nrow(y))
                    stop("At least one matrix in yearlySiteCovs has incorrect number of dimensions.")
            }
        yearlySiteCovs <- data.frame(lapply(yearlySiteCovs, function(x)
            as.vector(t(x))))
        }

    umf@yearlySiteCovs <- yearlySiteCovs
    umf
}





# This function constructs an unmarkedMultFrame object.
unmarkedFrameGMM <- function(y, siteCovs = NULL, obsCovs = NULL, numPrimary,
	yearlySiteCovs = NULL, type, obsToY, piFun)
{
    J <- ncol(y) / numPrimary
    if(!missing(type)) {
      if(!type %in% c("removal", "double"))
        stop("if specifying type, it should either be 'removal' or 'double'")
      switch(type,
        removal = {
          obsToY <- diag(J)
          obsToY[upper.tri(obsToY)] <- 1
          obsToY <- kronecker(diag(numPrimary), obsToY)
          piFun <- "removalPiFun"
          },
        double = {
          obsToY <- matrix(1, 2, 3)
          obsToY <- kronecker(diag(numPrimary), obsToY)
          piFun <- "doublePiFun"
          })
    } else {
        type <- "userDefined"
        if(missing(obsToY))
            stop("obsToY is required for gmultmix data with no specified type.")
        }

    umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = obsToY)
    umf <- as(umf, "unmarkedMultFrame")
    umf@numPrimary <- numPrimary
    if(class(yearlySiteCovs) == "list") {
        yearlySiteVars <- names(yearlySiteCovs)
        for(i in seq(length(yearlySiteVars))) {
            if(!(class(yearlySiteCovs[[i]]) %in% c("matrix","data.frame")))
                stop("At least one element of yearlySiteCovs is not a matrix or data frame.")
            if(ncol(yearlySiteCovs[[i]]) != numPrimary |
                nrow(yearlySiteCovs[[i]]) != nrow(y))
                    stop("At least one matrix in yearlySiteCovs has incorrect number of dimensions.")
            }
        if(is.null(obsNum)) obsNum <- ncol(obsCovs[[1]])
        yearlySiteCovs <- data.frame(lapply(yearlySiteCovs, function(x)
            as.vector(t(x))))
        }
    umf@yearlySiteCovs <- yearlySiteCovs
    umf <- as(umf, "unmarkedFrameGMM")
    umf@piFun <- piFun
    umf@samplingMethod <- type
    umf
}



# This function constructs an unmarkedMultFrame object.
unmarkedFrameGDS <- function(y, siteCovs, numPrimary,
	yearlySiteCovs, dist.breaks, survey, unitsIn, tlength)
{
    J <- ncol(y) / numPrimary
    obsToY <- matrix(1, 1, J)
    obsToY <- kronecker(diag(numPrimary), obsToY)
    if(missing(siteCovs))
        siteCovs <- NULL

    umf <- unmarkedFrame(y = y, siteCovs = siteCovs, obsToY = obsToY)
    umf <- as(umf, "unmarkedMultFrame")
    umf@numPrimary <- numPrimary
    if(missing(yearlySiteCovs))
        yearlySiteCovs <- NULL
    if(class(yearlySiteCovs) == "list") {
        yearlySiteVars <- names(yearlySiteCovs)
        for(i in seq(length(yearlySiteVars))) {
            if(!(class(yearlySiteCovs[[i]]) %in% c("matrix","data.frame")))
                stop("At least one element of yearlySiteCovs is not a matrix or data frame.")
            if(ncol(yearlySiteCovs[[i]]) != numPrimary |
                nrow(yearlySiteCovs[[i]]) != nrow(y))
                    stop("At least one matrix in yearlySiteCovs has incorrect number of dimensions.")
            }
        yearlySiteCovs <- data.frame(lapply(yearlySiteCovs, function(x)
            as.vector(t(x))))
        }
    if(identical(survey, "point")) {
        if(!missing(tlength))
            stop("tlength cannot be specified with point transect data")
        tlength <- rep(1, nrow(y))
        }

    umf@yearlySiteCovs <- yearlySiteCovs
    umf <- as(umf, "unmarkedFrameGDS")
    umf@dist.breaks <- dist.breaks
    umf@survey <- survey
    umf@unitsIn <- unitsIn
    umf@tlength <- tlength
    umf
}



# This function constructs an  object.
unmarkedFrameGPC <- function(y, siteCovs=NULL, obsCovs=NULL, numPrimary,
                             yearlySiteCovs=NULL) {
    if(numPrimary < 2)
        stop("numPrimary must be >1. Use pcount of numPrimary=1")
    J <- ncol(y) / numPrimary
    obsToY <- diag(J*numPrimary)
    if(missing(siteCovs))
        siteCovs <- NULL

    umf <- unmarkedFrame(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
                         obsToY = obsToY)
    umf <- as(umf, "unmarkedMultFrame")
    umf@numPrimary <- numPrimary
    if(missing(yearlySiteCovs))
        yearlySiteCovs <- NULL
    if(class(yearlySiteCovs) == "list") {
        yearlySiteVars <- names(yearlySiteCovs)
        for(i in seq(length(yearlySiteVars))) {
            if(!(class(yearlySiteCovs[[i]]) %in% c("matrix","data.frame")))
                stop("At least one element of yearlySiteCovs is not a matrix or data frame.")
            if(ncol(yearlySiteCovs[[i]]) != numPrimary |
                nrow(yearlySiteCovs[[i]]) != nrow(y))
                    stop("At least one matrix in yearlySiteCovs has incorrect number of dimensions.")
            }
        yearlySiteCovs <- data.frame(lapply(yearlySiteCovs, function(x)
            as.vector(t(x))))
        }

    umf@yearlySiteCovs <- yearlySiteCovs
    umf <- as(umf, "unmarkedFrameGPC")
    umf
}





unmarkedFramePCO <- function(y, siteCovs = NULL, obsCovs = NULL,
    yearlySiteCovs = NULL, mapInfo, numPrimary, primaryPeriod)
{
    M <- nrow(y)
    T <- numPrimary
    J <- ncol(y) / T
    if(missing(primaryPeriod))
        primaryPeriod <- matrix(1:T, M, T, byrow=TRUE)
    if(nrow(primaryPeriod) != M | ncol(primaryPeriod) != T)
        stop("Dimensions of primaryPeriod matrix should be MxT")
    if(any(primaryPeriod < 0, na.rm=TRUE))
        stop("Negative primaryPeriod values are not allowed.")
    if(any(is.na(primaryPeriod)))
        stop("Missing values are not allowed in primaryPeriod.")
    if(!identical(typeof(primaryPeriod), "integer")) {
        mode(primaryPeriod) <- "integer"
        warning("primaryPeriod values have been converted to integers")
        }
    ya <- array(y, c(M, J, T))
    yt.na <- apply(!is.na(ya), c(1,3), any)
    yt.na <- which(!yt.na)
    d.na <- which(is.na(primaryPeriod))
    if(!all(d.na %in% yt.na))
        stop("primaryPeriod values must be supplied for all non-missing values of y")
    increasing <- function(x) {
        x <- x[!is.na(x)]
        all(order(x) == 1:length(x))
        }
    if(!all(apply(primaryPeriod, 1, increasing)))
        stop("primaryPeriod values must increase over time for each site")
    if(class(obsCovs) == "list") {
        obsVars <- names(obsCovs)
        for(i in seq(length(obsVars))) {
            if(!(class(obsCovs[[i]]) %in% c("matrix", "data.frame")))
                stop("At least one element of obsCovs is not a matrix or data frame.")
            if(ncol(obsCovs[[i]]) != J*T | nrow(obsCovs[[i]]) != M)
                stop("At least one matrix in obsCovs has incorrect number of dimensions.")
            }
        obsCovs <- data.frame(lapply(obsCovs, function(x) as.vector(t(x))))
        }
    umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J*T))
    umf <- as(umf, "unmarkedMultFrame")
    umf@numPrimary <- numPrimary
    if(class(yearlySiteCovs) == "list") {
        yearlySiteVars <- names(yearlySiteCovs)
        for(i in seq(length(yearlySiteVars))) {
            if(!(class(yearlySiteCovs[[i]]) %in% c("matrix","data.frame")))
                stop("At least one element of yearlySiteCovs is not a matrix or data frame.")
            if(ncol(yearlySiteCovs[[i]]) != T |
                nrow(yearlySiteCovs[[i]]) != nrow(y))
                    stop("At least one matrix in yearlySiteCovs has incorrect number of dimensions.")
            }
        yearlySiteCovs <- data.frame(lapply(yearlySiteCovs, function(x)
            as.vector(t(x))))
        }
    umf@yearlySiteCovs <- yearlySiteCovs
    umf <- as(umf, "unmarkedFramePCO")
    umf@primaryPeriod <- primaryPeriod
    return(umf)
}






################ SHOW METHODS ############################################


setMethod("show", "unmarkedFrame", function(object)
{
    df <- as(object, "data.frame")
    cat("Data frame representation of unmarkedFrame object.\n")
    print(df)
})


setMethod("show", "unmarkedMultFrame",
    function(object)
{
    df <- as(object, "data.frame")
    ysc <- yearlySiteCovs(object)
    if(is.null(ysc)) {
        cat("Data frame representation of unmarkedFrame object.\n")
        print(df)
        }
    else {
        T <- object@numPrimary
        yscwide <- lapply(ysc, matrix, ncol=T, byrow=TRUE)
        df <- data.frame(df, yscwide)
        cat("Data frame representation of unmarkedFrame object.\n")
        print(df)
        }
})


############################ EXTRACTORS ##################################

# Extractor for site level covariates
setGeneric("siteCovs", function(object,...) standardGeneric("siteCovs"))

setMethod("siteCovs", "unmarkedFrame", function(object) {
    return(object@siteCovs)
})

setGeneric("yearlySiteCovs", function(object,...)
    standardGeneric("yearlySiteCovs"))
setMethod("yearlySiteCovs", "unmarkedMultFrame", function(object) {
    return(object@yearlySiteCovs)
})


setGeneric("obsCovs", function(object,...) standardGeneric("obsCovs"))
setMethod("obsCovs", "unmarkedFrame", function(object, matrices = FALSE) {
    M <- numSites(object)
    R <- obsNum(object)
    if(matrices) {
        value <- list()
        for(i in seq(length=length(object@obsCovs))){
            value[[i]] <- matrix(object@obsCovs[,i], M, R, byrow = TRUE)
        }
        names(value) <- names(object@obsCovs)
    } else {
        value <- object@obsCovs
    }
    return(value)
})


setGeneric("obsNum", function(object) standardGeneric("obsNum"))
setMethod("obsNum", "unmarkedFrame", function(object) nrow(object@obsToY))


setGeneric("numSites", function(object) standardGeneric("numSites"))
setMethod("numSites", "unmarkedFrame", function(object) nrow(object@y))


setGeneric("numY", function(object) standardGeneric("numY"))
setMethod("numY", "unmarkedFrame", function(object) ncol(object@y))


setGeneric("obsToY", function(object) standardGeneric("obsToY"))
setMethod("obsToY", "unmarkedFrame", function(object) object@obsToY)


setGeneric("obsCovs<-", function(object, value)
    standardGeneric("obsCovs<-"))
setReplaceMethod("obsCovs", "unmarkedFrame", function(object, value) {
    if(identical(class(object)[1], "unmarkedFrameDS"))
        stop("unmarkedFrameDS objects cannot have obsCovs")
    object@obsCovs <- as.data.frame(value)
    object
})


setGeneric("siteCovs<-", function(object, value)
    standardGeneric("siteCovs<-"))
setReplaceMethod("siteCovs", "unmarkedFrame", function(object, value) {
    object@siteCovs <- as.data.frame(value)
    object
})


setGeneric("yearlySiteCovs<-",
	function(object, value) standardGeneric("yearlySiteCovs<-"))
setReplaceMethod("yearlySiteCovs", "unmarkedMultFrame",
    function(object, value) {
        object@yearlySiteCovs <- as.data.frame(value)
        object
    })

setGeneric("obsToY<-", function(object, value) standardGeneric("obsToY<-"))
setReplaceMethod("obsToY", "unmarkedFrame", function(object, value) {
    object@obsToY <- value
    object
})



setGeneric("getY", function(object) standardGeneric("getY"))
setMethod("getY", "unmarkedFrame", function(object) object@y)


setGeneric("coordinates", function(object) standardGeneric("coordinates"))
setMethod("coordinates", "unmarkedFrame", function(object) {
    object@mapInfo@coordinates
})


setGeneric("projection", function(object) standardGeneric("projection"))
setMethod("projection", "unmarkedFrame", function(object) {
    object@mapInfo@projection
})

################################### SUMMARY METHODS ######################


setMethod("summary", "unmarkedFrame", function(object,...) {
    cat("unmarkedFrame Object\n\n")
    cat(nrow(object@y), "sites\n")
    cat("Maximum number of observations per site:",obsNum(object),"\n")
    mean.obs <- mean(rowSums(!is.na(getY(object))))
    cat("Mean number of observations per site:",round(mean.obs,2),"\n")
    cat("Sites with at least one detection:",
        sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
        "\n\n")
    cat("Tabulation of y observations:")
    print(table(object@y, exclude=NULL))
    if(!is.null(object@siteCovs)) {
        cat("\nSite-level covariates:\n")
        print(summary(object@siteCovs))
    }
    if(!is.null(object@obsCovs)) {
        cat("\nObservation-level covariates:\n")
        print(summary(object@obsCovs))
    }
})



setMethod("summary", "unmarkedFrameDS", function(object, ...)
{
    cat("unmarkedFrameDS Object\n\n")
    cat(object@survey, "-transect survey design", "\n", sep="")
    cat(paste("Distance class cutpoints (", object@unitsIn, "): ", sep=""),
        object@dist.breaks, "\n\n")
    cat(nrow(object@y), "sites\n")
    cat("Maximum number of distance classes per site:", ncol(getY(object)), "\n")
    mean.dc <- mean(rowSums(!is.na(getY(object))))
    cat("Mean number of distance classes per site:", round(mean.dc, 2), "\n")
    cat("Sites with at least one detection:",
        sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))), "\n\n")
    cat("Tabulation of y observations:")
    print(table(object@y, exclude=NULL))
    if(!is.null(object@siteCovs)) {
        cat("\nSite-level covariates:\n")
        print(summary(object@siteCovs))
    }
    if(!is.null(object@obsCovs)) {
        warning("Observation-level covariates cannot be used by distsamp()")
    }
})




setMethod("summary", "unmarkedMultFrame", function(object,...) {
    cat("unmarkedFrame Object\n\n")
    cat(nrow(object@y), "sites\n")
    cat("Maximum number of observations per site:",ncol(object@y),"\n")
    mean.obs <- mean(rowSums(!is.na(getY(object))))
    cat("Mean number of observations per site:",round(mean.obs,2),"\n")
    cat("Number of primary survey periods:", object@numPrimary, "\n")
    cat("Number of secondary survey periods:",
        obsNum(object) / object@numPrimary, "\n")
    cat("Sites with at least one detection:",
        sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
        "\n\n")
    cat("Tabulation of y observations:")
    print(table(object@y, exclude=NULL))
    if(!is.null(object@siteCovs)) {
        cat("\nSite-level covariates:\n")
        print(summary(object@siteCovs))
    }
    if(!is.null(object@obsCovs)) {
        cat("\nObservation-level covariates:\n")
        print(summary(object@obsCovs))
    }
    if(!is.null(object@yearlySiteCovs)) {
        cat("\nYearly-site-level covariates:\n")
        print(summary(object@yearlySiteCovs))
    }
})




################################# PLOT METHODS ###########################
# TODO:  come up with nice show/summary/plot methods for each data types.

setMethod("plot", c(x="unmarkedFrame", y="missing"),
	function (x, y, panels = 1, colorkey, strip=FALSE,
    ylab="Site", xlab="Observation", ...)
{
    y <- getY(x)
    ym <- max(y, na.rm=TRUE)
    M <- nrow(y)
    J <- ncol(y)
    y <- as.data.frame(y)
    colnames(y) <- paste("obs",1:J)
    y$site <- 1:M
    sites.per.panel <- M/panels
    y$group <- as.factor(round(seq(1,panels,length=M)))
    y2 <- melt(y, #measure.var = c("V1", "V2", "V3"),
        id.var=c("site","group"))
    if(missing(colorkey))
        colorkey <- list(at=0:(ym+1), labels=list(labels=as.character(0:ym),
            at=(0:ym)+0.5))
    levelplot(value ~ variable*site | group, y2,
        scales=list(relation="free", x=list(labels=1:J)),
        colorkey=colorkey, strip=strip, xlab=xlab, ylab=ylab, ...)
})


setMethod("hist", "unmarkedFrameDS", function(x, ...)
{
    y <- getY(x)
    dbreaks <- x@dist.breaks
    nb <- length(dbreaks)
    mids <- (dbreaks[-1] - dbreaks[-nb]) / 2 + dbreaks[-nb]
        distances <- rep(mids, times=colSums(y))
    hist(distances, breaks=dbreaks, ...)
})



################################# SELECTORS ##############################

# i is the vector of sites to extract

setMethod("[", c("unmarkedFrame", "numeric", "missing", "missing"),
    function(x, i)
{
    if(!require(reshape))
        stop("reshape package required")
    M <- numSites(x)
    if(length(i) == 0) return(x)
    if(any(i < 0) && any(i > 0))
        stop("i must be all positive or all negative indices.")
    if(all(i < 0)) { # if i is negative, then convert to positive
        i <- (1:M)[i]
        }
    y <- getY(x)[i,]
    if (length(i) == 1) {
        y <- t(y)
        }
    siteCovs <- siteCovs(x)
    obsCovs <- obsCovs(x)
    if (!is.null(siteCovs)) {
        siteCovs <- siteCovs(x)[i, , drop = FALSE]
        }
    if (!is.null(obsCovs)) {
        R <- obsNum(x)
        .site <- rep(1:M, each = R)
        obsCovs <- ldply(i, function(site) {
            subset(obsCovs, .site == site)
            })
        }
    umf <- x
    umf@y <- y
    umf@siteCovs <- siteCovs
    umf@obsCovs <- obsCovs
    umf
})


## remove obs only
### RBC: Why??? this doesn't allow umf[,c(1,1)]
setMethod("[", c("unmarkedFrame", "missing", "numeric", "missing"),
		function(x, i, j)
{
    if(!require(reshape))
        stop("reshape package required")
    y <- getY(x)
    obsCovs <- obsCovs(x)
    obsToY <- obsToY(x)
    obs.remove <- rep(TRUE, obsNum(x))
    obs.remove[j] <- FALSE
    y.remove <- t(obs.remove) %*% obsToY > 0
    y <- y[,!y.remove, drop=FALSE]
    obsCovs <- obsCovs[!rep(obs.remove, numSites(x)),, drop=FALSE]
    x@obsCovs <- obsCovs
    x@y <- y
    x@obsToY <- obsToY[!obs.remove,!y.remove, drop=FALSE]
    x
})


# i is as before and j is the obsNum to remove and corresponding y's
setMethod("[", c("unmarkedFrame","numeric", "numeric", "missing"),
		function(x, i, j)
{
    ## first remove sites
    umf <- x[i,]
    umf <- umf[,j]
    umf
})



### list is a ragged array of indices (y's) to include for each site.
### Typically useful for multilevel boostrapping.
setMethod("[", c("unmarkedFrame","list", "missing", "missing"),
    function(x, i, j)
{
    if(!require(reshape))
        stop("reshape package required")
    m <- numSites(x)
    J <- R <- obsNum(x)
    o2y <- obsToY(x)
    if (!identical(o2y, diag(R)))
        stop("Ragged subsetting of unmarkedFrames is only valid for diagonal obsToY.")
    J <- ncol(o2y)
    if (m != length(i)) stop("list length must be same as number of sites.")
    siteCovs <- siteCovs(x)
    y <- cbind(.site=1:m, getY(x))
    obsCovs <- cbind(.site=rep(1:m, each=R), obsCovs(x))

    obsCovs <- ddply(obsCovs, ~.site, function(df) {
        site <- df$.site[1]
        obs <- i[[site]]
        if (length(obs) > R)
            stop("All elements of list must be less than or equal to R.")
        obs <- c(obs, rep(NA, R-length(obs)))
        df[obs,]
        })
    obsCovs$.site <- NULL

    y <- apply(y, 1, function(row) {
        site <- row[1]
        row <- row[-1]
        obs <- i[[site]]
        obs <- c(obs, rep(NA, R-length(obs)))
        row[obs]
        })

    obsCovs(x) <- obsCovs
    x@y <- t(y)
    x
})




## for multframes, must remove years at a time
setMethod("[", c("unmarkedMultFrame", "missing", "numeric", "missing"),
		function(x, i, j)
{
    J <- obsNum(x)/x@numPrimary
    obs <- rep(1:x@numPrimary, each = J)
    years <- 1:x@numPrimary
    numPrimary <- length(j)
    obsj <- match(obs, j)
    j2 <- which(!is.na(obsj))
    u <- callNextMethod(x, i, j2)
    ysc <- yearlySiteCovs(x)
    if(!is.null(ysc)) {
        ysc <- ysc[rep(!is.na(match(years, j)), nrow(getY(x))),, drop=FALSE]
        u@yearlySiteCovs <- ysc
        }
    u@numPrimary <- numPrimary
    return(u)
})



## for multframes, must remove years at a time
setMethod("[", c("unmarkedMultFrame", "numeric", "missing", "missing"),
		function(x, i, j)
{
    if(!require(reshape))
        stop("reshape package required")
    M <- numSites(x)
    if(length(i) == 0) return(x)
    if(any(i < 0) && any(i > 0))
        stop("i must be all positive or all negative indices.")
    if(all(i < 0)) { # if i is negative, then convert to positive
        i <- (1:M)[i]
        }
    oldy <- getY(x)
    y <- oldy[i,]
    siteCovs <- siteCovs(x)
    obsCovs <- obsCovs(x)
    if (!is.null(siteCovs)) {
        siteCovs <- siteCovs(x)[i, , drop = FALSE]
        }
    if (!is.null(obsCovs)) {
        R <- obsNum(x)
        .site <- rep(1:M, each = obsNum(x)) #NULL     ## testing
        obsCovs <- ldply(i, function(site) {
            subset(obsCovs, .site == site)
            })
        }
    u <- unmarkedMultFrame(y=matrix(y, ncol=ncol(oldy)),
                           siteCovs=siteCovs,
                           obsCovs=obsCovs,
                           numPrimary=x@numPrimary)
    ysc <- x@yearlySiteCovs
    if(!is.null(ysc)) {
        T <- x@numPrimary
        sites <- rep(1:M, each=T)
        keep <- as.vector(sapply(i, function(x) which(sites %in% x)))
        ysc <- ysc[keep,, drop=FALSE]
        u@yearlySiteCovs <- ysc
        }
    u

})




setMethod("[", c("unmarkedFrameGMM", "numeric", "missing", "missing"),
		function(x, i, j)
{
    multf <- callNextMethod(x, i, j) # unmarkedMultFrame
    unmarkedFrameGMM(y=getY(multf), siteCovs=siteCovs(multf),
                     yearlySiteCovs=yearlySiteCovs(multf),
                     obsCovs=obsCovs(multf),
                     piFun=x@piFun, type=x@samplingMethod,
                     obsToY=multf@obsToY, numPrimary=multf@numPrimary)
})


setMethod("[", c("unmarkedFrameGPC", "numeric", "missing", "missing"),
		function(x, i, j)
{
    multf <- callNextMethod(x, i, j) # unmarkedMultFrame
    class(multf) <- "unmarkedFrameGPC"
    multf
})


setMethod("[", c("unmarkedFrameGPC", "missing", "numeric", "missing"),
		function(x, i, j)
{
    multf <- as(x, "unmarkedMultFrame")
    out <- callNextMethod(multf, i, j) # unmarkedMultFrame
    as(out, "unmarkedFrameGPC")
})





setMethod("[", c("unmarkedFrameGDS", "numeric", "missing", "missing"),
		function(x, i, j)
{
    multf <- callNextMethod(x, i, j) # unmarkedMultFrame
    sur <- x@survey
    if(sur=="line")
        unmarkedFrameGDS(y=getY(multf), siteCovs=siteCovs(multf),
                         yearlySiteCovs=yearlySiteCovs(multf),
                         numPrimary=x@numPrimary,
                         dist.breaks=x@dist.breaks,
                         tlength=x@tlength[i],
                         survey=sur,
                         unitsIn=x@unitsIn)
    else if(sur=="point")
        unmarkedFrameGDS(y=getY(multf), siteCovs=siteCovs(multf),
                         yearlySiteCovs=yearlySiteCovs(multf),
                         numPrimary=x@numPrimary,
                         dist.breaks=x@dist.breaks,
                         survey=sur,
                         unitsIn=x@unitsIn)
})



setMethod("[", c("unmarkedFramePCO", "numeric", "missing", "missing"),
		function(x, i, j)
{
    multf <- callNextMethod(x, i, j) # unmarkedMultFrame
    unmarkedFramePCO(y=getY(multf), siteCovs=siteCovs(multf),
                     yearlySiteCovs=yearlySiteCovs(multf),
                     obsCovs=obsCovs(multf),
                     numPrimary=x@numPrimary,
                     primaryPeriod=x@primaryPeriod[i,,drop=FALSE])
})


setMethod("[", c("unmarkedFramePCO", "missing", "numeric", "missing"),
		function(x, i, j)
{
    multf <- callNextMethod(x, i, j) # unmarkedMultFrame
    unmarkedFramePCO(y=getY(multf), siteCovs=siteCovs(multf),
                     yearlySiteCovs=yearlySiteCovs(multf),
                     obsCovs=obsCovs(multf),
                     numPrimary=length(j),
                     primaryPeriod=x@primaryPeriod[,j,drop=FALSE])
})






setMethod("head", "unmarkedFrame", function(x, n) {
    if(missing(n)) n <- 10
    umf <- x[1:n,]
    umf
})

############################### COERCION #################################

setAs("data.frame", "unmarkedFrame", function(from)
{
    umf <- formatWide(from)
    umf
})



setAs("unmarkedFrame", "data.frame", function(from)
{
    obsCovs <- obsCovs(from)
    siteCovs <- siteCovs(from)
    y <- getY(from)
    colnames(y) <- paste("y",1:ncol(y),sep=".")
    if(is.null(obsToY(from))) {
        obsNum <- ncol(y)
    } else {
        obsNum <- obsNum(from)
        }
    if(is.null(siteCovs)) siteCovs <- matrix(0,nrow(y),0)
    if(is.null(obsCovs)) {
        obsCovs <- matrix(0,nrow(y),0)
    } else {
        obsCovs <- data.frame(lapply(obsCovs,
            function(x) matrix(x, nrow(y), obsNum,byrow=T)))
        }
    df <- data.frame(y, siteCovs, obsCovs)
    df
})






