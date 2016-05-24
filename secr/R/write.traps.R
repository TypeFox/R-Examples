############################################################################################
## package 'secr'
## write.traps.R
## last changed 2009 06 11 2009 11 17 2010 04 30
## Write detector locations to text file in DENSITY format
## should remove conflict between row and ...
############################################################################################

write.traps <- function (object, file='', deblank = TRUE, header = TRUE, ndec = 2,
    covariates = FALSE, ...) {

    objectname <- ifelse (is.character(header),
        header, deparse(substitute(object), control=NULL))
    header <- ifelse (is.character(header), TRUE, header)

    if (!is(object, 'traps'))
        stop ("requires a 'traps' object")
    n <- nrow(object)
    object$x <- round(object$x,ndec)
    object$y <- round(object$y,ndec)

    # purge blanks from names
    if (deblank) row.names(object) <- gsub(' ','',row.names(object))

    poly <- detector(object) %in% c('polygon', 'polygonX')
    transect <- detector(object) %in% c('transect', 'transectX')
    if (poly) {
        temp <- cbind (polyID=polyID(object), x=object$x, y = object$y)
    }
    else if (transect) {
        temp <- cbind (transectID=transectID(object), x=object$x, y = object$y)
    }
    else {
        temp <- object
        if (!is.null(usage(object))) temp <- cbind(temp,usage(object))
    }

    covlist <- numeric(0)
    if (!is.null(covariates) & !is.null(covariates(object))) {
        covs <- covariates(object)
        if (is.character(covariates)) {
            covlist <- match(covariates, names(covs))
            covlist <- covlist[!is.na(covlist)]
        }
        else
            covlist <- names(covs)

        if (length(covlist)>0) {
            covnames <- paste(covlist, collapse=' ')
            covs <- covs[, covlist, drop=FALSE]
            ## assume order of levels of polyID matches order in object
            if (poly | transect)
                covs <- covs[as.numeric(polyID(object)), , drop=FALSE]
            for (i in 1:length(covlist))
                covs[,i] <- as.numeric(covs[,i])
            covs <- apply(covs,1,paste, collapse=' ')
            covs <- paste ('/',covs)
            temp <- cbind(temp, covs)
        }
    }

    if (header) {
        cat ("# Detector locations exported from '", objectname, "' \n",
             sep = "", file = file)
        cat ('#', format(Sys.time(), "%a %b %d %X %Y"), '\n', append = TRUE,
             file = file)
        if (poly)
            headtext <- '# polyID  x  y'
        else
        if (transect)
            headtext <- '# transectID  x  y'
        else
            headtext <- '# Detector  x  y'
        if (length(covlist)>0)
            headtext <- paste(headtext, covnames, sep = ' / ')

        cat (headtext, '\n', append = TRUE, file=file)
    }

    write.table(temp, file = file, append = header,
        row.names = !poly & !transect, col.names = FALSE, quote = FALSE, ...)

}
###############################################################################
