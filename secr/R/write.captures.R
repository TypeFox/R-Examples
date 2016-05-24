############################################################################################
## package 'secr'
## write.captures.R (was write.capthist - changed 2010-05-02)
## last changed 2009 03 30, 2009 06 11, 2009 07 08 2009 11 17
## revised 2010 04 01 with new argument append and character value allowed for header
## bug fixed 2010-06-07: use row names when single/multi
## 2011-09-09 trivial: ms()
## Write capture histories to text file in DENSITY format
############################################################################################

write.captures <- function (object, file='', deblank = TRUE, header = TRUE,
    append = FALSE, sess = '1', ndec = 2, covariates = FALSE, tonumeric = TRUE, ...)

{
    if (!is(object, 'capthist'))
        stop ("requires a 'capthist' object")

    if (ms(object)) {
        write.captures (object[[1]], file = file, deblank = deblank,
            header = deparse(substitute(object), control=NULL), append = append,
            sess = session(object)[1], ndec = ndec, covariates = covariates,
                        tonumeric = tonumeric, ...)
        for (i in 2:length(object)) {
            write.captures (object[[i]], file = file, deblank = deblank,
                header = FALSE, append = TRUE, sess = session(object)[i], ndec = ndec,
                covariates = covariates, tonumeric = tonumeric, ...)
        }
    }
    else {
        n <- nrow(object)
        S <- ncol(object)
        objectname <- ifelse (is.character(header),
            header, deparse(substitute(object), control=NULL))
        header <- ifelse (is.character(header), TRUE, header)
        det <- detector(traps(object))

        ID <- animalID(object)
        occ <- occasion(object)
        session <- rep(sess,length(ID))

        if (det %in% c('polygon','transect','polygonX','transectX')) {
            XY <- xy(object)
            temp <- data.frame (Session=session, ID=ID, Occasion=occ,
                x=round(XY$x,ndec), y=round(XY$y,ndec))
        }
        else if (det %in% c('signal')) {
            signal <- signal(object)
            trap <- trap(object)
            temp <- data.frame (Session=session, ID=ID, Occasion=occ, Detector=trap, Signal=signal)
        }
        else {
            trap <- trap(object)
            temp <- data.frame (Session=session, ID=ID, Occasion=occ, Detector=trap)
        }

        if (!is.null(covariates) & !is.null(covariates(object))) {
            covs <- covariates(object)
            if (is.character(covariates)) {
                covlist <- match(covariates, names(covs))
                covlist <- covlist[!is.na(covlist)]
            }
            else
                covlist <- names(covs)
            if (length(covlist)>0) {
##                for (i in 1:length(covlist))
##                covs[,i] <- as.numeric(covs[,i])
##                temp <- cbind(temp, covs[ID, covlist, drop=FALSE])
## 2014-04-05
                if (tonumeric) {
                    for (i in 1:length(covlist))
                        covs[,i] <- as.numeric(covs[,i])
                }
                temp <- cbind(temp, covs[match(ID, rownames(object)), covlist, drop=FALSE])
            }
        }

        if (header) {
            cat ("# Capture histories exported from '", objectname, "' \n", sep="", file=file)
            cat ('#', format(Sys.time(), "%a %b %d %X %Y"), '\n', append = TRUE, file=file)
            cat ('#', names(temp), '\n', append = TRUE, file=file)
            append <- TRUE
        }
        if (deblank) temp$Session <- gsub(' ','', temp$Session)
        if (deblank) temp$Session <- gsub(',','', temp$Session)
        if (any(nchar(temp$Session)>17)) {
            warning ("truncating long session names")
            temp$Session <- substring(temp$Session,1,17)
        }
        write.table(temp, file = file, row.names = FALSE, col.names = FALSE,
            append = append, quote = FALSE, ...)

    }
}
############################################################################################

