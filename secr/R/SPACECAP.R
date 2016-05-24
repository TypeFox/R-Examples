############################################################################################
## package 'secr'
## SPACECAP.R
## Write captures and detector locations to text files SPACECAP format
## Read captures and detector locations from text files in SPACECAP format
## 2010 04 10, 2010 04 30, 2012 12 18
############################################################################################

write.SPACECAP <- function (object, mask = NULL, buffer = 100, ndec = 2, filestem='') {

    objectname <-  deparse(substitute(object), control=NULL)
    if (!is(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (detector(traps(object)) != 'proximity')
        stop ("requires proximity detectors")
    if (inherits(object, 'list'))
        stop ("multiple sessions not compatible with SPACECAP")

    AC <- ifelse (filestem == '', '', paste(filestem,'AC.csv',sep='') )
    TD <- ifelse (filestem == '', '', paste(filestem,'TD.csv',sep='') )
    SS <- ifelse (filestem == '', '', paste(filestem,'SS.csv',sep='') )

    nanimal <- nrow(object)
    nocc <- ncol(object)
    trps <- traps(object)
    ntrap <- nrow(trps)
    trps$x <- round(trps$x,ndec)
    trps$y <- round(trps$y,ndec)

    ## traps
    temp <- trps
    if (!is.null(usage(trps))) {
        usage(traps) <- (usage(traps)>0) * 1  ## force to binary for SPACECAP 2012-12-17
        temp <- cbind(temp, usage(trps))
    }
    else
        temp <- cbind(temp, matrix(1, nrow = ntrap, ncol = nocc))
    cat (c('LOC_ID', 'X_Coord','Y_Coord',1:nocc), sep=',', file=TD)
    cat ('\n', append = TRUE, file=TD)
    write.table(temp, file = TD, append = TRUE, quote = FALSE, sep=',',
        row.names = 1:ntrap, col.names = FALSE)

    ## animals
    trap <- rep(rep( 1:ntrap, rep(nanimal * nocc, ntrap)), abs(object))
    ID <- rep(rep( 1:nanimal, nocc * ntrap ), abs(object))
    occ <- rep(rep (rep (1:nocc, ntrap), rep(nanimal, nocc*ntrap)), abs(object))
    temp <- data.frame (trap=trap, ID=ID, occ=occ)
    cat (c('LOC_ID','ANIMAL_ID','SO'), sep=',', file=AC)
    cat ('\n', append = TRUE, file=AC)
    write.table(temp, file = AC, row.names = FALSE, col.names = FALSE,
        append = TRUE, sep = ',', quote = FALSE)

    ## mask
    if  (is.null(mask))
        mask <- make.mask(trps, buffer=buffer)
    else {
        if (!inherits(mask, 'mask'))
            stop ("requires valid 'mask'")
        if (inherits(mask, 'list'))
            stop ("multi-session 'mask' is not suitable")
    }
    mask$habitat <- rep(1,nrow(mask))
    cat (c('X_COORD','Y_COORD','HABITAT'), sep=',', file=SS)
    cat ('\n', append = TRUE, file=SS)
    write.table(mask, file = SS, append = TRUE, quote = FALSE, sep=',',
        row.names = FALSE, col.names = FALSE)

}

############################################################################################

read.SPACECAP <- function (AC, TD, detector = 'proximity', session = '1') {

    capt <- read.csv(AC)
    trps <- read.csv(TD)

    if (any (!(c('LOC_ID','ANIMAL_ID','SO') %in% names(capt))))
        stop ("AC not valid")
    if (any (!(c('LOC_ID', 'X_Coord','Y_Coord') %in% names(trps))))
        stop ("TD not valid")

    used <- trps[4:ncol(trps)]
    nocc <- ncol(used)
    traps <- trps[, c('X_Coord','Y_Coord')]
    names(traps) <- c('x','y')
    class(traps) <- c("traps", "data.frame")
    detector(traps) <- detector
    covariates(traps) <- NULL
    if (any(used != 1))  ## don't bother if all detectors always used
        usage(traps) <- as.matrix(used)   ## 2012-12-18
    ux <- unique(traps$x)
    uy <- unique(traps$y)
    attr(traps,'spacex') <- ifelse (length(ux)>1, min(dist(ux)), NA)
    attr(traps,'spacey') <- ifelse (length(uy)>1, min(dist(uy)), NA)
    spacing(traps) <- spacing(traps)   ## !!
    capt$session <- rep(session,nrow(capt))
    capt <- capt[,c('session','ANIMAL_ID','SO','LOC_ID')]
    make.capthist(capt, traps, fmt = 'trapID', noccasions = nocc)
}

############################################################################################

