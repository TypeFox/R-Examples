############################################################################################
## package 'secr'
## DENSITY.R
## 2010 04 29, 2010 05 01
## 2010-11-21 tweaked msg
## 2010-11-21 added cutval to read.capthist
## 2010-11-26 trapcovnames passed to read.traps covnames
## 2010-12-01 suppressed warning re- detectors not recognised by Density
## 2012-10-19 telemetry detector type
## 2013-01-11 adjust call to count.fields etc to allow binary.usage
##            to be passed through to read.traps
## Write capture histories and traps to text files in DENSITY format
############################################################################################
 
write.capthist <- function (object, filestem = deparse(substitute(object)),
   sess = '1', ndec = 2, covariates = FALSE, tonumeric = TRUE, ...)

{
    sep <- list(...)$sep
    suffix <- '.txt'
    if(!is.null(sep))
        if (sep == ',')
            suffix <- '.csv'

    if (!is(object, 'capthist'))
        stop ("requires a 'capthist' object")
    ## following assignment unused 2011-09-26
    ##    det <- detector(traps(object))

    objectname <- deparse(substitute(object), control=NULL)
    if (filestem=='') {
        captfile <- ''
        trapfile <- ''
    }
    else {
        captfile <- paste(filestem, 'capt', suffix, sep='')
        trapfile <- paste(filestem, 'trap', suffix, sep='')
    }

    if (inherits(object, 'list')) {

        same <- all(sapply(traps(object)[-1], identical, traps(object)[[1]]))

        if (filestem != '') {
            if (same & !(filestem==''))
                trapfile <- paste(filestem, 'trap.txt', sep='')
            else
                trapfile <- paste(filestem, 'trap', session(object)[1], '.txt', sep='')
        }
        tempname <- ifelse(same, paste('traps(',objectname,')',sep=''),
            paste('traps(',objectname,') Session ',session(object)[1], sep=''))
        write.traps (traps(object)[[1]], file = trapfile, deblank = TRUE,
                     header = tempname, ndec = ndec, covariates = covariates, ...)
        write.captures (object[[1]], file = captfile, deblank = TRUE,
            header = deparse(substitute(object), control=NULL), append = FALSE,
            sess = session(object)[1], ndec = ndec, covariates = covariates,
                        tonumeric = tonumeric,
                        ...)
        for (i in 2:length(object)) {
            if (!same) {
                trapfile <- paste(filestem, 'trap', session(object)[i], suffix, sep='')
                write.traps (traps(object)[[i]], file = trapfile, deblank = TRUE,
                             header = FALSE, ndec = ndec, covariates = covariates,...)
            }

            write.captures (object[[i]], file = captfile, deblank = TRUE,
                header = FALSE, append = TRUE, sess = session(object)[i], ndec = ndec,
                            covariates = covariates, tonumeric = tonumeric, ...)
        }
    }
    else {
        tempname <- paste('traps(',objectname,')',sep='')
        write.captures (object, file = captfile, ..., deblank = TRUE,
            header = deparse(substitute(object), control=NULL), append = FALSE,
            sess = session(object), ndec = ndec, covariates = covariates, tonumeric = tonumeric)
        write.traps (traps(object), file = trapfile, header = tempname, covariates = covariates,  ...)
    }
}
############################################################################################

## 2015-10-02
## form joint trap layout for marking and sighting data
MRtraps <- function (marktraps, sighttraps, qvec) {
    qvec <- qvec>0
    alltraps <- rbind(marktraps, sighttraps,renumber = FALSE)
    rown <- c( rownames(marktraps), rownames(sighttraps))
    if (any(duplicated(rown)))
        rown <- c(paste('M', rownames(marktraps), sep='.'), 
                                paste('R', rownames(sighttraps), sep='.'))
    rownames(alltraps) <- rown
    tmat <- matrix(0, nrow = nrow(alltraps), ncol = length(qvec))
    tmat[1:nrow(marktraps), qvec] <- 1
    tmat[(nrow(marktraps)+1) : nrow(alltraps), !qvec] <- 1
    usage(alltraps) <- tmat
    alltraps
}
    
read.capthist <- function (captfile, trapfile, detector = 'multi', fmt = c('trapID','XY'),
                          noccasions = NULL, covnames = NULL, trapcovnames = NULL,
                          cutval = NULL, verify = TRUE, noncapt = 'NONE', tol = 0.01, 
                          snapXY = FALSE, markocc = NULL, ...) {

    fmt <- match.arg(fmt)
    dots <- match.call(expand.dots = FALSE)$...

    if ((detector %in% .localstuff$polydetectors) & !(fmt == 'XY'))
        stop ("polygon-like detectors require fmt = XY")
    if (length(captfile) != 1)
        stop ("requires single 'captfile'")
    filetype <- function(x) {
        nx <- nchar(x)
        tolower(substring(x, nx-3, nx))
    }
    if (missing(trapfile) & (detector=='telemetry'))
        trapfile <- 0  ## dummy value; not used

    countargs <- dots[names(dots) %in% names(formals(count.fields))]
    if (filetype(captfile) == '.csv')
        countargs$sep <- ','
    countargs$file <- captfile
    nfield <- max(do.call(count.fields, countargs))

    if (fmt == 'trapID') {
        nvar <- 4
        colcl <- c('character','character','character','character', rep(NA,nfield-nvar))
    }
    else {
        nvar <- 5
        colcl <- c('character','character','character',NA,NA, rep(NA,nfield-nvar))
    }

    defaultargs <- list(sep = '', comment.char = '#')
    if (filetype(captfile)=='.csv') defaultargs$sep <- ','
    captargs <- replacedefaults (defaultargs, list(...))
    captargs <- captargs[names(captargs) %in% names(formals(read.table))]
    capt <- do.call ('read.table', c(list(file = captfile, as.is = TRUE,
        colClasses = colcl), captargs) )

    ## let's be clear about this...
    if (fmt =='trapID')
        names(capt)[1:4] <- c('Session','AnimalID','Occ','Trap')
    else
        names(capt)[1:5] <- c('Session','AnimalID','Occ','X','Y')
    if (any(is.na(capt[,1:nvar])))
        stop ("missing values not allowed")

    ## allow injections 2014-07-26
    injected <- substring(capt$Occ,1,1) == '+'
    inject.time <- ifelse(injected, as.numeric(capt$Occ), 0)
## need to retrospectively zero this detection...
## and perhaps set prior to NA?
    capt$Occ <- as.numeric(capt$Occ)

    ## assumes file= is first argument of read.traps
    ## allows for multiple trap files
    defaultdots <- list(sep = '', comment.char = '#')
    if (filetype(trapfile[1])=='.csv') defaultdots$sep <- ','
    dots <- replacedefaults (defaultdots, list(...))
    readtraps <- function (x, mo) {
        if (detector == 'telemetry') {
            buffer <- c(-10,10)
            trps <- expand.grid(x = range(capt$X)+buffer,
                                y = range(capt$Y)+buffer)[c(1,3,4,2,1),]
            rownames(trps) <- 1:5
            class(trps) <- c('traps', 'data.frame')
            attr(trps, 'polyID') <- factor(rep(1,nrow(trps)))
            attr(trps, 'detector') <- 'telemetry'
            trps
        }
        else {
            
            do.call ('read.traps', c(list(file = x), detector = detector,
                                 list(covnames = trapcovnames), 
                                 list(markocc = mo), dots) )
        }
    }
   
    ## use mapply to carry session-varying markocc?
    molist <- inherits(markocc, 'list')
    if (molist) {
        trps <- mapply(readtraps, trapfile, markocc, SIMPLIFY = FALSE)
    }
    else
        trps <- sapply(trapfile, readtraps, mo = markocc, simplify = FALSE)
    
    if (length(trps)==1) trps <- trps[[1]]

    temp <- make.capthist(capt, trps, fmt = fmt,  noccasions = noccasions,
        covnames = covnames, sortrows = TRUE, cutval = cutval,
        noncapt = noncapt, tol = tol, snapXY = snapXY)

    ## temporary (?) way to attach injection times, not changng make.capthist
    ## vector of occasion after inj at which first available for detection
    ## may be all zero
    ## consider NA instead of zero for non-injected animals
    ## include in verify.capthist checks that
    ## -- animal not detected before injected
    ## -- length inject.time equals number of animals
    ## -- animal not injected twice
    ## 2014-07-27
    attr(temp, 'inject.time') <- inject.time

    if (verify) verify(temp)
    temp
}

############################################################################################

