#######################################################################################
## utility.R
#######################################################################################

## 2012-04-13 pointsInPolygon extended to allow mask 'polygon'
## 2012-07-25 nclusters() added
## 2012-09-04 leadingzero moved here
## 2012-10-21 HN, HZ etc moved here
## 2012-10-28 model.string moved here
## 2012-11-03 gradient moved here
## 2012-11-03 several other functions moved here from functions.r
## 2012-12-03 make.lookup moved here
## 2013-04-12 getknownclass()
## 2013-04-12 getnmix()
## 2013-04-13 h.levels()
## 2013-04-20 new detection functions
## 2013-06-06 fullbeta function (was part of secrloglik)
## 2013-06-08 get.nmix extended to allow hcov and not h2/h3 model for g0,sigma
## 2013-06-15 inflate()
## 2013-07-19 a0 in valid.detectpar
## 2013-10-28 fixed mlogit.untransform
## 2013-10-29 fixpmix moved here: code was in model.average.R and methods.R

## 2013-06-17 I have so far resisted the temptation to add HCU hazard cumulative uniform df
##            based on Horne & Garton 2006
## 2013-11-09 pointsInPolygon bug fix
## 2013-11-09 getbinomN moved from pdot
## 2013-11-16 xy2CH telemetry; patched for covariates 2013-12-02
## 2013-11-20 getcoord from SpatialPolygons
## 2014-03-18 complete.beta functions for fixedbeta
## 2014-08-19 moved secr.lpredictor from secrloglik.R
## 2014-08-21 smooths function
## 2014-08-29 masklength()
## 2014-09-01 valid.userdist() revised 2014-10-13
## 2014-09-09 make.lookup uses values rounded to 10 dp
## 2014-09-10 getnoneucnames()
## 2014-09-10 getcellsize()
## 2014-12-04 group.levels error check now precedes attempt to extract by groups
## 2015-01-29 improved error message in secr.lpredictor
## 2015-03-31 nparameters
## 2015-04-03 mapbeta moved from score.test.R
## 2015-11-02 xyinpoly moved from verify.R
## 2015-11-17 improved robustness of mapbeta when real parameter missing from new model 
## 2016-01-08 addzerodf function used by join()
## 2016-01-08 addzeroCH function used by sim.resight()
#######################################################################################

# Global variables in namespace
#
## define a local environment (in namespace?) for temporary variables e.g. iter
## e.g. Roger Peng https://stat.ethz.ch/pipermail/r-devel/2009-March/052883.html

.localstuff <- new.env()
.localstuff$validdetectors <- c('single','multi','proximity','count',
    'polygonX', 'transectX', 'signal', 'signalnoise', 'polygon', 'transect', 'times',
                                'unmarked','presence','telemetry')
.localstuff$simpledetectors <- c('single','multi','proximity','count')
.localstuff$individualdetectors <- c('single','multi','proximity','count',
    'polygonX', 'transectX', 'signal', 'signalnoise', 'polygon', 'transect', 'times',
                                     'telemetry')
.localstuff$pointdetectors <- c('single','multi','proximity','count',
    'signal', 'signalnoise', 'unmarked','presence')
.localstuff$polydetectors <- c('polygon','transect','polygonX','transectX','telemetry')
.localstuff$exclusivedetectors <- c('single','multi','polygonX','transectX')
## 'signal' is not a count detector 2011-02-01
.localstuff$countdetectors <- c('count','polygon','transect','unmarked','telemetry')
.localstuff$detectors3D <- c('proximity','count','signal','signalnoise','polygon',
                             'transect','times','unmarked','presence','telemetry')
.localstuff$iter <- 0
.localstuff$iter2 <- 0
.localstuff$detectionfunctions <-
        c('halfnormal',
      'hazard rate',
      'exponential',
      'compound halfnormal',
      'uniform',
      'w exponential',
      'annular normal',
      'cumulative lognormal',
      'cumulative gamma',
      'binary signal strength',
      'signal strength',
      'signal strength spherical',
      'signal-noise',
      'signal-noise spherical',
      'hazard halfnormal',
      'hazard hazard rate',
      'hazard exponential',
      'hazard annular normal',
      'hazard cumulative gamma')

.localstuff$DFN <- c('HN', 'HR', 'EX', 'CHN', 'UN', 'WEX', 'ANN', 'CLN', 'CG',
                     'BSS', 'SS', 'SSS', 'SN', 'SNS', 'HHN', 'HHR', 'HEX', 'HAN', 'HCG')

detectionfunctionname <- function (fn) {
    .localstuff$detectionfunctions[fn+1]
}

detectionfunctionnumber <- function (detname) {
    dfn <- match (toupper(detname), .localstuff$DFN)
    if (is.na(dfn))
        dfn <- match (tolower(detname), .localstuff$detectionfunctions)
    if (is.na(dfn))
        stop ("unrecognised detection function ", detname)
    dfn-1
}
parnames <- function (detectfn) {
    parnames <- switch (detectfn+1,
        c('g0','sigma'),   ## 0
        c('g0','sigma','z'),
        c('g0','sigma'),
        c('g0','sigma','z'),
        c('g0','sigma'),
        c('g0','sigma','w'),
        c('g0','sigma','w'),
        c('g0','sigma','z'),
        c('g0','sigma','z'),
        c('b0','b1'),
        c('beta0','beta1', 'sdS'),    ## include cutval?
        c('beta0','beta1', 'sdS'),    ## include cutval?
        c('beta0','beta1', 'sdS','muN','sdN'),
        c('beta0','beta1', 'sdS','muN','sdN'),
        c('lambda0','sigma'),
        c('lambda0','sigma','z'),
        c('lambda0','sigma'),
        c('lambda0','sigma','w'),
        c('lambda0','sigma','z'),
        ,
        c('g0','sigma')    ## 20
    )
}
getdfn <- function (detectfn) {
    switch (detectfn+1, HN, HR, EX, CHN, UN, WEX, ANN, CLN, CG, BSS, SS, SSS,
                       SN, SNS, HHN, HHR, HEX, HAN, HCG)
}

valid.detectfn <- function (detectfn, valid = c(0:3,5:18)) {
# exclude 4 uniform: too numerically flakey
    if (is.null(detectfn))
        stop ("requires 'detectfn'")
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)
    if (!(detectfn %in% valid))
        stop ("invalid detection function")
    detectfn
}

valid.detectpar <- function (detectpar, detectfn) {
    if (is.null(detectpar) | is.null(detectfn))
        stop ("requires valid 'detectpar' and 'detectfn'")

    ## 2013-07-19, 2013-10-22
    ## replace a0 with g0 or lambda0 as appropriate to detectfn
    if ('a0' %in% names(detectpar)) {
        aname <- if (detectfn %in% 0:8) 'g0' else 'lambda0'
        lambda0 <- detectpar[['a0']] / 2 / pi / detectpar[[2]]^2 * 10000
        detectpar[[aname]] <- if (detectfn %in% 0:8) 1-exp(-lambda0) else lambda0
    }

    if (!all(parnames(detectfn) %in% names(detectpar)))
        stop ("requires 'detectpar' ", paste(parnames(detectfn), collapse=','),
            " for ", detectionfunctionname(detectfn), " detectfn")
    detectpar[parnames(detectfn)]
}

valid.model <- function(model, CL, detectfn, hcov, userdist, sessioncovnames) {
    ## 2014-08-25
    if (any(sapply(model, badsmooths)))
        warning("smooth term may be unsuitable for secr: does not specify k or fx where required")
}

getuserdistnames <- function (userdist) {
    ## return the names of any supplementary arguments of user-provided function
    ## for non-euclidean distance computations
    if (is.function(userdist)) {
        distnames <- try(userdist(), silent = TRUE)
        if (!is.character(distnames))
            stop("invalid userdist function - ",
                 "should return parameter names when called with no arguments")
        distnames
    }
    else
        character(0)
}

valid.pnames <- function (details, CL, detectfn, alltelem, sighting, nmix) {
    ## modelled parameters
    pnames <- switch (detectfn+1,
        c('g0','sigma'),           # 0 halfnormal
        c('g0','sigma','z'),       # 1 hazard rate
        c('g0','sigma'),           # 2 exponential
        c('g0','sigma','z'),       # 3
        c('g0','sigma'),           # 4
        c('g0','sigma','w'),       # 5
        c('g0','sigma','w'),       # 6
        c('g0','sigma','z'),       # 7
        c('g0','sigma','z'),       # 8
        c('b0','b1'),              # 9
        c('beta0','beta1','sdS'),  # 10
        c('beta0','beta1','sdS'),  # 11
        c('beta0','beta1','sdS'),  # 12  cf parnames() in utility.R: muN, sdN?
        c('beta0','beta1','sdS'),  # 13  cf parnames() in utility.R: muN, sdN?
        c('lambda0','sigma'),      # 14 hazard halfnormal
        c('lambda0','sigma','z'),  # 15 hazard hazard rate
        c('lambda0','sigma'),      # 16 hazard exponential
        c('lambda0','sigma','w'),  # 17
        c('lambda0','sigma','z'))  # 18

    if (details$param %in% c(2,6))
        pnames[1] <- 'esa'
    if (details$param %in% c(3,5))
        pnames[1] <- 'a0'
    if (details$param %in% 4:6) {
        pnames[2] <- 'sigmak'
        pnames <- c(pnames, 'c')
    }
    if (!CL)
      pnames <- c('D', pnames)
    if ('noneuc' %in% getuserdistnames(details$userdist))
      pnames <- c(pnames, 'noneuc')
    if (sighting)
      pnames <- c(pnames, 'pID')
    if (alltelem) {
        rnum <- match(c('D','lambda0','a0','esa','g0'), pnames)
        rnum[is.na(rnum)] <- 0
        pnames <- pnames[-rnum]
    }
    if (nmix>1)
        pnames <- c(pnames, 'pmix')
    pnames
}
#-------------------------------------------------------------------------------

valid.userdist <- function (userdist, detector, xy1, xy2, mask) {
    ## 2014-09-01, 2014-09-10, 2014-10-13
    if (is.null(userdist)) {
        ## default to Euclidean distance
        result <- edist(xy1, xy2)  ## as.matrix(dist(rbind(xy1, xy2)))[1:nr1, (nr1+1):(nr1+nr2)]
    }
    else {
        if (detector %in% .localstuff$polydetectors) {
            stop ("userdist cannot be used with polygon detector types;")
        }
        if (is.function(userdist))
        {
            OK <- getuserdistnames(userdist) %in% names(covariates(mask))
            if ((length(OK)>0) & !all(OK))
                stop ("covariates required by userdist function not in mask : ",
                      paste(getuserdistnames(userdist)[!OK], collapse=','))
            result <- do.call(userdist, c(list(xy1, xy2, mask)))
        }
        else {
            result <- userdist
        }
        if (!all(dim(result) == c(nrow(xy1), nrow(xy2))))
            stop ("invalid distance matrix dim = ", dim(result)[1], ',', dim(result)[2])
        baddist <- (!is.finite(result)) | (result<0) | is.na(result)
        if (any(baddist)) {
            warning ("replacing infinite, negative and NA userdist values with 1e10")
            result[baddist] <- 1e10
        }
    }
    result
}
#-------------------------------------------------------------------------------

new.param <- function (details, model, CL) {
    esa <- 'esa' %in% names(model)
    a0 <- 'a0' %in% names(model)
    sigmak <- 'sigmak' %in% names(model)
    newparam <- details$param
    if (esa & !sigmak) {
        newparam <- 2
    }
    if (a0 & !sigmak) {
        newparam <- 3
    }
    if (sigmak) {
        if (esa) {
            newparam <- 6
        }
        else {
            if (CL)
                stop ("sigmak parameterization requires full model, not CL, unless also 'esa'")
            newparam <- ifelse(a0, 5, 4)
        }
    }
    if (newparam  != details$param)
        warning ("Using parameterization details$param = ", newparam)
    newparam
}

detectorcode <- function (object, MLonly = TRUE) {
    ## numeric detector code from traps object
    detcode <- switch (detector(object),
        single = -1,
        multi = 0,
        proximity = 1,
        count = 2,
        polygonX = 3,
        transectX = 4,
        signal = 5,
        polygon = 6,
        transect = 7,
        times = 8,
        ## cue = 9, defunct 2015-10-01 secr 2.10.0
        unmarked = 10,
        presence = 11,
        signalnoise = 12,
        telemetry = 13,
        -2)
    if (MLonly) {
        detcode <- ifelse (detcode==-1, 0, detcode)
        if (detcode<0)
            stop ("Unrecognised detector type")
    }
    detcode
}

## 2013-06-16
ndetectpar <- function (detectfn) {
    length(parnames(detectfn))
}

replacedefaults <- function (default, user) replace(default, names(user), user)

discreteN <- function (n, N) {
    tN <- trunc(N)
    if (N != tN) tN + sample (x = c(1,0), prob = c(N-tN, 1-(N-tN)),
        replace = T, size = n)
    else rep(tN,n)
}

ndetector <- function (traps) {
    if (detector(traps) %in% .localstuff$polydetectors)
        length(levels(polyID(traps)))
    else
        nrow(traps)
}

memo <- function (text, trace) {
    ## could use message(text), but does not immediately flush console
    if (trace) { cat (text, '\n')
    flush.console() }
}

insertdim <- function (x, dimx, dims) {
  ## make vector of values
  ## using x repeated so as to fill array
  ## with dim = dims and the x values occupying dimension(s) dimx
  olddim <- 1:length(dims)
  olddim <- c(olddim[dimx], olddim[-dimx])
  temp <- array (dim=c(dims[dimx], dims[-dimx]))
  tempval <- array(dim=dims[dimx])
  if (length(x) > length(tempval))
      tempval[] <- x[1:length(tempval)]
  else
      tempval[] <- x     ## repeat as needed
  temp[] <- tempval  ## repeat as needed
  if (is.factor(x))
    factor(levels(x), levels=levels(x))[aperm(temp, order(olddim))]   ## 2010 02 25
  else
    as.vector(aperm(temp, order(olddim)))
}

pad1 <- function (x, n) {
## pad x to length n with dummy (first value)
    if (is.factor(x)) {
        xc <- as.character(x)
        xNA <- c(xc, rep(xc[1], n-length(xc)))
        out <- factor(xNA, levels=levels(x))
    }
    else out <- c(x, rep(x[1], n-length(x)))
    out
}

padarray <- function (x, dims) {
    temp <- array(dim=dims)
    dimx <- dim(x)
    if (length(dimx)<2 | length(dimx)>3)
        stop ("invalid array")
    if (length(dimx)>2) temp[1:dimx[1], 1:dimx[2], 1:dimx[3]] <- x
    else temp[1:dimx[1], 1:dimx[2]] <- x
    temp
}

## regularize a list of formulae
## added 2009 08 05
stdform <- function (flist) {
    LHS <- function (form) {
        trms <- as.character (form)
        if (length(trms)==2) '' else trms[2]
    }
    RHS <- function (form) {
        trms <- as.character (form)
        if (length(trms)==3) as.formula(paste(trms[c(1,3)])) else form
    }
    lhs <- sapply(flist, LHS)
    temp <- lapply(flist, RHS)
    if (is.null(names(flist))) names(temp) <- lhs
    else names(temp) <- ifelse(names(flist) == '', lhs, names(flist))
    temp
}

## Start of miscellaneous functions

invlogit <- function (y) 1/(1+exp(-y))   # plogis(y)
logit    <- function (x) log(x/(1-x))    # qlogis(x), except for invalid argument
sine     <- function (x) asin (x*2-1)
invsine  <- function (y) (sin(y)+1) / 2
odds     <- function (x) x / (1-x)
invodds  <- function (y) y / (1+y)

lnbinomial <- function (x,size,prob) {
  lgamma (size+1) - lgamma (size-x+1) - lgamma (x+1) +
      x * log(prob) + (size-x) * log (1-prob)
}
############################################################################################
## moved from methods.r 2012-10-28

model.string <- function (model, userDfn) {
    if (!is.null(userDfn)) {
        if (!is.null(model$D))
            model$D <- paste('~userD', userDfn('name'), sep='.')
    }
    temp <- paste (names(model), as.character(model), collapse=' ', sep='')
    temp
}
fixed.string <- function (fixed) {
    if (is.null(fixed) | length(fixed)==0) 'none'
    else paste (names(fixed), as.character(fixed), collapse=', ', sep=' = ')
}
############################################################################################

var.in.model <- function(v,m) v %in% unlist(lapply(m, all.vars))

get.nmix <- function (model, capthist, hcov) {
    model$D <- NULL  ## ignore density model
    model$pmix <- NULL ## pmix alone cannot make this a mixture model
    nmix <- 1
    if (any(var.in.model('h2', model))) {
        nmix <- 2
        if (any(var.in.model('h3', model)))
            stop ("do not combine h2 and h3")
    }
    if (any(var.in.model('h3', model))) {
        nmix <- 3
    }
    if ((nmix == 1) & (!is.null(hcov))) {
        if (ms(capthist))
            capthist <- capthist[[1]]
        if (is.factor(covariates(capthist)[,hcov]))
            lev <- levels(covariates(capthist)[,hcov])
        else
            lev <- levels(factor(covariates(capthist)[,hcov]))
        if (all(is.na(covariates(capthist)[,hcov])))
            stop ("hcov missing for all individuals, but detection model invariant")
        if (length(lev) < 2)
            stop ("hcov covariate not found or has fewer than 2 levels")
        if (length(lev) > 2)
            warning ("hcov covariate has more than 2 levels; using only first two")
        nmix <- 2
    }
    nmix
}
############################################################################################

fixpmix <- function(x, nmix) {

    ## x is a list with component pmix that is a matrix (dataframe)
    ## with columns 'estimate' and 'se' (and possibly others)
    ## referring to the linear predictor of pmix (i.e. on mlogit
    ## scale) and rows corresponding to rows in newdata
    ## (i.e. arbitrary combinations of predictors, including mixture
    ## class h2 or h3)

    ####################################################
    ## It is necessary that newdata include all levels
    ## of the mixture class.
    ####################################################

    ## 2013-10-29
    ## assuming mixture is always last dimension...

    ## previously used in collate, model.average and predict.secr
    ## 2015-09-30 incorporated in secr.lpredictor

    temp <- matrix(x$pmix[,'estimate'], ncol = nmix)
    if (nmix==2) temp[,x$pmix[,'h2']] <- x$pmix[,'estimate']
    if (nmix==3) temp[,x$pmix[,'h3']] <- x$pmix[,'estimate']
    temp2 <- apply(temp, 1, clean.mlogit)
    x$pmix[,'estimate'] <- as.numeric(t(temp2))
    if (nmix==2)
        x$pmix[as.numeric(x$pmix$h2)==1,'se'] <- x$pmix[as.numeric(x$pmix$h2)==2,'se']
    else
        x$pmix[,'se'] <- rep(NA, nrow(x$pmix))   ## don't know how
    x
}

############################################################################################

add.cl <- function (df, alpha, loginterval, lowerbound = 0) {

## add lognormal or standard Wald intervals to dataframe with columns
## 'estimate' and 'SE.estimate'
## lowerbound added 2011-07-15
    z <- abs(qnorm(1-alpha/2))
    if (loginterval) {
        delta <- df$estimate - lowerbound
        df$lcl <- delta / exp(z * sqrt(log(1 + (df$SE.estimate /
                        delta)^2))) + lowerbound
        df$ucl <- delta * exp(z * sqrt(log(1 + (df$SE.estimate /
                        delta)^2))) + lowerbound
    }
    else {
        df$lcl <- pmax(lowerbound, df$estimate - z * df$SE.estimate)
        df$ucl <- df$estimate + z * df$SE.estimate
    }
    df
}

###############################################################################

spatialscale <- function (object, detectfn, session = '') {
    if (inherits(object, 'secr')) {
        if (ms(object))
            detpar <- detectpar(object)[[session]]
        else
            detpar <- detectpar(object)
        cutval <- object$details$cutval
    }
    else {
        detpar <- object
        cutval <- object$cutval
    }
    if (!is.null(detpar$sigma)) detpar$sigma
    else if (detectfn == 10) {
        (cutval - detpar$beta0) / detpar$beta1
    }
    else if (detectfn == 11) {
        d11 <- function(d, beta0, beta1, c) beta0 +
            beta1 * (d-1) - 10 * log10(d^2) - c
        interval <- c(0,10 * (cutval - detpar$beta0) / detpar$beta1)
        uniroot (d11, interval, detpar$beta0, detpar$beta1, cutval)$root
    }
    else if (detectfn == 9) {
#        (0.5 - detpar$b0) / detpar$b1
        - 1 / detpar$b1   ## 2010-11-01
    }
    else stop ("unrecognised detectfn")
}

###############################################################################

pointsInPolygon <- function (xy, poly, logical = TRUE) {
    xy <- matrix(unlist(xy), ncol = 2)  ## in case dataframe
    if (inherits(poly, 'SpatialPolygonsDataFrame')) {
        xy <- SpatialPoints(xy)
        ## 2013-04-20 update for deprecation of 'overlay'
        ## OK <- overlay (xy, poly)
        ## 2014-12-11
        proj4string(poly) <- CRS()
        OK <- sp::over (xy, poly)
        ## bug fix 2013-11-09
        if (!is.null(dim(OK)))
            OK <- OK[,1]
        !is.na(OK)
    }
    else if (inherits(poly, 'mask')) {  # 2012-04-13
        if (ms(poly))
            stop ("multi-session masks not supported")
        sp <- spacing(poly)
        minx <- min(poly$x, na.rm = TRUE)
        miny <- min(poly$y, na.rm = TRUE)
        mask <- sweep(poly, MARGIN = 2, FUN = '+', STATS = c(-minx, -miny))
        mask <- round(mask/sp) + 1
        xy <- sweep(xy, MARGIN = 2, FUN = '+', STATS = c(-minx, -miny))
        xy <- round(xy/sp) + 1
        ## 2013-03-06 tweak
        ## xy[xy<0] <- NA
        xy[xy<=0] <- NA
        xy[,1][xy[,1]>max(mask$x, na.rm = TRUE)] <- NA
        xy[,2][xy[,2]>max(mask$y, na.rm = TRUE)] <- NA

        maskmatrix <- matrix(0, ncol = max(mask$y, na.rm = TRUE), nrow = max(mask$x, na.rm = TRUE))
        maskmatrix[as.matrix(mask)] <- 1:nrow(mask)
        inside <- maskmatrix[as.matrix(xy)]
        inside[is.na(inside)] <- 0
        if (logical)
            inside <- inside > 0
        inside
    }
    else {
        checkone <- function (xy1) {
            temp <- .C('inside',  PACKAGE = 'secr',
                as.double (xy1),
                as.integer (0),
                as.integer (np-1),
                as.integer (np),
                as.double (poly),
                result = integer(1))
            as.logical(temp$result)
        }
        poly <- matrix(unlist(poly), ncol = 2)  ## in case dataframe
        np <- nrow(poly)
        apply(xy, 1, checkone)
    }
}
###############################################################################
## logical for whether object specifies userDfn

userD <- function (object) {
    if (!inherits(object, 'secr'))
        stop ("requires secr fitted model")
    !is.null(object$details$userDfn)
}

###############################################################################

## mean and SD if x numeric
## was statfn 2011-11-08
getMeanSD <- function(xy) {
    MeanSD <- function (x) {
        if (is.numeric(x))
            c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
        else
            c(NA,NA)
    }
    as.data.frame (apply(xy, 2, MeanSD))
}
###############################################################################

nclusters <- function (capthist) {
    if (ms(capthist)) {
	lapply(capthist, nclusters)
    }
    else 	{
        nmash <- attr(capthist, 'n.mash')
        ifelse (is.null(nmash), 1, length(nmash))
    }
}
###############################################################################
## moved here from make.grid 2012-09-04

# leadingzero <- function (x) {
#     formatC(x, width=max(nchar(x)), flag='0')  ## returns character value
# }

## clunky but effective re-write 2012-09-04
leadingzero <- function (x) {
    x <- as.character(x)
    w <- max(nchar(x))
    n0 <- function(n) paste(rep('0',n), collapse='')
    paste(sapply(w-nchar(x), n0), x, sep='')
}

###############################################################################

## Detection functions

HN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]
    g0 * exp (-r^2 / 2 / sigma^2)
}
HR <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    g0 * (1 - exp (-(r / sigma)^-z))
}
EX <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]
    g0 * exp (-r / sigma)
}
UN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]
    ifelse (r<=sigma, g0, 0)
}
CHN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    g0 * ( 1 - (1 - exp (-r^2 / 2 / sigma^2)) ^ z )
}
WEX <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
    ifelse(r<=w, g0, g0*exp (-(r-w) / sigma))
}
ANN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
    g0 * exp (-(r-w)^2 / 2 / sigma^2)
}
CLN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    CV2 <- (z/sigma)^2
    sdlog <- log(1 + CV2)^0.5
    meanlog <- log(sigma) - sdlog^2/2
    g0 * plnorm(r, meanlog, sdlog, lower.tail = FALSE)
}
CG <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    g0 * pgamma(r, shape=z, scale=sigma/z, lower.tail = FALSE)
}
CN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    x <- z * (r - sigma)
    g0 * (1 + (1 - exp(x)) / (1 + exp(x)))/2
}
BSS <- function (r, pars, cutval) {
    b0 <- pars[1]; b1 <- pars[2]
    gam <- -(b0 + b1 * r);
    pnorm (gam, mean=0, sd=1, lower.tail=FALSE)
}
SS <- function (r, pars, cutval) {
    beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3]
    if (is.null(cutval))
        stop ("require 'details$cutval' for signal strength plot")
    mu <- beta0 + beta1 * r
    1 - pnorm (q=cutval, mean=mu, sd=sdS)
}
SSS <- function (r, pars, cutval) {
    beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3]
    if (is.null(cutval))
        stop ("require 'details$cutval' for signal strength plot")
    ## spherical so assume distance r measured from 1 m
    mu <- beta0 - 10 * log ( r^2 ) / 2.302585 + beta1 * (r-1)
    mu[r<1] <- beta0
    1 - pnorm (q=cutval, mean=mu, sd=sdS)
}
SN <- function (r, pars, cutval) {
    beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3];
    muN <- pars[4]; sdN <- pars[5]
    muS <- beta0 + beta1 * r
    1 - pnorm (q=cutval, mean=muS-muN, sd=sqrt(sdS^2+sdN^2))
}
SNS <- function (r, pars, cutval) {
    beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3];
    muN <- pars[4]; sdN <- pars[5]
    ## spherical so assume distance r measured from 1 m
    muS <- beta0 - 10 * log ( r^2 ) / 2.302585 + beta1 * (r-1)
    muS[r<1] <- beta0
    1 - pnorm (q=cutval, mean=muS-muN, sd=sqrt(sdS^2+sdN^2))
}
HHN <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]
    1 - exp(-lambda0 * exp (-r^2 / 2 / sigma^2))
}
HHR <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    1 - exp(-lambda0 * ( 1 - exp (-(r / sigma)^-z)))
}
HEX <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]
    1 - exp(-lambda0 * exp (-r / sigma))
}
HAN <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
    lambda0 * exp (-(r-w)^2 / 2 / sigma^2)
}
HCG <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    lambda0 * pgamma(r, shape=z, scale=sigma/z, lower.tail = FALSE)
}

############################################################################################

gradient <- function (pars, fun, eps=0.001, ...)
## quick & dirty 2009 09 14
## used by plot.secr for delta method limits
{
  est <- pars
  g   <- pars
  for (i in 1:length(est))
  {
      temp     <- est[i]
      if (temp != 0.0) delta <- eps * abs(temp)
      else             delta <- eps
      est[i]  <- temp - delta
      fminus  <- fun (est, ...)
      est[i]  <- temp + delta
      fplus   <- fun (est, ...)
      g[i]    <- (fplus - fminus) / (2.0 * delta)
      est[i]  <- temp;
  }
  g
}
############################################################################################

distancetotrap <- function (X, traps) {
    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coor in col 2
    X <- matrix(unlist(X), ncol = 2)

    nxy <- nrow(X)
    ## 2011-10-14
    detecttype <- detector(traps)
    detecttype <- ifelse (is.null(detecttype), '', detecttype)
    if (detecttype %in% .localstuff$polydetectors) {
        ## approximate only
        traps <- split(traps, polyID(traps))
        trpi <- function (i, n=100) {
            intrp <- function (j) {
                tmp <- data.frame(traps[[i]][j:(j+1),])
                if (tmp$x[1] == tmp$x[2])
                    data.frame(x=rep(tmp$x[1], n),
                               y=seq(tmp$y[1], tmp$y[2], length=n))
                else
                    data.frame(approx(tmp, n=n))
            }
            tmp <- lapply(1:(nrow(traps[[i]])-1),intrp)
            do.call(rbind, tmp)
        }
        trps <- do.call(rbind, lapply(1:length(traps), trpi))
        trps <- matrix(unlist(trps), ncol = 2)
    }
    else
        ## 2015-10-18 added protection
        trps <- matrix(unlist(traps), ncol = 2)

    if (inherits(trps, 'SpatialPolygonsDataFrame')) {
        trps <- coordinates(trps@polygons[[1]]@Polygons[[1]])
        warning("using only first polygon of SpatialPolygonsDataFrame")
    }

    temp <- .C('nearest',  PACKAGE = 'secr',
        as.integer(nxy),
        as.double(X),
        as.integer(nrow(trps)),
        as.double(unlist(trps)),
        index = integer(nxy),
        distance = double(nxy)
    )
    if (detecttype %in% c('polygon', 'polygonX')) {
        inside <- lapply(traps, pointsInPolygon, xy=X)
        inside <- do.call(rbind, inside)
        temp$distance [apply(inside,2,any)] <- 0
    }
    temp$distance
}

nearesttrap <- function (X, traps) {
    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coord in col 2
    X <- matrix(unlist(X), ncol = 2)
    nxy <- nrow(X)
    if (inherits(traps, 'SpatialPolygonsDataFrame')) {
        traps <- coordinates(traps@polygons[[1]]@Polygons[[1]])
        warning("using only first polygon of SpatialPolygonsDataFrame")
    }
    temp <- .C('nearest',  PACKAGE = 'secr',
        as.integer(nxy),
        as.double(X),
        as.integer(nrow(traps)),
        as.double(unlist(traps)),
        index = integer(nxy),
        distance = double(nxy)
    )
    temp$index
}

transform <- function (x, link) {
  switch (link,
          identity = x,
          log = log(x),
          neglog = log(-x),
          logit = logit(x),
          odds = odds(x),
          sin = sine(x)
  )
}

untransform <- function (beta, link) {
  switch (link,
          identity = beta,
          log = exp(beta),
          neglog = -exp(beta),
          logit = invlogit(beta),
          odds = invodds(beta),
          sin = invsine(beta))
}

se.untransform <- function (beta, sebeta, link) {
  switch (link,
          identity = sebeta,
          log = exp(beta) * sqrt(exp(sebeta^2)-1),
          neglog = exp(beta) * sqrt(exp(sebeta^2)-1),
          logit = invlogit(beta) * (1-invlogit(beta)) * sebeta,
          sin = NA)         ####!!!!
}

# mlogit.untransform <- function (beta, mix) {
#     ## beta should include values for all classes (mixture components)
#     nmix <- max(mix)    ## assume zero-based
#     b <- beta[2:nmix]    ## 2010 02 26
#     pmix <- numeric(nmix)
#     pmix[2:nmix] <- exp(b) / (1+sum(exp(b)))
#     pmix[1] <- 1 - sum(pmix[2:nmix])
#     pmix[mix]   ## same length as input
# }

mlogit.untransform <- function (beta, latentmodel) {
    if (!missing(latentmodel)) {
        ## this old code disordered the returned values; 2013-10-28
        ## tmp <- split(beta, latentmodel)
        ## unlist(lapply(tmp, mlogit.untransform))
        for (i in unique(latentmodel))
            beta[latentmodel==i] <- mlogit.untransform(beta[latentmodel==i])
        beta
    }
    else {
        ## beta should include values for all classes (mixture components)
        nmix <- length(beta)
        if (sum(is.na(beta)) != 1) {
            ## replaced 2013-06-06
            ## stop ("require NA for a single reference class in mlogit.untransform")
            rep(NA, length(beta))
        }
        else {
            nonreference <- !is.na(beta)   # not reference class
            b <- beta[nonreference]
            pmix <- numeric(nmix)
            pmix[nonreference] <- exp(b) / (1+sum(exp(b)))
            pmix[!nonreference] <- 1 - sum(pmix[nonreference])
            pmix
        }
    }
}

clean.mlogit <- function(x) {
## bad line removed 2013-05-09
##    x <- x/sum(x)
    ## 2014-08-19 for robustness...
    if (is.na(x[2])) x[2] <- 1-x[1]
    x[1] <- NA   ## assumed reference class
    logit(mlogit.untransform(x))
}

mlogit <- function (x) {
    ## return the mlogit of an unscaled vector of positive values
    ## 2013-04-14
    logit(x/sum(x))
}

# vector version of transform()
Xtransform <- function (real, linkfn, varnames) {
  out <- real
  for (i in 1:length(real)) {
      vn <- varnames[i]
      out[i] <- switch (linkfn[[vn]],
                  identity = real[i],
                  log = log(real[i]),
                  neglog = log(-real[i]),
                  logit = logit(real[i]),
                  odds = odds(real[i]),
                  sin = sine(real[i]))
  }
  out
}
se.Xtransform <- function (real, sereal, linkfn, varnames) {
  out <- real
  for (i in 1:length(real)) {
      vn <- varnames[i]
      out[i] <- switch (linkfn[[vn]],
                  identity = sereal[i],
                  log = log((sereal[i]/real[i])^2 + 1)^0.5,
                  neglog = log((sereal[i]/-real[i])^2 + 1)^0.5,
                  logit = sereal[i] / real[i] / (1 - real[i]),
                  sin = NA)
  }
  out
}

# vector version of untransform()
Xuntransform <- function (beta, linkfn, varnames) {
  out <- beta
  for (i in 1:length(beta)) {
      vn <- varnames[i]
      out[i] <- switch (linkfn[[vn]],
                  identity = beta[i],
                  log = exp(beta[i]),
                  neglog = -exp(beta[i]),
                  logit = invlogit(beta[i]),
                  odds = invodds(beta[i]),
                  sin = invsine(beta[i]))
  }
  out
}

se.Xuntransform <- function (beta, sebeta, linkfn, varnames)
# Approximate translation of SE to untransformed scale
# Delta method cf Lebreton et al 1992 p 77
{
  out <- beta
  if (length(beta)!=length(sebeta))
      stop ("'beta' and 'sebeta' do not match")
  if (!all(varnames %in% names(linkfn)))
      stop ("'linkfn' component missing for at least one real variable")
  for (i in 1:length(beta)) {
      vn <- varnames[i]
      out[i] <- switch (linkfn[[vn]],
                  identity = sebeta[i],
                  log = exp(beta[i]) * sqrt(exp(sebeta[i]^2)-1),
                  neglog = exp(beta[i]) * sqrt(exp(sebeta[i]^2)-1),
                  logit = invlogit(beta[i]) * (1-invlogit(beta[i])) * sebeta[i],
                  sin = NA)         ####!!!!
  }
  out
}

## End of miscellaneous functions
############################################################################################

group.levels <- function (capthist, groups, sep='.') {
    if (inherits(capthist, 'list')) {
        temp <- lapply(capthist, group.levels, groups, sep)   ## sep added 2010 02 24
        sort(unique(unlist(temp)))  ## vector of global levels
    }
    else {
        if (is.null(groups)) 0
        else {
            ## 2014-12-04 error check precedes attempt to extract by groups
            if (!all(groups %in% names(covariates(capthist))))
                stop ("one or more grouping variables is missing ",
                      "from covariates(capthist)")
            temp <- as.data.frame(covariates(capthist)[,groups])
##            if (ncol(temp) != length(groups))
##                stop ("one or more grouping variables is missing ",
##                      "from covariates(capthist)")
            sort(levels(interaction(temp, drop=T, sep=sep)))  # omit null combinations, sort as with default of factor levels
        }
    }
}
############################################################################################

h.levels <- function (capthist, hcov, nmix) {
    ## determine the first nmix levels of a factor individual covariate
    if (is.null(hcov))
        as.character(1:nmix)
    else {
        if (ms(capthist)) {
            ## take first session as we can assume factor covariates have same levels in
            ## all sessions
            capthist <- capthist[[1]]
        }
        hcov <- covariates(capthist)[,hcov]
        if (!is.factor(hcov)) {
            warning("hcov was coerced to a factor")
            hcov <- factor(hcov)
        }
        levels(hcov)[1:nmix]
    }
}
############################################################################################

n.occasion <- function (capthist) {
## return the number of sampling occasions for each session in capthist
    if (inherits(capthist, 'list')) {
        sapply(capthist, n.occasion)
    }
    else {
        ncol(capthist)
    }
}

############################################################################################

group.factor <- function (capthist, groups, sep='.')
## convert a set of grouping factors to a single factor (g)
## levels common to all sessions
{
    if (inherits(capthist, 'list')) {
        temp <- lapply(capthist, group.factor, groups)  ## recursive call
        grouplevels <- group.levels(capthist, groups)
        if (length(grouplevels)<2)
            temp
        else
            # list; force shared factor levels on each component
            lapply (temp, factor, levels=grouplevels)
    }
    else {
        if (is.null(groups) | (length(groups)==0) )
            return (factor(rep(1, nrow(capthist))))
        temp <- as.data.frame(covariates(capthist)[,groups])
        if (ncol(temp) != length(groups))
            stop ("one or more grouping variables is missing from ",
                  "covariates(capthist)")
        temp <- interaction(temp, drop=T, sep=sep)  # omit null combinations
        temp
    }
}
############################################################################################

make.lookup <- function (tempmat) {

    ## should add something to protect make.lookup from bad data...
    nrw <- nrow(tempmat)
    ncl <- ncol(tempmat)
    nam <- colnames(tempmat)

    df <- is.data.frame(tempmat)
    if (df) {
       lev <- lapply(tempmat, levels)
       tempmat <- unlist(sapply(tempmat, as.numeric, simplify = FALSE))
    }
    dimnames(tempmat) <- NULL

    temp <- .C('makelookup', PACKAGE = 'secr',
        ## as.double(tempmat),
        ## 2014-09-09
        ## compare values rounded to 10 dp
        ## required because floating point orthogonal polynomials sometimes differ
        as.double(round(tempmat, 10)),
        as.integer(nrw),
        as.integer(ncl),
        unique = integer(1),
        y      = double(nrw * ncl),
        index  = integer(nrw),
        result = integer(1))

    if (temp$result != 0)
        stop ("error in external function 'makelookup'; ",
              "perhaps problem is too large")
    lookup <- matrix(temp$y[1:(ncl*temp$unique)], nrow = temp$unique, byrow = T)
    colnames(lookup) <- nam
    if (df) {
        lookup <- as.data.frame(lookup)
        ## restore factors
        for (i in 1: length(lev))
            if (!is.null(lev[[i]]))
                lookup[,i] <- factor(lev[[i]][lookup[,i]], levels = lev[[i]])
    }
    list (lookup=lookup, index=temp$index)
}
###############################################################################

## 2013-04-12, 2013-06-05
## Return an integer vector of class membership defined by a categorical
## individual covariate in a capthist object. Individuals of unknown
## class (including those with class exceeding nmix) are coded 1,
## others as (class number + 1). When no mixture is specified (nmix == 1)
## all are coded as unknown.
getknownclass <- function(capthist, nmix, hcov) {
    if (ms(capthist)) {
        lapply(capthist, getknownclass, nmix = nmix, hcov = hcov)
    }
    else {
        if ((nmix>1) & (!is.null(hcov))) {
            tmp <- as.numeric(factor(covariates(capthist)[,hcov])) + 1
            tmp[is.na(tmp) | (tmp>(nmix+1))] <- 1
            attr(tmp,'levels') <- levels(factor(covariates(capthist)
                                                [,hcov]))[1:nmix]
            tmp
        }
        else
            rep(1,nrow(capthist))
    }
}

###############################################################################

getnmix <- function (details) {
    if (is.null(details$nmix))
       1
    else
       details$nmix
}

###############################################################################

## expand beta parameter vector using template of 'fixed beta'
## fixed beta fb input is missing (NA) for estimated beta parameters
fullbeta <- function (beta, fb) {
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta  ## partial beta (varying only)
        beta <- fb             ## complete beta
    }
    beta
}
###############################################################################

## inflate a convex outline along all radii by linear factor 'rmult'
## 2013-06-15
inflate <- function (xy, rmult = 1) {
    xy <- as.matrix(xy)
    centre <- apply(xy, 2, mean)
    xy <- sweep(xy, MARGIN = 2, STATS = centre, FUN = '-')
    r <- apply(xy, 1, function(z) sqrt(sum(z^2)))
    theta <- atan2 (xy[,2], xy[,1])
    r <- r * rmult
    xy <- cbind(r * cos(theta), r * sin(theta))
    sweep(xy, MARGIN = 2, STATS = centre, FUN = '+')
}
###############################################################################
## moved from pdot.R 2013-11-09
getbinomN <- function (binomN, detectr) {
    if (detectr %in% .localstuff$countdetectors) {
        if (is.null(binomN))
            return(0)
        else if (binomN == 'usage')
            return(1)
        else
            return(binomN)
    }
    else
        return(1)
}
###############################################################################

## convert xylist attribute of a combined dataset into a standalone capthist
xy2CH <- function (CH, inflation = 1e-8) {
    xylist <- telemetryxy(CH)
    if (is.null(xylist))
        stop ("requires xylist attribute")
    n <- length(xylist)
    neach <- sapply(xylist, nrow)
    allxy <- do.call(rbind, xylist)
    trps <-  allxy[chull(allxy),]
    trps <- rbind(trps, trps[1,,drop=F])
    trps <- inflate(trps, 1 + inflation)  ## see also telemetry.R

    trps <- as.data.frame(trps)
    dimnames(trps) <- list(1:nrow(trps), c('x','y'))
    class(trps) <- c("traps","data.frame")
    detector(trps) <- "telemetry"
    polyID(trps) <- factor(rep(1,nrow(trps)))

    rown <- rep(names(xylist), neach)
    newCH <- array(neach, dim = c(n, 1, 1))
    attr(newCH, "detectedXY") <- allxy
    if (!is.null(covariates(CH))) {
        rowlookup <- match(names(xylist), rownames(CH))
        covariates(newCH) <- covariates(CH)[rowlookup,, drop=FALSE]
    }
    class(newCH) <- "capthist"
    traps(newCH) <- trps
    newCH
}
###############################################################################
## return coordinates from simple SpatialPolygons object
## returns list
getcoord <- function(obj){
    if (!inherits(obj, 'SpatialPolygons'))
        stop ("requires SpatialPolygons object")
    if (length(obj@polygons) > 1)
        warning ("using only first 'polygons'")
    Polygons <- obj@polygons[[1]]
    lapply(Polygons@Polygons, coordinates)
}


###############################################################################
## moved from mask.check.r 2014-08-28

inflatechull <- function (poly, r, ntheta = 60) {
    theta <- (2*pi) * (1:ntheta) / ntheta
    ## add supernumerary vertices
    temp  <- data.frame(x = apply(expand.grid(poly$x, r * cos(theta)),1,sum),
                   y = apply(expand.grid(poly$y, r * sin(theta)),1,sum))
    hull <- chull(temp)
    temp[c(hull,hull[1]), ]
}

###############################################################################
## used by sim.capthist to update telemetry boundary polygon 2013-11-21
refreshMCP <- function (CH, tol) {
    if (detector(traps(CH)) %in% c('polygon','polygonX','telemetry'))
        allxy <- xy(CH)
    else
        stop ("requires polygon or telemetry detector type")
    trps <-  allxy[chull(allxy),]
    trps <- inflatechull(trps, r = tol)   ## 2014-08-28
    class(trps) <- c("traps","data.frame")
    names(trps) <- c('x','y')
    detector(trps) <- detector(traps(CH))
    polyID(trps) <- rep(1,nrow(trps))
    traps(CH) <- trps
    CH
}
###############################################################################

## 2015-01-14 sess -> sessnum
## moved from derivedMS 2013-12-15
maskarea <- function (mask, sessnum = 1) {
    if (!ms(mask)) nrow(mask) * attr(mask,'area')
    else nrow(mask[[sessnum]]) * attr(mask[[sessnum]],'area')
}
## 2014-08-29
masklength <- function (mask, sessnum = 1) {
    if (!ms(mask)) nrow(mask) * attr(mask,'spacing')/1000
    else nrow(mask[[sessnum]]) * attr(mask[[sessnum]],'spacing')/1000
}
###############################################################################

complete.beta <- function (object) {
    fb <- object$details$fixedbeta
    if (!is.null(fb)) {
        nbeta <- length(fb)
        fb[is.na(fb)] <- object$fit$par
        beta <- fb
    }
    else {
        beta <- object$fit$par
    }
    beta
}
###############################################################################

complete.beta.vcv <- function (object) {
    fb <- object$details$fixedbeta
    if (!is.null(fb)) {
        nbeta <- length(fb)
        beta.vcv <- matrix(NA, nrow = nbeta, ncol = nbeta)
        beta.vcv[is.na(fb[row(beta.vcv)]) & is.na(fb[col(beta.vcv)])] <- object$beta.vcv
    }
    else {
        beta.vcv <- object$beta.vcv
    }
    beta.vcv
}
###############################################################################

smooths <- function (formula) {
    ## which terms in formula are smooths?
    ## returns logical vector
    labels <- attr(terms(formula), "term.labels")
    if (length(labels) > 0)
        sapply(labels, function (x) any(sapply(c("s\\(", "te\\(", "poly\\("), grepl, x)))
    else
        logical(0)
}
############################################################################################

polys <- function (formula) {
    ## which terms in formula are orthogonal polynomials?
    ## returns logical vector
    labels <- attr(terms(formula), "term.labels")
    if (length(labels) > 0)
        sapply(labels, grepl, pattern = "poly\\(")
    else
        logical(0)
}
############################################################################################

badsmooths <- function (formula) {
    ## does smooth specification conform to secr requirements?
    ## returns TRUE/FALSE
    labels <- attr(terms(formula), "term.labels")
    if (length(labels) > 0) {
        smoothterms <- sapply(labels, function (x)
                              any(sapply(c("s\\(", "te\\("), grepl, x)))
        labels <- labels[smoothterms]
        any(sapply(labels, function(x)
               grepl("s\\(", x) & !grepl("k =", x))) |
        any(sapply(labels, function(x)
               grepl("te\\(", x) & (!grepl("fx = TRUE", x) | !grepl("k =", x))))
    }
    else
        FALSE
}
############################################################################################

gamsetup <- function(formula, data, ...) {
    ## use 'session' column as dummy LHS so gam does not gag
    ## (cf secrgam:::make.density.design.matrix)
    ## session is always present in detection data, must be added for D
    if (is.null(data$session)) data$session <- rep(1,nrow(data))
    formula <- update.formula(formula, session ~ .)
    setup <- gam(formula, data = data, fit = FALSE, ...)
    colnames(setup$X) <- setup$term.names
    setup
}
############################################################################################

general.model.matrix <- function (formula, data, gamsmth = NULL, ...) {

    ## A function to compute the design matrix for the model in
    ## 'formula' given the data in 'data'. This is merely the result
    ## of model.matrix() unless 'formula' includes smooth terms -- s()
    ## or te() as described in mgcv ?formula.gam.

    ## If smooth terms are present then the matrix may be based on a
    ## previous gam setup (provided in the argument 'gamsmth') or
    ## computed de novo with gam(..., fit = FALSE)

    ## note 2014-08-24
    ## orthogonal polynomials e.g. poly(x,2) are handled by model.matrix,
    ## but the information needed for prediction at new data is not
    ## saved by secr.fit, so predict.secr generally fails with message
    ## "'degree' must be less than number of unique points"

    ##  head(eval(parse(text = attr(terms(~ poly(x,y, degree=2)),
    ##  'term.labels')[1]), env=possummask))

    ## 2014-08-24, 2014-09-09

    if (any(polys(formula)))
        stop ("orthogonal polynomials are temporarily blocked")  ## 2014-09-12
    if (any(smooths(formula))) {
        if (is.null(gamsmth)) {
            ## setup knots etc from scratch
            gamsetup(formula, data, ...)$X
        }
        else {
            ## fool predict.gam into generating the necessary
            ## predictor matrix from previous setup
            class (gamsmth) <- 'gam'
            gamsmth$coefficients <- rep(NA, ncol(gamsmth$X))
            mat <- mgcv::predict.gam(gamsmth, newdata = data, type = 'lpmatrix')
            colnames(mat) <- colnames(gamsmth$X)
            mat
        }
    }
    else {
        ## model.matrix(formula, data, ...)
        model.matrix(formula, data)
    }
}
###############################################################################

secr.lpredictor <- function (formula, newdata, indx, beta, field, beta.vcv=NULL,
                             smoothsetup = NULL) {
    ## form linear predictor for a single 'real' parameter
    ## smoothsetup should be provided whenever newdata differs from
    ## data used to fit model and the model includes smooths from gam
    vars <- all.vars(formula)
    ## improved message 2015-01-29
    OK <- vars %in% names(newdata)
    if (any(!OK)) {
        missingvars <- paste(vars[!OK], collapse = ', ')
        if (sum(!OK) == 1)
            stop ("model covariate ", missingvars, " not found in 'newdata'")
        else
            stop ("model covariates ", missingvars, " not found in 'newdata'")
    }
    newdata <- as.data.frame(newdata)
    lpred <- matrix(ncol = 2, nrow = nrow(newdata), dimnames = list(NULL,c('estimate','se')))

    ## 2014-08-19
    ## mat <- model.matrix(model, data = newdata)

    mat <- general.model.matrix(formula, data = newdata, gamsmth = smoothsetup)
    if (nrow(mat) < nrow(newdata))
        warning ("missing values in predictors?")

    ## drop pmix beta0 column from design matrix (always zero)
    nmix <- 1
    if (field=='pmix') {
        mat <- mat[,-1,drop=FALSE]
        if ('h2' %in% names(newdata)) nmix <- 2
        if ('h3' %in% names(newdata)) nmix <- 3
    }
    lpred[,1] <- mat %*% beta[indx]

    ## 2015-09-30 new code for pmix based on old fixpmix function in utility.R
    ## deals with mlogit link always used by pmix
    if ((nmix > 1) & (field == 'pmix')) {
        temp <- matrix(lpred[,1], ncol = nmix)
        if (nmix==2) temp[,newdata[,'h2']] <- lpred[,1]
        if (nmix==3) temp[,newdata[,'h3']] <- lpred[,1]
        temp2 <- apply(temp, 1, clean.mlogit)
        lpred[,1] <- as.numeric(t(temp2))
        if (nmix==2) {
            h2.1 <- as.numeric(newdata$h2)==1
            h2.2 <- as.numeric(newdata$h2)==2
            lpred[h2.1,2] <- lpred[h2.2,2]
        }
        else
            lpred[,2] <- rep(NA, nrow(lpred))   ## don't know how
    }
    ## 2015-09-30 end of new code

    if (is.null(beta.vcv) | (any(is.na(beta[indx])))) return ( cbind(newdata,lpred) )
    else {
        vcv <- beta.vcv[indx,indx, drop = FALSE]
        nrw <- nrow(mat)
        vcv <- apply(expand.grid(1:nrw, 1:nrw), 1, function(ij)
            mat[ij[1],, drop=F] %*% vcv %*% t(mat[ij[2],, drop=F]))  # link scale
        vcv <- matrix (vcv, nrow = nrw)
        ## 2015-09-30
        if (field=='pmix') {
            if (nmix==2)
                vcv[h2.1,h2.1] <- vcv[h2.2,h2.2]
            else
                vcv[,] <- NA
        }
        lpred[,2] <- diag(vcv)^0.5
        temp <- cbind(newdata,lpred)
        attr(temp, 'vcv') <- vcv
        return(temp)
    }
}

############################################################################################

## 2014-10-06
edist <- function (xy1, xy2) {
  nr <- nrow(xy1)
  nc <- nrow(xy2)
  x1 <- matrix(xy1[,1], nr, nc)
  x2 <- matrix(xy2[,1], nr, nc, byrow=T)
  y1 <- matrix(xy1[,2], nr, nc)
  y2 <- matrix(xy2[,2], nr, nc, byrow=T)
  sqrt((x1-x2)^2 + (y1-y2)^2)
}

############################################################################################

## 2015-01-14
## least cost paths from mask including barriers to movement
## use edist for equivalent Euclidean distances
## requires raster package

nedist <- function (xy1, xy2, mask, inf = Inf, ...) {
    newargs <- list(...)
    if (missing(mask)) mask <- xy2
    noneuc <- covariates(mask)$noneuc
    if (is.null(noneuc)) noneuc <- rep(1, nrow(mask))
    defaultargs <- list(transitionFunction = mean, directions = 16)
    args <- replace(defaultargs, names(newargs), newargs)
    args$x <- raster(mask, values = noneuc)
    if (requireNamespace('gdistance')) {    ## 2015-01-23
        tr <- do.call(gdistance::transition, args)
        tr <- gdistance::geoCorrection(tr, type = "c", multpl = FALSE)
        out <- gdistance::costDistance(tr, as.matrix(xy1), as.matrix(xy2))
    }
    else stop ("package gdistance is required for nedist")
    if (is.finite(inf)) out[!is.finite(out)] <- inf
    out
}

############################################################################################

getcellsize <- function (mask) {
    if (inherits(mask, 'linearmask'))
        cell <- attr(mask, 'spacing') / 1000  ## per km
    else
        cell <- attr(mask, 'area')            ## per ha
    if (is.null(cell))
        stop ("mask lacks valid cell size (area or spacing)")
    cell
}
############################################################################################

## 2014-10-25
## intercept and fix certain ,models with bad defaults
updatemodel <- function (model, detectfn, detectfns, oldvar, newvar) {
    if (detectfn %in% detectfns) {
        for (i in 1:length(oldvar)) {
            if (oldvar[i] %in% names(model)) {
                names(model)[names(model) == oldvar[i]] <- newvar[i]
                warning ("replacing ", oldvar[i], " by ", newvar[i], " in model for detectfn ",
                         detectfn)
            }
        }
    }
    model
}
############################################################################################

## Manually remove some mask points

deleteMaskPoints <- function (mask, onebyone = TRUE, add = FALSE, poly = NULL,
                              poly.habitat = FALSE, ...) {
    ## interface does not work properly in RStudio

    if (ms(mask)) {         ## a list of mask objects
        if (inherits(poly, 'list') & (!is.data.frame(poly)))
            stop ("lists of polygons not implemented in 'make.mask'")
        temp <- lapply (mask, deleteMaskPoints, onebyone = onebyone, add = add,
                        poly = poly, poly.habitat = poly.habitat, ...)
        class (temp) <- c('list', 'mask')
        temp
    }
    else {
        plot(mask, add = add, ...)
        if (!is.null(poly)) {
            SPDF <- inherits(poly, "SpatialPolygonsDataFrame")
            if (!SPDF) {
                poly <- matrix(unlist(poly), ncol = 2)
                poly <- rbind (poly, poly[1,])  # force closure of poly
            }
            if (poly.habitat)
                pointstodrop <- (1:nrow(mask))[!pointsInPolygon(mask, poly)]
            else
                pointstodrop <- (1:nrow(mask))[pointsInPolygon(mask, poly)]
        }
        else if (onebyone) {
            cat ('Click to select points; right-click to stop\n')
            flush.console()
            xy <- locator(type = 'p', pch=1, col='red')
            pointstodrop <- if (length(xy$x)==0)
                numeric(0)
            else
                nearesttrap(xy, mask)
        }
        else {
            cat ('Click to select polygon vertices; right-click to stop\n')
            flush.console()
            xy <- locator(type = 'l', col='red')
            xy <- as.data.frame(xy)
            xy <- rbind(xy, xy[1,])
            if (poly.habitat)
                pointstodrop <- (1:nrow(mask))[!pointsInPolygon(mask, xy)]
            else
                pointstodrop <- (1:nrow(mask))[pointsInPolygon(mask, xy)]
        }
        npts <- length(pointstodrop)
        if (npts>0) {
            points(mask[pointstodrop,], pch = 16, col = 'red')
            if(.Platform$OS.type == "windows") {
                pl <- if (npts>1) 's' else ''
                msg <- paste ('Delete ', npts, ' red point',pl, '?', sep='')
                response <-  utils::winDialog(type = "okcancel", msg)
            } else {
                response <- 'OK'
            }
            if (response == 'OK') {
                mask <- subset(mask, -pointstodrop)
            if (npts==1)
                message("1 point deleted")
            else
                message(npts, " points deleted")
            }
        else
            message ("point(s) not deleted")
        }
        else
            message ("no points to delete")
        plot(mask, col='green')
        mask
    }
}
############################################################################################

nparameters <- function (object) {
    Npar <- max(unlist(object$parindx))
    Npar <- Npar + length(object$details$miscparm)
    ## allow for fixed beta parameters
    if (!is.null(object$details$fixedbeta))
        Npar <- Npar - sum(!is.na(object$details$fixedbeta))
    Npar
}
############################################################################################

mapbeta <- function (parindx0, parindx1, beta0, betaindex)

    ## Extend beta vector from simple model (beta0) to a more complex (i.e. general)
    ## model, inserting neutral values (zero) as required.
    ## For each real parameter, a 1:1 match is assumed between
    ## beta values until all beta values from the simpler model are
    ## used up. THIS ASSUMPTION MAY NOT BE JUSTIFIED.
    ## betaindex is a user-controlled alternative.

{
    ## list of zeroed vectors, one per real parameter
    beta1 <- lapply(parindx1, function (x) {x[]<-0; x})

    if (!is.null(betaindex)) {
        beta1 <- unlist(beta1)
        if (sum(betaindex>0) != length(beta0))
            stop ("invalid 'betaindex'")
        beta1[betaindex] <- beta0
        beta1
    }
    else {
        ## indx is within-parameter rather than absolute index
        ## for each _original_ real parameter
        indx <- lapply(parindx0, function(x) x-x[1]+1)
        ## for (j in 1:length(beta1))
        ## improved replace by name2015-11-17
        for (j in names(beta1)) {
            if (j %in% names(beta0))
                beta1[[j]][indx[[j]]] <- beta0[parindx0[[j]]]
        }
        unlist(beta1)
    }
}
############################################################################################

xyinpoly <- function (xy, trps) {
    ptinside <- function (i,k) {
        ## is point i inside poly k?
        polyxy <- as.matrix(lxy[[k]])
        polyxy <- rbind(polyxy, polyxy[1,])   ## close 2014-08-28
        nr <- nrow(polyxy)
        temp <- .C('inside',  PACKAGE = 'secr',
                   as.double (xy[i,]),
                   as.integer (0),
                   as.integer (nr-1),
                   as.integer (nr),
                   as.double (polyxy),
                   result = integer(1))
        as.logical(temp$result)
    }
    lxy <- split (trps, levels(polyID(trps)))
    
    firstinside <- function (i) {
        frstk <- 0
        for (k in 1:length(lxy)) {
            if (ptinside(i,k)) {
                frstk <- k
                break
            }
        }
        frstk
    }
    sapply(1:nrow(xy), firstinside)
}
############################################################################################
##
## settings for mark-resight
markresight <- function (capthist, mask, CL, fixed, chat, sessnum) {
    markocc <- markocc(traps(capthist))
    s <- ncol(capthist)
    m <- nrow(mask)
    if (is.null(markocc)) {
        markocc <- rep(1, s)
        Tu <- c(-1,0)
        Tm <- c(-1,0)
        allsighting <- FALSE
    }
    else {
        allsighting <- !any(markocc>0)

        Tu <- if (CL) NULL else Tu(capthist)
        Tm <- Tm(capthist)
        if (!is.null(fixed$pID)) {
            if (fixed$pID == 1) Tm <- NULL
        }
        if (!any(markocc==0))
            Tm <- NULL
        ## second element is number of values; '1' indicates summed counts
        Tu <- if (is.null(Tu)) c(-1,0) else c(length(unlist(Tu)), Tu)
        Tm <- if (is.null(Tm)) c(-1,0) else c(length(unlist(Tm)), Tm)
    }
    if (!is.null(chat)) {
        if (is.matrix(chat))
            chat <- chat[sessnum,]
        else {
            chat <- unlist(chat)
            if (length(chat)==1) {
                chat <- c(chat,1)
                warning("assuming chat 1.0 for Tm")
            }
            if (chat[1]<1) {
                warning("chat(Tu) = ", chat[1], ", setting to 1.0")
            }
            if (chat[2]<1) {
                warning("chat(Tm) = ", chat[2], ", setting to 1.0")
            }
            chat <- pmax(chat, 1)
        }
    }
    else
        chat <- c(1,1)


    pi.mask <- -1      ## signals pimask not used
    if (allsighting) {
        ## pi.mask is Pr(marked animal is from pixel m)
        ## i.e. pdf(x) * area
        pi.mask <- rep(1/nrow(mask), nrow(mask))
        if (!is.null(maskcov <- covariates(mask))) {
            if ('marking' %in% names (maskcov)) {
                if (any(is.na(maskcov$marking)) | any (maskcov$marking<0))
                    stop ("invalid marking covariate in mask")
                pi.mask <- maskcov$marking / sum (maskcov$marking)
            }
        }
    }

    list(markocc = markocc, Tu = Tu, Tm = Tm, allsighting = allsighting, chat = chat,
         pi.mask = pi.mask)
}

############################################################################################

addzerodf <- function (df, oldCH, sess) {
    ## add dummy detection records to dataframe for 'all-zero' case
    ## that arises in sighting-only mark-resight with known marks
    allzero <- apply(oldCH,1,sum)==0
    naz <- sum(allzero)
    if (naz > 0) {
        df0 <- expand.grid(newID = rownames(oldCH)[allzero], newocc = NA, 
                           newtrap = trap(oldCH)[1], alive = TRUE, sess = sess, 
                           stringsAsFactors = FALSE)
        df <- rbind(df,df0)
        if (!is.null(xy(oldCH))) {
            df$x <- c(xy(oldCH)$x, rep(NA, naz))
            df$y <- c(xy(oldCH)$y, rep(NA, naz))
        }
        if (!is.null(signal(oldCH)))  {
            df$signal <- c(signal(oldCH), rep(NA, naz))
        }
    }
    df
}
############################################################################################

## including pre-marked animals never sighted
## cov is optional dataframe of covariates
addzeroCH <- function (CH, nzero, cov = NULL) {
    if (nzero == 0)
        return(CH)
    else {
        
        nc <- nrow(CH)
        chdim <- dim(CH)
        chdim[1] <- nzero
        extra <- array(0, dim=chdim)
        dimnames(extra) <- c(list(paste('Z', 1:nzero, sep='')), dimnames(CH)[2:3])
        CH2 <- abind(CH, extra, along = 1)
        class(CH2) <- 'capthist'
        traps(CH2) <- traps(CH)
        xy(CH2) <- xy(CH)  ## order is not affected by adding zero histories
        if (!is.null(covariates(CH)) & (nrow(CH)>0)) {
            if (is.null(cov)) {
                cov <- covariates(CH)[rep(1,nzero),]
                cov[,] <- NA   ## covariates are unknown
            }
            covariates(CH2) <- rbind(covariates(CH), cov[1:nzero,])
        }
        ## ... and other essential attributes?
        CH2
    }
}
############################################################################################

