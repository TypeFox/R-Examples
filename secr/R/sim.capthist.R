## package 'secr'
## sim.capthist.R
## simulate capture histories
## 2009 10 08 sim.resight
## 2010 07 01 allow alphanumeric detection functions
## 2010 10 09 annular normal and cumulative lognormal detection functions
## 2011 03 19 allow zero detections
## 2011 03 27 multiple sessions
## 2011 09 09 p.available; session numbers
## 2011 09 30 presence, unmarked treated as special case of proximity
## 2013 04 04 extend to k-specific g0
## 2013 05 10 extend to general learned response (recapfactor)
## 2013 06 28 explicit entry points to C code (not simfunctionname)
## 2014-08-28 revamped for userdist and distmat (pre-computed distance matrix)
## 2015-02-23 warning when ms(popn) and renumber = TRUE
## 2015-04-01 session-specific detection
## 2015-04-13 sim.capthist passes nsessions to sim.popn
## 2015-11-02 sim.resight, expands re-written
## 2015-12-15 sim.resight delivers for unresolved sightings (markocc = -1)
###############################################################################

expands <- function (x, s, default = 1) {
    if (is.null(x))
        rep(default, s)
    else {
        y <- numeric(s)
        y[] <- x
        y
    }
}

expandsk <- function (x, s, k, default = 1) {
    if (is.null(x))
        matrix(default, nrow = s, ncol = k)
    else
        matrix(x, nrow = s, ncol = k)
}

sim.capthist <- function (
    traps,
    popn = list(D = 5, buffer = 100, Ndist = 'poisson'),
    detectfn = 0,
    detectpar = list(),
    noccasions = 5,
    nsessions = 1,
    binomN = NULL,
    exactN = NULL,
    p.available = 1,
    renumber = TRUE,
    seed = NULL,
    maxperpoly = 100,
    chulltol = 0.001,
    userdist = NULL,
    savepopn = FALSE
    )

#
# Simulate detection of known population
#

## A note on sort order of records  2009 12 01

## Simulation routines return the primary detection data in 'ski' order,
## meaning that occasion (s) changes fastest and individual (i) slowest -
## this allows efficient sequential output as new animals are detected.
## In R the order is changed to 'isk' as this is pictorially more natural.

## Secondary data (xy locations, signal strength) are generated in the order
## 'kis' (detector (k) changing fastest), because the simulation routines use -
## for (s=0; s<*ss; s++)
##   for (i=0; i<*N; i++)
##     for (k=0; k<*kk; k++) {
##     ...
##     }
## or similar loops.
## For consistency, all data are sorted to 'isk' order before returning to R.
## Secondary simulated data are re-sorted into 'ksi' order by creating the
## index 'start' within C as required for prwipolygon & prwitransect.
## The functions 'animalID', 'occasion', and 'trap' are the safest way to
## retrieve values for detections in isk order.

{
    poplist <- inherits(popn,'popn') & inherits(popn,'list')
    ##------------------------------------------------------------------------
    ## multi-session loop
    if (poplist | (nsessions > 1)) {

        if (poplist & (nsessions>1) & (length(popn) != nsessions))
            stop ("incompatible use of popn list and nsessions>1")

        nsessions <- ifelse (poplist, length(popn), nsessions)
        if (ms(traps) & (length(traps) != nsessions))
            stop ("incompatible use of traps list and nsessions>1")

        output <- vector(nsessions, mode='list')
        nocc <- numeric(nsessions)
        nocc[] <- noccasions

        ##################
        ## set random seed
        ## copied from simulate.lm
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1)
        if (is.null(seed))
            RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        else {
            R.seed <- get(".Random.seed", envir = .GlobalEnv)
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }
        #####################################################################
        ## generate population list if not provided
        if (!inherits(popn,'popn')) {
            popn <- replacedefaults(list(D = 5, buffer = 100,
                Ndist = 'poisson'), popn)
            ## will fail with multiple traps
            popn <- sim.popn (popn$D, core = traps, buffer = popn$buffer,
                covariates = NULL, Ndist = popn$Ndist, nsessions = nsessions)
            poplist <- TRUE
        }

        #####################################################################
        ## availability preliminaries
        if (poplist) {
            if (any(p.available != 1))
                warning ("incomplete availability not implemented ",
                         "for population lists")
            if(renumber)
                warning ("typically use renumber = FALSE for multiple sessions")
        }
        else {
            if (!(length(p.available) %in% 1:2))
                stop ("p.available must be vector of 1 or 2 probabilities")
            availability <- 'random'
            if (length(p.available) == 1)
                ## random temporary emigration
                available <- runif(nrow(popn)) < p.available
            else {
                ## Markovian temporary emigration
                availability <- 'Markov'
                equilibrium.p <- (p.available[2] / (1-p.available[1]+p.available[2]))
                available <- runif(nrow(popn)) < equilibrium.p
            }
        }

        #####################################################################
        ## session-specific detection 2015-04-01
        if (detectfn %in% c(0:7, 14:18)) {
            df0name <- if (detectfn %in% (0:7)) 'g0' else 'lambda0'
            dfzname <- if (detectfn %in% (5:6)) 'x' else 'z'
            df0 <- expands(detectpar[[df0name]], nsessions, default = NULL)
            sigma <- expands(detectpar[['sigma']], nsessions, default = NULL)
            dfz <- expands(detectpar[[dfzname]], nsessions, default = NULL)
        }

        ########################################################################
        ## loop over sessions
        for (t in 1:nsessions) {
            ## if (R > 1)  ## modified 2015-10-29

            if (poplist)
                temppop <- popn[[t]]
            else {
                temppop <- subset(popn, available)
                ##-------------------------------------------------------------
                ## update availability in preparation for next session
                if (availability == 'random') {
                    available <- runif(nrow(popn)) < p.available
                }
                else {
                    p.vect <- ifelse(available, p.available[1], p.available[2])
                    available <- runif(nrow(popn)) < p.vect
                }
                #--------------------------------------------------------------
            }

            ## select session-specific traps if necessary
            trps <- if (ms(traps)) traps[[t]] else traps

            ## session-specific detection parameters
            if (detectfn %in% c(0:7, 14:18)) {
                detectpar[[df0name]] <- df0[t]
                detectpar[['sigma']] <- sigma[t]
                detectpar[[dfzname]] <- dfz[t]
            }

            ## recursive call for each session
            ## we have set the seed above
            output[[t]] <- sim.capthist(trps, temppop, detectfn, detectpar,
                noccasions = nocc[t], nsessions = 1,
                binomN = binomN, exactN = exactN, p.available = 1,
                renumber = renumber, seed = NULL, maxperpoly)

            session( output[[t]] ) <- t   ## added 2011-09-09
        }
        ########################################################################

        class(output) <- c('list','capthist')
        names(output) <- 1:nsessions
        output
    }   ## end of multi-session

    else {

        ##################
        ## set random seed
        ## copied from simulate.lm
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1)
        if (is.null(seed))
            RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        else {
            R.seed <- get(".Random.seed", envir = .GlobalEnv)
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }
        ##################

        if (is.null(detector(traps)))
            stop ("'traps' lacks detector type")

        usge <- usage(traps)
        if (is.null(usge))
            usge <- matrix (1, nrow = ndetector(traps), ncol = noccasions)
        else {
            if (nrow(usge) != ndetector(traps))
                stop ("invalid usage matrix; number of rows ",
                      "must match number of detectors")
            if (ncol(usge) != noccasions) {
                noccasions <- ncol(usge)
                warning ("'noccasions' does not match usage ",
                         "attribute of 'traps'; ignored")
            }
        }


        validatepar <- function (x, xrange) {
            xname <- deparse(substitute(x))
            if (is.null(x))
                stop ("no value for ", xname)
            if (any(is.na(x)))
                stop ("NA is not a valid value for ", xname)
            if (any(x < xrange[1]))
                warning ("value for ", xname, " is less than minimum ",
                    xrange[1])
            if (any(x > xrange[2]))
                warning ("value for ", xname, " is greater than maximum ",
                    xrange[2])
        }

        ## added 2010-07-01
        if (is.character(detectfn))
            detectfn <- detectionfunctionnumber(detectfn)
        if (detector(traps) %in% c('signal')) {
            if (!(detectfn %in% c(10,11))) {
                warning ("forcing detection function = 10 for signal detectors")
                detectfn <- 10
            }
        }
        else if (detector(traps) %in% c('signalnoise')) {
            if (!(detectfn %in% c(12,13))) {
                warning ("forcing detection function = 12 for signalnoise detectors")
                detectfn <- 12
            }
        }

        ## Detection function parameters

        ##    0  halfnormal
        ##    1  hazard rate
        ##    2  exponential
        ##    3  compound halfnormal
        ##    4  uniform
        ##    5  w-exponential
        ##    6  annular normal
        ##    7  cumulative lognormal
        ##    8  cumulative gamma??
        ##    9  binary signal strength (b0 = (beta0-c)/sdS, b1 = beta1/sdS)
        ##   10  signal strength (signal detectors only)
        ##   11  signal strength with spherical spreading (signal detectors only)
        ##   12  signal-noise (signalnoise detectors only)
        ##   13  signal-noise with spherical spreading (signalnoise detectors only)
        ##   14  hazard halfnormal
        ##   15  hazard hazard-rate
        ##   16  hazard exponential
        ##   17  annular normal
        ##   18  cumulative gamma??

        ## 2012-12-23
        if (!is.null(usage(traps))) {
            if (!is.null(noccasions)) {
                if (noccasions != ncol(usage(traps)))
                    warning ("noccasions differs from ncol(usage); using latter")
            }
            noccasions <- ncol(usage(traps))
        }

        ## extended for uniform (detectfn=4) 2010-06-13
        ## extended for hazard XX detection functions 2013-04-24
        if (detectfn %in% c(0:4))  defaults <- list(g0 = 0.2, sigma = 25, z = 1)
        if (detectfn %in% c(5))    defaults <- list(g0 = 0.2, sigma = 25, w = 10)
        if (detectfn %in% c(6))    defaults <- list(g0 = 0.2, sigma = 25, w = 10)
        if (detectfn %in% c(7,8))    defaults <- list(g0 = 0.2, sigma = 25, z = 5)
        if (detectfn %in% c(9))  defaults <- list(b0 = 1, b1=-0.1, cutval = 60,
            tx = 'identity')
        if (detectfn %in% c(10,11))  defaults <- list(beta0 = 90, beta1=-0.2,
            sdS = 2, cutval = 60, sdM = 0, tx = 'identity')
        if (detectfn %in% c(14:16))  defaults <- list(lambda0 = 0.2, sigma = 25, z = 1)
        if (detectfn %in% c(17))  defaults <- list(lambda0 = 0.2, sigma = 25, w =10)
        if (detectfn %in% c(18))  defaults <- list(lambda0 = 0.2, sigma = 25, z = 5)

        if (detectfn %in% c(12,13))  defaults <- list(beta0 = 90, beta1=-0.2,
            sdS = 2, cutval = 10, muN = 40, sdN = 2, sdM = 0, tx = 'identity')
        else defaults <- c(defaults, list(truncate = 1e+10, recapfactor = 1.0))

        defaults$binomN <- 0
        if (detector(traps) == 'proximity') defaults$binomN <- 1
        if (detector(traps) %in% c('signal','signalnoise')) defaults$binomN <- 1
        if (!is.null(binomN)) {
            if ((detector(traps) == 'count') & (tolower(binomN) == 'usage'))
                binomN <- 1   ## code for 'binomial size from usage' 2012-12-22
            detectpar$binomN <- binomN
        }

        detectpar <- replacedefaults(defaults, detectpar)

        if (detectfn %in% c(0:7, 14:18)) {
            # g0    <- expands(detectpar$g0, noccasions)
            # changed to matrix for s- and k-specific g0 2013-04-04
            if (detectfn %in% (0:7)) {
                g0    <- expandsk(detectpar$g0, s = noccasions, k = ndetector(traps))
                if (detector(traps) %in% .localstuff$countdetectors) {
                    if (detectpar$binomN == 0)
                        validatepar(g0, c(0,Inf))  ## Poisson lambda
                    else
                        validatepar(g0, c(0,1))    ## Binomial p
                }
                df0 <- g0
            }
            else  {
                lambda0 <- expandsk(detectpar$lambda0, s = noccasions, k = ndetector(traps))
                validatepar(lambda0, c(0,Inf))  ## Poisson lambda
                df0 <- lambda0
            }

            sigma <- expands(detectpar$sigma, noccasions)
            z     <- expands(ifelse(detectfn %in% c(5,6),
                detectpar$w, detectpar$z), noccasions)
            validatepar(sigma, c(1e-10,Inf))
            validatepar(z, c(0,Inf))
        }

        # Acoustic detection function parameters
        if (detectfn %in% c(10,11,12,13)) {
            tx <- detectpar$tx
            cutval <- detectpar$cutval
            sdM <- detectpar$sdM
            beta0 <- expands(detectpar$beta0, noccasions)
            beta1 <- expands(detectpar$beta1, noccasions)
            sdS   <- expands(detectpar$sdS, noccasions)
            muN   <- expands(detectpar$muN, noccasions)
            sdN   <- expands(detectpar$sdN, noccasions)
            validatepar(beta0, c(-Inf,Inf))
            validatepar(beta1, c(-Inf,Inf))
            validatepar(sdS, c(0,Inf))
            validatepar(sdN, c(0,Inf))
        }
        else if (detectfn %in% c(9)) {
            df0 <- expandsk(detectpar$b0, s = noccasions, k = ndetector(traps))
            sigma <- expands(detectpar$b1, noccasions)
            z <- 0
            cutval <- detectpar$cutval
        }
        else {
            cutval <- NULL
            truncate <- ifelse(is.null(detectpar$truncate),
                1e+10, detectpar$truncate)
            validatepar(truncate, c(1e-10, Inf)) ## must be positive
        }

        # Allow general learned response for traps
        if (detector(traps) %in% c('single','multi')) {
            rfl <- length(detectpar$recapfactor)
            if (!(rfl %in% 1:2))
                stop ("invalid recapfactor")
            if (rfl==1)
                detectpar$recapfactor <- c(detectpar$recapfactor,1)
            df0 <- c(df0, df0 * detectpar$recapfactor[1])
            sigma <- c(sigma, sigma * detectpar$recapfactor[2])
        }
        else {
            if (any(detectpar$recapfactor != 1.0))
                stop("learned response available only for 'single', 'multi' detectors")
        }

        #-----------------------------------------------------------------------------#
        if (!inherits(popn,'popn')) # generate if not provided
        {
            popn <- replacedefaults(list(D = 5, buffer = 100,
                Ndist = 'poisson'), popn)
            popn <- sim.popn (popn$D, core = traps, buffer = popn$buffer,
                covariates = NULL, Ndist = popn$Ndist)
        }
        #-----------------------------------------------------------------------------#

        N <- nrow(popn)
        animals <- as.matrix(popn)
        k <- nrow(traps)
        simfunctionname <- paste('trapping', detector(traps), sep='')

        ##-----------------------------------------------------------------
        ## user-provided distances

        if (is.null(userdist))
            distmat <- -1
        else {
            ## 2014-10-09
            ## move towards IHP/k simulations
            ## mask is NULL unless IHP or linear
            mask <- attr(popn, 'mask')
            ## ASSUME MASK HAS NONEUC REQUIREMENTS 2014-10-30
            distmat <- valid.userdist(userdist,
                                  detector(traps),
                                  xy1 = traps,
                                  xy2 = animals,
                                  mask = mask)
        }
        #-----------------------------------------------------------------

        if (detector(traps) %in% c('single','multi')) {
            if (detector(traps) == 'single')
                temp <- .C('trappingsingle', PACKAGE = 'secr',
                           as.double(df0),
                           as.double(sigma),
                           as.double(z),
                           as.integer(noccasions),
                           as.integer(k),
                           as.integer(N),
                           as.double(animals),
                           as.double(unlist(traps)),
                           as.double(distmat),
                           as.double(usge),
                           as.integer(detectfn),
                           as.double(truncate^2),
                           n = integer(1),
                           caught = integer(N),
                           value=integer(N*noccasions),
                           resultcode = integer(1)
                           )
            else
                temp <- .C('trappingmulti', PACKAGE = 'secr',
                           as.double(df0),
                           as.double(sigma),
                           as.double(z),
                           as.integer(noccasions),
                           as.integer(k),
                           as.integer(N),
                           as.double(animals),
                           as.double(unlist(traps)),
                           as.double(distmat),
                           as.double(usge),
                           as.integer(detectfn),
                           as.double(truncate^2),
                           n = integer(1),
                           caught = integer(N),
                           value=integer(N*noccasions),
                           resultcode = integer(1)
                           )
            if (temp$resultcode != 0)
                stop ("call to '", simfunctionname, "' failed")
            w <- matrix(ncol = temp$n, nrow = noccasions, dimnames =
                 list(1:noccasions, NULL))
            if (temp$n > 0) w[,] <- temp$value[1:(temp$n*noccasions)]
            w <- t(w)
        }

        else
        if (detector(traps) %in% c('polygonX','transectX')) {
            if (detector(traps) == 'polygonX') {
                nk <- length(levels(polyID(traps)))
                k <- table(polyID(traps))
                temp <- .C('trappingpolygonX', PACKAGE = 'secr',
                           as.double(df0),
                           as.double(sigma),
                           as.double(z),
                           as.integer(noccasions),
                           as.integer(nk),
                           as.integer(k),
                           as.integer(N),
                           as.double(animals),
                           as.double(unlist(traps)),
                           as.double(usge),
                           as.integer(detectfn),
                           as.double(truncate^2),
                           n = integer(1),
                           caught = integer(N),
                           detectedXY = double (N*noccasions),
                           value = integer(N*noccasions),
                           resultcode = integer(1)
                           )
            }
            else {
                nk <- length(levels(transectID(traps)))
                k <- table(transectID(traps))
                temp <- .C('trappingtransectX', PACKAGE = 'secr',
                           as.double(df0),
                           as.double(sigma),
                           as.double(z),
                           as.integer(noccasions),
                           as.integer(nk),
                           as.integer(k),
                           as.integer(N),
                           as.double(animals),
                           as.double(unlist(traps)),
                           as.double(usge),
                           as.integer(detectfn),
                           as.double(truncate^2),
                           n = integer(1),
                           caught = integer(N),
                           detectedXY = double (N*noccasions),
                           value = integer(N*noccasions),
                           resultcode = integer(1)
                           )
            }
            if (temp$resultcode != 0)
                stop ("call to '", simfunctionname, "' failed")

            w <- matrix(ncol = temp$n, nrow = noccasions, dimnames =
                list(1:noccasions, NULL))
            if (temp$n > 0) {
                w[,] <- temp$value[1:(temp$n*noccasions)]
            }
            w <- t(w)

            if (temp$n > 0) {
                ## put XY coordinates in attribute
                nd <- sum(abs(w) > 0)
                detectedXY <- data.frame(matrix(ncol = 2,
                    temp$detectedXY[1:(2*nd)]))
                names(detectedXY) <- c('x','y')
                attr(w, 'detectedXY') <- detectedXY
            }
            else
                attr(w, 'detectedXY') <- NULL

        }

        else
        if (detector(traps) %in% c('proximity', 'presence','unmarked')) {
            binomN <- 1

            temp <- .C("trappingproximity", PACKAGE = 'secr',
                as.double(df0),
                as.double(sigma),
                as.double(z),
                as.integer(noccasions),
                as.integer(k),
                as.integer(N),
                as.double(animals),
                as.double(unlist(traps)),
                as.double(distmat),
                as.double(usge),
                as.integer(detectfn),
                as.double(truncate^2),
                as.integer(binomN),
                n = integer(1),
                caught = integer(N),
                value = integer(N*noccasions*k),
                resultcode = integer(1)
            )
            if (temp$resultcode != 0)
                stop ("call to 'trappingproximity' failed")
            w <- array(dim=c(noccasions, k, temp$n), dimnames =
                list(1:noccasions,NULL, NULL))
            if (temp$n>0) {
                w[,,] <- temp$value[1:(temp$n*noccasions*k)]
            }
            w <- aperm(w, c(3,1,2))
        }
        else
        ## includes presence 2011-09-26
        if (detector(traps) %in% c('count')) {
            binomN <- detectpar$binomN
            temp <- .C("trappingcount", PACKAGE = 'secr',
                as.double(df0),
                as.double(sigma),
                as.double(z),
                as.integer(noccasions),
                as.integer(k),
                as.integer(N),
                as.double(animals),
                as.double(unlist(traps)),
                as.double(distmat),
                as.double(usge),
                as.integer(detectfn),
                as.double(truncate^2),
                as.integer(binomN),
                n = integer(1),
                caught = integer(N),
                value = integer(N*noccasions*k),
                resultcode = integer(1)
            )
            if (temp$resultcode != 0)
                stop ("call to 'trappingcount' failed")
            w <- array(dim=c(noccasions, k, temp$n), dimnames =
                list(1:noccasions,NULL, NULL))
            if (temp$n>0) {
                w[,,] <- temp$value[1:(temp$n*noccasions*k)]
            }
            w <- aperm(w, c(3,1,2))
            
        }
        else

        if (detector(traps) %in% c('signal','signalnoise')) {
            if (detectpar$binomN != 1)
                stop ("binomN != 1 not yet implemented for signal detectors")
            temp <- .C("trappingsignal", PACKAGE = 'secr',
                as.double(beta0),
                as.double(beta1),
                as.double(sdS),
                as.double(cutval),
                as.double(muN),       # used only if signalnoise detector type
                as.double(sdN),       # used only if signalnoise detector type
                as.double(sdM),
                as.integer(noccasions),
                as.integer(k),
                as.integer(N),
                as.double(animals),
                as.double(unlist(traps)),
                as.double(distmat),
                as.double(usge),
                as.integer(detectfn),
                n = integer(1),
                caught = integer(N),
                signal = double(N*noccasions*k),
                noise = double(N*noccasions*k),
                value = integer(N*noccasions*k),    # detected/not detected
                resultcode = integer(1)             # 0,1,2
            )
            if (temp$resultcode != 0)
                stop ("call to 'trappingsignal' failed with resultcode ", temp$resultcode)
            w <- array(dim=c(noccasions, k, temp$n), dimnames =
                 list(1:noccasions,NULL,NULL))
            if (temp$n>0)  {
                w[,,] <- temp$value[1:(temp$n * noccasions * k)]
            }
            w <- aperm(w, c(3,1,2))
        }
        else
        if (detector(traps) == 'times') {
            temp <- .C("trappingtimes", PACKAGE = 'secr',
                as.double(df0),
                as.double(sigma),
                as.double(z),
                as.integer(noccasions),
                as.integer(k),
                as.integer(N),
                as.double(animals),
                as.double(unlist(traps)),
                as.double(distmat),
                as.double(usge),
                as.integer(detectfn),
                as.double(truncate^2),
                n = integer(1),
                caught = integer(N),
                times = double(N*noccasions*k),
                value = integer(N*noccasions*k),
                resultcode = integer(1)
            )
            if (temp$resultcode != 0)
                stop ("call to 'trappingtimes' failed")
            w <- array(dim=c(noccasions, k, temp$n), dimnames =
                list(1:noccasions,NULL,NULL))
            if (temp$n>0)  {
                w[,,] <- temp$value[1:(temp$n * noccasions * k)]
            }
            w <- aperm(w, c(3,1,2))
            if (temp$n>0)  {
                attr(w, 'times') <- temp$times[1:sum(w)]
            }
        }
        else if (detector(traps) %in% c('polygon','transect','telemetry')) {
            if (detector(traps) == 'polygon') {
                nk <- length(levels(polyID(traps)))
                k <- table(polyID(traps))
                temp <- .C('trappingpolygon', PACKAGE = 'secr',
                           as.double(df0),
                           as.double(sigma),
                           as.double(z),
                           as.integer(noccasions),
                           as.integer(nk),
                           as.integer(k),
                           as.integer(N),
                           as.double(animals),
                           as.double(unlist(traps)),
                           as.double(usge),
                           as.integer(detectfn),
                           as.double(truncate^2),
                           as.integer(detectpar$binomN),
                           as.integer(maxperpoly),
                           n = integer(1),
                           caught = integer(N),
                           detectedXY = double (N*noccasions*nk*maxperpoly*2),
                           value = integer(N*noccasions*nk),
                           resultcode = integer(1)
                           )
                polynames <- levels(polyID(traps))
            }
            if (detector(traps) == 'transect') {
                nk <- length(levels(transectID(traps)))
                k <- table(transectID(traps))
                temp <- .C('trappingtransect', PACKAGE = 'secr',
                           as.double(df0),
                           as.double(sigma),
                           as.double(z),
                           as.integer(noccasions),
                           as.integer(nk),
                           as.integer(k),
                           as.integer(N),
                           as.double(animals),
                           as.double(unlist(traps)),
                           as.double(usge),
                           as.integer(detectfn),
                           as.double(truncate^2),
                           as.integer(detectpar$binomN),
                           as.integer(maxperpoly),
                           n = integer(1),
                           caught = integer(N),
                           detectedXY = double (N*noccasions*nk*maxperpoly*2),
                           value = integer(N*noccasions*nk),
                           resultcode = integer(1)
                           )
                polynames <- levels(polyID(traps))
            }
            if (detector(traps) == 'telemetry') {
                nk <- 1
                polynames <- '1'
                if (is.null(exactN)) exactN <- 0
                temp <- .C('trappingtelemetry', PACKAGE = 'secr',
                           as.double(df0),
                           as.double(sigma),
                           as.double(z),
                           as.integer(noccasions),
                           as.integer(N),
                           as.double(animals),
                           as.integer(detectfn),
                           as.double(truncate^2),
                           as.integer(detectpar$binomN),
                           as.integer(exactN),
                           as.integer(maxperpoly),
                           n = integer(1),
                           caught = integer(N),
                           detectedXY = double (N*noccasions*50*maxperpoly*2),
                           value = integer(N*noccasions),
                           resultcode = integer(1)
                           )
            }
            if (temp$resultcode != 0) {
                if (temp$resultcode == 2)
                    stop ("more than ", maxperpoly, "  detections per animal",
                          " per polygon per occasion")
                else
                    stop ("call to ",simfunctionname, " failed")
            }
            w <- array(dim=c(noccasions, nk, temp$n),
                dimnames = list(1:noccasions, polynames, NULL))

            if (temp$n > 0) {
                w[,,] <- temp$value[1:prod(dim(w))]
            }
            w <- aperm(w, c(3,1,2))

            if (temp$n > 0) {
                ## put XY coordinates in attribute
                nd <- sum(abs(w))
                detectedXY <- data.frame(matrix(ncol = 2,
                    temp$detectedXY[1:(2*nd)]))
                names(detectedXY) <- c('x','y')
                attr(w, 'detectedXY') <- detectedXY
            }
            else
                attr(w, 'detectedXY') <- NULL
        }
        else stop ('Unrecognised detector type')

        if (!is.null(covariates(popn))) {
            covariates(w) <- covariates(popn)[as.logical(temp$caught),, drop=F]
        }

        class(w)             <- 'capthist'    ## NOT data.frame
        traps(w)             <- traps
        attr(w, 'cutval')    <- cutval
        attr(w, 'seed')      <- RNGstate      ## save random seed
        attr(w, 'detectpar') <- detectpar
        session(w)           <- '1'           ## dummy session values for now

        ## optionally save population, whether simulated or input  2014-04-27
        if (savepopn)
            attr(w, 'popn') <- popn

        if (detector(traps) %in% c('signal','signalnoise')) {
            if (temp$n>0)  {
                signal(w) <- temp$signal[1:sum(w)]
                if (detector(traps) %in% c('signalnoise'))
                    noise(w) <- temp$noise[1:sum(w)]
            }
        }
        if (detector(traps) %in% 'telemetry')  ## 2013-11-21
            if (!is.na(chulltol))
             if (chulltol >= 0) w <- refreshMCP(w, chulltol)

        if (renumber && (temp$n>0)) rownames(w) <- 1:temp$n
        else {
            rown <- rownames(popn)[temp$caught > 0]
            caught <- temp$caught[temp$caught>0]
            rownames(w) <- rown[order(caught)]
        }

        w
    }
}
############################################################################################

## 2015-10-02 sim.resight substantially revised
## 2015-10-02 markocc (previously q) may be vector (0/1/2) length noccasions
##            0 = resighting, 1 = marking, 2 = marking and resighting
## pmark is probability an individual is marked, across the entire population, if all occasions
## are sighting occasions

sim.resight <- function (traps, popn = list(D = 5, buffer = 100, Ndist = 'poisson'), ...,
                         pID = 1, unmarked = TRUE, nonID = TRUE, unsighted = TRUE,
                         pmark = 0.5, Nmark = NULL, markingmask = NULL) {

    ############################################################################
    ## checks
    if (!(detector(traps) %in% c('proximity', 'count', 'polygon', 'transect')))
        stop ("only for binary or count proximity detectors or polygon or transect searches")

    if (is.null(markocc(traps))) {
        warning ("using default markocc 1 0 0 0 0")
        markocc <- c(1,0,0,0,0)
    }
    else
        markocc <- markocc(traps)
    allsighting <- !any(markocc>0)

    ## markocc(traps) <- markocc
    markocc(traps) <- NULL
    
    ## special case: Nmark animals distributed according to mask covariate
    distributedmarking <- FALSE
    if (!is.null(markingmask) & !is.null(Nmark)) {
        if (!is.null(covariates(markingmask)))
            distributedmarking <- 'marking' %in% names(covariates(markingmask))
    }
    
    ############################################################################
    ## make complete, identified capthist including all detections
    ## marking and sighting

    dots <- list(...)
    dots$traps <- traps

    if (!is.null(dots$seed)) {
        set.seed(dots$seed)
        dots$seed <- NULL
    }
    dots$noccasions <- length(markocc)  ## override noccasions
    dots$renumber <- FALSE  ## to match animalID

    if (!inherits(popn,'popn'))         ## generate popn if not provided
    {
        popn <- replacedefaults(list(D = 5, buffer = 100,
                                     Ndist = 'poisson'), popn)
        if (!(distributedmarking & allsighting))
            popn <- sim.popn (popn$D, core = traps, buffer = popn$buffer,
                              covariates = NULL, Ndist = popn$Ndist)
    }
    #---------------------------------------------------------------------------
    if (allsighting) {
        if (distributedmarking) {
            Nunmark <- if (popn$Ndist == 'fixed') {
                if (!is.null(popn$Nbuffer)) popn$Nbuffer - Nmark
                else popn$D * maskarea(markingmask) - Nmark
            }
            else { 
                rpois(1, popn$D * maskarea(markingmask)) - Nmark
            }
            if (Nunmark<0) {
                warning ("Nunmark < 0; setting to 0")
                Nunmark <- 0
            }
            covariates(markingmask)$marking <- covariates(markingmask)$marking / sum(covariates(markingmask)$marking)
            covariates(markingmask)$unmarking <- 1 - covariates(markingmask)$marking
            popnM <- sim.popn(D = 'marking', core = markingmask, Ndist = 'fixed', Nbuffer = Nmark, model2D='IHP')
            popnU <- sim.popn(D = 'unmarking', core = markingmask, Ndist = 'fixed', Nbuffer = Nunmark, model2D='IHP')
        }
        else {
            gotten <- rep(TRUE, nrow(popn))
            
            ## circumscribe marked animals
            if (!is.null(markingmask)) {
                gotten <- pointsInPolygon(popn, markingmask)
                if (!is.null(covariates(markingmask)))
                    if ('marking' %in% names(covariates(markingmask)))
                        warning("uniform distribution over marking mask;",
                                " covariate 'marking' ignored")
            }
            
            if (!is.null(Nmark)) {
                if (Nmark > sum(gotten)) {
                    warning("Nmark exceeds population size, marking all")
                }
                else {
                    gottenN <- sample.int(sum(gotten), size = Nmark, replace = FALSE)
                    gotten[gotten] <- (1:sum(gotten)) %in% gottenN
                }
            }
            else {
                gotten <- gotten & (runif(length(gotten)) < pmark)
            }
            popnM <- subset(popn, gotten)     ## inside A0 and marked
            popnU <- subset(popn, !gotten)    ## outside A0 or not marked
        }
        dots$popn <- popnM
    }
    else
        dots$popn <- popn
    #---------------------------------------------------------------------------

    capthist <- do.call('sim.capthist', dots)

    ############################################################################

    ## transform simulated capthist object into resighting data
    S <- length(markocc)
    K <- ndetector(traps)

    session <- rep(1,sum(abs(capthist)))
    df <- data.frame(session, animalID=animalID(capthist), occasion=occasion(capthist),
                     trapID = trap(capthist), stringsAsFactors = FALSE)
    if (detector(traps) %in% c('polygon', 'transect')) {
        df <- cbind(df[,-4], xy(capthist))
        df$trapID <- xyinpoly(df[,4:5], traps(capthist))
    }

    if (nrow(df) == 0) {
        dfID <- data.frame (session=1, animalID="NONE", occasion = 1, trapID=1)
        dfnonID <- dfID[-1,]
    }
    else {
        # data frame of detections on marking occasions
        dfmarking <- df[markocc[df$occasion]>0,,drop=FALSE]
        # when was each animal first marked?
        df$markingoccasion <- dfmarking$occasion[match(df$animalID, dfmarking$animalID)]
        df$markingoccasion[is.na(df$markingoccasion)] <- S+1    ## never marked
        if (allsighting) df$markingoccasion <- rep(0, nrow(df)) ## all marked
        # separate data frames for detections of marked and unmarked animals
        dfmarked <- df[df$occasion >= df$markingoccasion,, drop = FALSE]
        dfunmarked <- df[df$occasion < df$markingoccasion,, drop = FALSE]
        # which of the marked detections were identified?
        dfmarked$ID <- (dfmarked$occasion == dfmarked$markingoccasion) |
            (runif(nrow(dfmarked)) < pID)
        # separate data frames for detections of ID and nonID marked animals
        dfnonID <- dfmarked[!dfmarked$ID, , drop = FALSE]
        dfID    <- dfmarked[dfmarked$ID, , drop = FALSE]
    }

    if (detector(traps) %in% c('polygon', 'transect')) {
        CH <- make.capthist(dfID[,1:5, drop=FALSE], traps, fmt = 'XY', noccasions = S)
    }
    else {
        CH <- make.capthist(dfID[,1:4, drop=FALSE], traps, fmt = 'trapID', noccasions = S)
    }

    ## marked but never sighted
    if (allsighting & unsighted) {
        nzero <- nrow(popnM) - nrow(CH)
        CH <- addzeroCH(CH, nzero)
    }
    ############################################################################

    if (detector(traps(capthist)) %in% c('polygon', 'transect'))
        trapnames <- levels(polyID(traps(capthist)))
    else
        trapnames <- row.names(traps(capthist))

    ##------------------------------------------------------------------------------

    unresolved <- markocc == -1
    markocc(traps(CH)) <- markocc
    if (unmarked) {

        countfn <- function(x) {
            x <- x * col(x)   ## x 0/1
            tabulate(x[x>0], nbins = K)
        }

        if (allsighting) {
            ## apply same sampling to unmarked fraction of population
            dots$popn <- popnU
            capthistU <- do.call('sim.capthist', dots)
            Tu <- as.matrix(table(factor(trap(capthistU), levels = trapnames),
                                  factor(occasion(capthistU), levels = 1:S)))
        }
        else {
            Tu <- as.matrix(table(factor(dfunmarked$trapID, levels = trapnames),
                                  factor(dfunmarked$occasion, levels = 1:S)))
        }
        Tu[,markocc>0] <- 0  ## just in case...
        
        ## 2015-12-15
        if (sum(unresolved)>0) {
            markedsightings <- t(apply(CH[,unresolved,], 2:3, sum, drop = FALSE))
            CH[,unresolved,] <- 0
            Tu[,unresolved] <- Tu[,unresolved] + markedsightings
        }    
        
        Tu(CH) <- Tu
    }
    ##------------------------------------------------------------------------------

    if (nonID) {
        Tm <- as.matrix(table(factor(dfnonID$trapID, levels = trapnames),
                                          factor(dfnonID$occasion, levels = 1:S)))
        ## 2015-12-15
        if (sum(unresolved)>0) {
            ## add to unmarked sightings
            Tu(CH)[,unresolved] <- Tu[,unresolved] + Tm[,unresolved]
            Tm[,unresolved] <- 0
        }
        Tm(CH) <- Tm 
    }
    ##------------------------------------------------------------------------------

    ## if savepopn = TRUE...
    if (!is.null(attr(capthist, 'popn'))) {
        if (allsighting) {
            popn <- rbind(popnM, popnU)
            covariates(popn) <- data.frame(marked = rep(c(TRUE, FALSE), c(nrow(popnM), nrow(popnU))))
            attr(CH, 'popn') <- popn
        }
        else {
            ## 2015-12-18; this could be moved to sim.capthist
            covariates(popn) <- data.frame(marked = rownames(popn) %in% dfmarked$animalID)   
            attr(CH, 'popn') <- popn   ## NULL unless 
        }
    }
    CH
}
############################################################################################

