###############################################################################
## package 'secr'
## sim.popn.R
## simulate spatially distributed population
## transferred from methods.R 2010-06-08
## last changed
## 2010-06-14  added Ndist = 'specified', using D to carry N
## 2011-03-27  conditional attachment of rownames (correct bug when N = 0)
## 2011-04-06  coastal beta option
## 2011-10-20  poly argument
## 2012-04-10  IHP D may be covariate; more checks
## 2012-04-10  MRC
## 2013-11-23 IHP for Ndist fixed
## 2014-04-18 session-specific density
## 2014-09-03 model2D = "linear" option, and some recoding of IHP
## 2014-12-29 buffertypes concave, convex
## 2015-01-13 debugged non-rectangular buffers; see sim.popn.test.R in testing
## 2015-02-18 multisession sim.popn updated for Nbuffer
###############################################################################

toroidal.wrap <- function (pop) {
    bb <- attr(pop, 'boundingbox')
    xmin <- min(bb$x)
    xmax <- max(bb$x)
    ymin <- min(bb$y)
    ymax <- max(bb$y)
    xrange <- xmax-xmin
    yrange <- ymax-ymin
    remainder <- function (x,y) x - (x %/% y * y)
    pop$x <- ifelse (pop$x>xmax, xmin + remainder (pop$x-xmax, xrange), pop$x)
    pop$x <- ifelse (pop$x<xmin, xmax - remainder (xmin - pop$x, xrange), pop$x)
    pop$y <- ifelse (pop$y>ymax, ymin + remainder (pop$y-ymax, yrange), pop$y)
    pop$y <- ifelse (pop$y<ymin, ymax - remainder (ymin - pop$y, yrange), pop$y)
    pop
}
tile <- function (popn, method = "reflect") {
    bbox <- attr(popn, 'boundingbox')
    if (method== "reflect") {
        p2 <- rbind.popn(popn, flip(popn,lr=min(bbox$x)), flip(popn,lr=max(bbox$x)))
        rbind.popn(p2, flip(p2,tb=min(bbox$y)), flip(p2,tb=max(bbox$y)))
    }
    else if (method == "copy") {
        ht <- max(bbox$y) - min(bbox$y)
        wd <- max(bbox$x) - min(bbox$x)
        p2 <- rbind.popn(popn, shift(popn,c(-wd,0)), shift(popn, c(wd,0)))
        rbind.popn(p2, shift(p2,c(0,-ht)), shift(p2, c(0,ht)))
    }
    else
        stop ("unrecognised method")
}

sim.popn <- function (D, core, buffer = 100, model2D = c("poisson",
    "cluster", "IHP", "coastal", "hills", "linear"), buffertype =
    c("rect", "concave", "convex"), poly = NULL,
    covariates = list(sex = c(M = 0.5,F = 0.5)), number.from = 1, Ndist
    = c('poisson','fixed','specified'), nsessions = 1, details = NULL,
    seed = NULL, keep.mask = model2D %in% c('IHP','linear'), Nbuffer = NULL,
    ...)  {

    model2D <- match.arg(model2D)
    Ndist <- match.arg(Ndist)
    buffertype <- match.arg(buffertype)
    if (buffertype %in% c('convex','concave') & (model2D != 'poisson'))
        stop ("buffertype incompatible with model2D")

    if (nsessions > 1) {
        discrete <- function(x) {
            fr <- x-trunc(x)
            sample (c(trunc(x), trunc(x)+1), size=1, prob=c(1-fr, fr))
        }
        session.popn <- function (s, D, Nbuffer) {
            ## independent population
            if (s > 1) seed <- NULL   ## 2015-02-18
            if (!is.null(Nbuffer))
                if (is.na(Nbuffer)) Nbuffer <- NULL
            sim.popn (D[1], core, buffer, model2D, buffertype, poly,
                covariates, number.from, Ndist, nsessions = 1, details, seed,
                keep.mask, Nbuffer[1])
        }
        turnover <- function (oldpopn, t) {
            ## project existing population
            ## assume normal movement kernel
            ## assume lambda lacks process variance
            ## ideally lambda lognormal
            ## need 'wrap' option for toroidal wrapping of 'rect' locations
            newstart <- max(as.numeric(rownames(oldpopn))) + 1
            if (turnoverpar$survmodel=='binomial') {
                survive <- sample (c(FALSE, TRUE), nrow(oldpopn), replace = TRUE,
                    c(1-turnoverpar$phi[t],turnoverpar$phi[t]))
                nsurv <- sum(survive)
            }
            else {   ## assume 'discrete'
                nsurv <- discrete (turnoverpar$phi[t] * nrow(oldpopn))
                survive <- sample (nrow(oldpopn), replace = FALSE, size = nsurv)
                survive <- sort(survive)   ## numeric indices
            }
            newpopn <- subset.popn(oldpopn, subset=survive)
            if (turnoverpar$sigma.m[t] > 0) {
                newpopn[,] <- newpopn[,] + rnorm (2*nsurv, mean = 0,
                                                  sd = turnoverpar$sigma.m[t])
                if (turnoverpar$wrap)
                    newpopn <- toroidal.wrap(newpopn)
            }
            gam <- turnoverpar$lambda[t] - turnoverpar$phi[t]
            if (gam<0)
                stop ("invalid gamma in turnover")
            nrecruit <- switch (turnoverpar$recrmodel,
                constantN = nrow(oldpopn) - nsurv,
                discrete = discrete(gam * nrow(oldpopn)),
                binomial = rbinom(1, nrow(oldpopn), gam),
                poisson = rpois (1, gam * nrow(oldpopn)))
                ## under Pradel model members of superpopulation have binomial
                ## probability of entry at this point?
                ## cf Schwarz & Arnason betas
            if (nrecruit>0) {
#               recruits <- sim.popn(D = nrecruit, core = core, buffer = buffer,
#                                    model2D = model2D, buffertype = buffertype, poly = poly,
#                                    covariates = covariates, number.from = newstart,
#                                    Ndist = 'specified', nsessions = 1, details = details)
              ## 2015-04-06 using Nbuffer
              recruits <- sim.popn(D = D, core = core, buffer = buffer,
                                   model2D = model2D, buffertype = buffertype, poly = poly,
                                   covariates = covariates, number.from = newstart,
                                   Ndist = 'specified', Nbuffer = nrecruit,
                                   nsessions = 1, details = details)
              ## danger: resets random seed
                newpopn <- rbind.popn(newpopn, recruits, renumber = FALSE)
            }
            class(newpopn) <- class(MSpopn[[1]])
            attr(newpopn, 'mask') <- attr(MSpopn[[1]], 'mask')
            ## minor attributes are neglected... beware in future
            attr(newpopn, 'losses') <- nrow(oldpopn)-nsurv
            attr(newpopn, 'recruits') <- nrecruit
            newpopn
        }
        turnoverpar <- list(lambda = NULL, phi = 0.7, sigma.m = 0, wrap = TRUE,
                            survmodel = 'binomial', recrmodel = 'poisson')
        expands <- function (param, s) {
            if (!is.null(param)) {
                param <- rep(param,s)[1:s]
            }
            param
        }
        turnoverpar <- replace (turnoverpar, names(details), details)
        
        if (is.null(details$lambda)) {
            ## independent
            ## MSpopn <- lapply (1:nsessions, session.popn)
            ## 2014-04-18, 2015-02-18
            if (missing(D))
                D <- rep(NA, nsessions)
            else if (length(D) == 1)
                D <- rep(D, nsessions)
            else
                if (length(D) != nsessions) stop ("length(D) should equal nsessions")
            if (is.null(Nbuffer))
                Nbuffer <- rep(NA, nsessions)
            else if (length(Nbuffer) == 1)
                Nbuffer <- rep(Nbuffer, nsessions)
            else
                if (length(Nbuffer) != nsessions) stop ("length(Nbuffer) should equal nsessions")

            MSpopn <- mapply (session.popn, 1:nsessions, D, Nbuffer, SIMPLIFY = FALSE)
        }
        else {
            ## projected
            MSpopn <- vector(nsessions, mode = 'list')
            MSpopn[[1]] <- session.popn(1, D, Nbuffer)
            turnoverpar$lambda  <- expands(turnoverpar$lambda, nsessions)
            turnoverpar$phi     <- expands(turnoverpar$phi, nsessions)
            turnoverpar$sigma.m <- expands(turnoverpar$sigma.m, nsessions)
            for (i in 2:nsessions) {
                MSpopn[[i]] <- turnover(MSpopn[[i-1]], i-1)
            }
        }
        if (model2D == 'linear')
            class(MSpopn) <- c('list', 'linearpopn', 'popn')
        else
            class(MSpopn) <- c('list','popn')
        names(MSpopn) <- 1:nsessions
        MSpopn
    }
    else {
        ##########################
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
        ##########################

        getnm <- function (scaleattribute = 'area', scale = 1, D) {
            
            ## 2014-09-03 for "IHP" and "linear"
            if ((length(D) == 1) & (is.character(D)))
                D <- covariates(core)[,D]
            if ((length(D) == 1) & (is.numeric(D)))
                D <- rep(D, nrow(core))
            if (any(is.na(D))) {
                D[is.na(D)] <- 0
                warning ("NA values of D set to zero")
            }
            if (any(D<0)) {
                D <- pmax(0,D)
                warning ("negative D set to zero")
            }
            ## D vector, 1 per cell
            D <- D * attr(core, scaleattribute) * scale
            ## 2015-04-06 using Nbuffer
            if (!is.null(Nbuffer)) {  ## includes Ndist == 'specified'
              N <- round(Nbuffer)
            }
            else {
              N <- sum(D)
            }
            if (Ndist == 'poisson') {
                N <- rpois(1,N)
            }
            rmultinom (1, N, D)
        }
        ##########################

        if (model2D %in% c('IHP')) {
            if (!inherits(core, 'mask'))
                stop ("for model2D = IHP, 'core' should be a habitat mask")
            nm <- getnm('area', 1, D)
            jitter <- matrix ((runif(2*sum(nm))-0.5) * attr(core,'spacing'), ncol = 2)
            animals <- core[rep(1:nrow(core), nm),] + jitter
            animals <- as.data.frame(animals)
            xl <- range(animals[,1])
            yl <- range(animals[,2])
        }
        else  if (model2D == 'linear') {
            if (!inherits(core, 'linearmask'))
                stop ("for model2D = linear, 'core' should be a linear mask")
            nm <- getnm('spacing', 0.001, D)
            animals <- core[rep(1:nrow(core), nm),] * 1  ## * 1 to shed attributes...
            animals <- as.data.frame(animals)
            xl <- range(animals[,1])
            yl <- range(animals[,2])

        }
        else {
            ## 2014-12-29, 2015-01-11
            # population in arena +/- buffer from traps
            buff <- c(-buffer,+buffer)
            xl <- range(core$x) + buff
            yl <- range(core$y) + buff

            area <- diff(xl) * diff(yl) * 0.0001  # ha not sq metres

            bufferpoly <- switch(buffertype,
                                 rect = NA,
                                 convex = buffer.contour(core, buffer = buffer,
                                                         convex = TRUE, plt = FALSE)[[1]],
                                 concave = buffer.contour(core, buffer = buffer,
                                                          convex = FALSE, plt = FALSE))

            bufferarea <- switch (buffertype,
                                  rect = area,
                                  convex = polyarea (bufferpoly),
                                  concave = sum(sapply(bufferpoly, polyarea)))

            if ((buffertype == 'concave') & (is.null(Nbuffer)))
                warning("automatic Nbuffer unreliable with concave buffer")

            ## If target number (Nbuffer) not specified, get from density x area
            ## N is the number of centres to simulate in the unbuffered (rectangular) area
            if (is.null(Nbuffer)) {
                N  <- switch (Ndist,
                              poisson = rpois(1, lambda = D[1] * area),
                              fixed = discreteN (1, D[1] * area),
                              specified = round(D[1]))
                Nbuffer  <- switch (Ndist,
                                    poisson = rpois(1, lambda = D[1] * bufferarea),
                                    fixed = discreteN (1, D[1] * bufferarea),
                                    specified = round(D[1]))
            }
            else {
                Nbuffer  <- switch (Ndist,
                                    poisson = rpois(1, lambda = Nbuffer),
                                    fixed = discreteN (1, Nbuffer),
                                    specified = round(Nbuffer))
                N <- Nbuffer
            }

#             if (buffertype %in% c('convex','concave')) {
#                 N <- N * area / bufferarea
#             }

            if (model2D == 'poisson') {
                animals <- data.frame (
                    x = runif(N)*diff(xl)+xl[1],
                    y = runif(N)*diff(yl)+yl[1])

                ## 2014-12-29, 2015-01-11, 2015-01-13 allow buffering / poly
                if (buffertype %in% c('convex','concave')) {

#                     if (Ndist == 'fixed') {
                        maxtries <- 10
                        tries <- 1
                        animals <- switch(buffertype,
                                      convex = animals[pointsInPolygon(animals, bufferpoly),],
                                      concave = animals[distancetotrap(animals, core)<= buffer,])
                        ## repeat if not enough
                        while ((nrow(animals) < Nbuffer) & (tries < maxtries)) {
                            animals <- rbind(animals, data.frame (
                                                                  x = runif(N)*diff(xl)+xl[1],
                                                                  y = runif(N)*diff(yl)+yl[1]))
                            animals <- switch(buffertype,
                                        convex = animals[pointsInPolygon(animals, bufferpoly),],
                                        concave = animals[distancetotrap(animals, core)<= buffer,])
                            tries <- tries + 1
                        }
                        if (tries >= maxtries)
                            warning("exceeded maxtries in sim.popn")
                        animals <- animals[1:Nbuffer,]
#                     }
#                     else {
#                         animals <- switch(buffertype,
#                                       convex = animals[pointsInPolygon(animals, bufferpoly),],
#                                       concave = animals[distancetotrap(animals, core)<= buffer,])
#                     }
                }
            }
            else if (model2D == 'coastal') {
                if (is.null(details$Beta))
                    details$Beta <- c(1,1.5,5,1)
                a1 <- details$Beta[1]
                b1 <- details$Beta[2]
                a2 <- details$Beta[3]
                b2 <- details$Beta[4]
                animals <- data.frame (x = rbeta(N,a1,b1)*diff(xl)+xl[1],
                                       y = rbeta(N,a2,b2)*diff(yl)+yl[1])
            }
            else if (model2D=='hills') {
                hills <- details$hills
                if (is.null(hills)) hills <- c(1,1)
                hills <- c(hills, rep(0,4-length(hills)))
                nhillx <- abs(hills[1])
                nhilly <- abs(hills[2])
                offset <- any(hills[1:2]<0)
                dx <- hills[3]
                dy <- hills[4]
                dx <- ifelse (dx<0, runif(1), dx)
                dy <- ifelse (dy<0, runif(1), dy)
                xhill <- sample(0:(nhillx-1), N, replace=TRUE)
                yhill <- sample(0:(nhilly-1), N, replace=TRUE)
                x <- asin(runif(N)*2 - 1)/pi + 0.5 + xhill
                y <- asin(runif(N)*2 - 1)/pi + 0.5 + yhill
                if (offset) x <- x + (yhill %% 2) * 0.5
                x <- x/nhillx + dx
                y <- y/nhilly + dy
                x <- ifelse (x>1, x-1, x)
                y <- ifelse (y>1, y-1, y)
                animals <- data.frame (x = x * diff(xl)+xl[1],
                                       y = y * diff(yl)+yl[1])
            }
            else if (model2D=='cluster') {
                ## Neyman-Scott distribution with wrapping
                xrange <- diff(xl)
                yrange <- diff(yl)
                if (details$mu<=0) {
                    nparent <- N   ## not clustered
                    offspr <- sweep(matrix(runif(2*nparent), ncol = 2), 2, c(xrange,yrange), '*')
                }
                else {
                     nparent <- switch (Ndist,
                         poisson = rpois(1, lambda=D[1] * area/details$mu),
                         fixed = discreteN (1, D[1] * area / details$mu),
                         specified = discreteN (1, D[1] / details$mu))  ## here arg D is N
                     N <- nparent * details$mu
                     if (nparent==0)
                         warning ("zero clusters")
                     parent <-  sweep(matrix(runif(2*nparent), ncol = 2), 2, c(xrange,yrange), '*')

                     offspr <- matrix(rnorm(2*N), ncol = 2) * details$hsigma
                     parentn <- rep(1:nparent, details$mu)
                     offspr <- offspr + parent[parentn,]
                     while (any ((offspr[,1]<0) | (offspr[,1]>xrange) | (offspr[,2]<0) |
                                 (offspr[,2]>yrange))) {
                       offspr[,1] <- ifelse (offspr[,1]<0, offspr[,1]+xrange, offspr[,1])
                       offspr[,1] <- ifelse (offspr[,1]>xrange, offspr[,1]-xrange, offspr[,1])
                       offspr[,2] <- ifelse (offspr[,2]<0, offspr[,2]+yrange, offspr[,2])
                       offspr[,2] <- ifelse (offspr[,2]>yrange, offspr[,2]-yrange, offspr[,2])
                     }
                }
                animals <- as.data.frame(sweep(offspr,2,c(xl[1],yl[1]),'+'))
            }
            else stop ("unrecognised 2-D distribution")
        }
        names(animals) <- c('x','y')
        attr(animals,'covariates') <- NULL
        if (!is.null(covariates)) {
            tempcov <- list()
            for (i in 1:length(covariates)) {
               covi <- sample (names(covariates[[i]]), replace = T, size= nrow(animals),
                               prob=covariates[[i]])
               temptxt <- paste ('tempcov$', names(covariates[i]), '<- covi',
                               sep = '')
               eval(parse(text=temptxt))
            }
            attr(animals,'covariates') <- as.data.frame(tempcov)
        }
        if (nrow(animals) > 0)   ## condition added 2011-03-27
            row.names (animals) <- number.from : (nrow(animals)+number.from-1)

        if (keep.mask) {
            if (model2D %in% c('IHP','linear'))
                attr(animals, 'mask') <- core
        }
        if (model2D == 'linear') {
            class(animals) <- c('linearpopn', 'popn', 'data.frame')
        }
        else {
            class(animals) <- c('popn', 'data.frame')
            ##-------------------------
            ## restrict to a polygon
            ## added 2011-10-20
            ## NOTE 2014-12-29 this breaks Ndist = 'fixed'
            if (!is.null(poly)) {
                animals <- subset(animals, poly = poly, ...)
            }
            ##-------------------------
        }
        attr(animals, 'seed') <- RNGstate   ## save random seed
        attr(animals, 'Ndist') <- Ndist
        attr(animals, 'model2D') <- model2D
        attr(animals, 'buffertype') <- buffertype
        attr(animals, 'boundingbox') <- expand.grid (x=xl,y=yl)[c(1,3,4,2),]
        animals
    }
}

