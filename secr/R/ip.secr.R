###############################################################################
## package 'secr'
## ip.secr.R
## last changed 2009 09 12, 2009 09 14
## 2010 07 01 alpha detection function
## 2011 06 15 tidy up
## 2011 12 20 maxtries argument
## 2012-11-02 ncores
## 2012-11-02 proctime[3]
## 2012-12-30 usage OK
## 2015-05-17 tweaked failure condition
## 2015-11-24 replaced argument boxsize with boxsize1, boxsize2
###############################################################################

ip.secr <- function (capthist, predictorfn = pfn, predictortype = 'null',
    detectfn = 0, mask = NULL, start = NULL, boxsize = 0.2, boxsize2 = boxsize, centre = 3,
    min.nsim = 10, max.nsim = 2000, CVmax = 0.002, var.nsim = 1000,
    maxbox = 5, maxtries = 2, ncores = 1, seed = NULL, trace = TRUE, ...) {

    ## ... passed to sim.popn e.g. buffer = 100, Ndist = 'fixed'
    ## boxsize may be vector of length np
    ## pfn defined below
    if (maxbox<2) stop("ip.secr maxbox >= 2")
    if (is.list(capthist))
        stop ("'ip.secr' not implemented for multi-session 'capthist'")
    ptm  <- proc.time()
    cl <- match.call(expand.dots = TRUE)
    cl <- paste(names(cl)[-1],cl[-1], sep=' = ', collapse=', ' )
    cl <- paste('ip.secr(', cl, ')')

    ## added 2010-07-01
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)

    if (!(detectfn %in% 0:11))
        stop ("unrecognised detectfn")
    if (!(detectfn %in% c(0,2,4)))
        stop (detectionfunctionname(detectfn),
            " detection function not implemented in ip.secr")
    pnames <- c('D', parnames(detectfn))
    np <- length(pnames)
    traps <- traps(capthist)
    noccasions <- ncol(capthist)
    if (length(boxsize)==1) boxsize1 <- rep(boxsize, np)
    else if (length(boxsize1) != np)
        stop ("invalid boxsize vector")
    else boxsize1 <- boxsize
    if (length(boxsize2)==1) boxsize2 <- rep(boxsize2, np)
    else if (length(boxsize2) != np)
        stop ("invalid boxsize vector")

    if (is.null(mask))
        core <- expand.grid(x = range(traps$x), y = range(traps$y))
    else
        core <- NULL

    ## added 2012-11-02
    if (ncores > 1) {
        clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
        clusterSetRNGStream(clust, seed)
        clusterEvalQ(clust, requireNamespace('secr'))
        clusterExport(clust, c("capthist", "predictorfn", "predictortype",
                               "maxtries", "negloglikM0", "negloglikMb",
                               "jack.est", "M0", "Mb", "Mh", "traps", "pnames",
                               "noccasions", "mask", "core", "odds","invodds",
                               "RPSV", "pfn"),
                      environment())
    }
    else {
        clust <- NULL
        set.seed(seed)
    }

    # to simulate one realization
    simfn <- function (parval, ...) {
        D <- parval[1]
        detectpar <- as.list(parval[-1])
        names(detectpar) <- pnames[-1]
        detectpar[['g0']] <- invodds(detectpar[['g0']])
        attempts <- 0
        allOK <- FALSE
        repeat {
            if (is.null(mask))
                popn <- sim.popn (D = D, core = core, model2D='poisson', ...)
            else
                popn <- sim.popn (D = D, core = mask, model2D='IHP', ...)
            simcapthist <- sim.capthist(traps, popn, detectfn, detectpar,
                                        noccasions)
            if (nrow(simcapthist)==0) {
                warning ("ip.secr: no captures in simulation", call. = FALSE)
            }
            else
                if ((sum(abs(simcapthist)>0) - nrow(simcapthist)) < 1)
                warning ("ip.secr: no re-captures in simulation", call. = FALSE)
            predicted <- try (predictorfn (simcapthist, predictortype), silent = TRUE)
            if (inherits(predicted, 'try-error')) {
                predicted <- NA
            }
            attempts <- attempts+1

            ## exit loop if exceeded maxtries or all OK
            allOK <- !any(is.na(predicted)) & all(is.finite(predicted))
            if ((attempts >= maxtries) | allOK)
                break
        }
        if (!allOK) {
            ## 2015-05-17 replaced with allOK if (attempts >= maxtries) {
            stop ("ip.secr: no successful simulation after ", maxtries,
                  " attempts", call. = FALSE)
        }
        if (attempts > 1)
            warning ("ip.secr: simulation repeated", call. = FALSE)
        predicted
    }

    ## to test if current solution is within box
    within <- function (i) (par[i] >= vertices[[i]][1]) & (par[i] <=
                               vertices[[i]][2])

    ## target values of predictor
    y <- predictorfn(capthist, predictortype)
    if (length(y) != np)
        stop ("need one predictor for each parameter ",
             paste(pnames, collapse=" "))

    if (is.null(start)) {
        if (trace) {
            cat('\nFinding starting values ...\n')
            flush.console()
        }
        if (is.null(mask)) {
            automask <- make.mask(traps, buffer=3*RPSV(capthist, CC = TRUE))
            start <- unlist(autoini(capthist, automask))
        }
        else
            start <- unlist(autoini(capthist, mask))
        ## ad hoc bias adjustment
        if (detector(traps)=='single')
            start[2] <- invodds(odds(start[2]) * 1.4)
        if (detectfn %in% 1) start <- c(start,5)  ## z
    }
    par <- start
    names(par) <- pnames
    par['g0'] <- odds(par['g0'])

    ####################################################################
    ## keep trying until within box or exceeds maxbox

    for (m in 1:maxbox) {
        if (trace) {
            cat('\nFitting box', m, '...   (g0 on odds scale) \n')
        }
        names(par) <- pnames
        boxsize <- if (m == 1) boxsize1 else boxsize2
        vertices <- sweep (1 + outer(c(-1,1), boxsize), MARGIN = 2,
                           FUN = '*', STATS = par)
        vertices <- data.frame(vertices)
        names(vertices) <- pnames
        rownames(vertices) <- c('min','max')
        if (trace) {
            print(vertices)
            cat('\n')
            flush.console()
        }
        design <- as.matrix(expand.grid (as.list(vertices)))
        centrepoints <- matrix(par, nrow = centre, ncol = np, byrow=T)
        design <- rbind(design, centrepoints)
        basedesign <- design[rep(1:nrow(design), min.nsim),]
        sim <- NULL
        design <- NULL
        ## baseindices <- rep(c(1:2^np, rep(0,centre)), min.nsim)
        ## indices <- numeric(0)

        repeat {

            if (ncores > 1) {
                list(...) # evaluate any promises cf boot
                newsim <- parRapply(clust, basedesign, simfn, ...)   ## ...added 2012-11-08
                newsim <- t(matrix(newsim, ncol = nrow(basedesign)))
            }
            else {
                newsim <- t(apply(basedesign,1,simfn, ...))   ## ...added 2012-11-080
            }

            OK <- (newsim[,3]>0) & (!is.na(newsim[,3]))  ## require valid RPSV
            sim <- rbind(sim, newsim[OK,])
            ## indices <- c(indices, baseindices[OK])
            design <- rbind(design,basedesign[OK,])
            if ((nrow(design) > max.nsim)) break
            sim.lm <- lm ( sim ~ design )
            CV <- sapply(summary(sim.lm), function(x) x$sigma) / y /
                sqrt(nrow(sim))
            ## following is almost identical in effect; requires 'indices'
            ## CV <- apply (sim,2, function (x) tapply(x, indices, function(y)
            ##       sd(y)/mean(y)/sqrt(nrow(sim))))
            if (all(CV <= CVmax)) break
        }
        if (nrow(design)/(2^np+centre) > max.nsim) {
            warning ("exceeded maximum allowable replicates ",
                     "without achieving 'CVmax'")
            return (NA)
        }
        else {
            B <- coef(sim.lm)[-1,]
            B <- solve(t(B))  ## invert
            lambda <- coef(sim.lm)[1,]   ## intercepts
            par <- as.numeric(B %*% matrix((y - lambda), ncol = 1))
            ## 2015-11-24, 2015-12-02 if (all(sapply(1:np, within))) break
            ## only break on second or later box if differ boxsize
            if (all(sapply(1:np, within)) & (all(boxsize == boxsize2) | (m>1))) break
        }
    }

    if (!all(sapply(1:np, within)))
        warning ("solution not found after ", maxbox, " attempts")
    ####################################################################

    names(par) <- pnames

    if (var.nsim>1) {
        if (trace) {
            cat('Simulating for variance ...\n')
            flush.console()
            cat('\n')
        }
        vardesign <- matrix(par, nrow = var.nsim, ncol = np, byrow = T)
        colnames(vardesign) <- pnames

        if (ncores > 1) {
            list(...) # evaluate any promises cf boot
            newsim <- parRapply(clust, vardesign, simfn)
            newsim <- t(matrix(newsim, ncol = nrow(vardesign)))
        }
        else {
            newsim <- t(apply(vardesign,1,simfn))
        }

        V <- var(newsim)  ## additional simulations for var-covar matrix
        vcov <- B %*% V %*% t(B)

        ## compare estimates to parametric bootstrap
        n <- apply(newsim, 2, function(x) sum(!is.na(x)))
        ymean <- apply(newsim, 2, mean, na.rm=T)
        yse <- apply(newsim, 2, function(x) sd(x, na.rm=T) / sum(!is.na(x)))

        bootstrap <- data.frame (target = y, nsim = n, simulated = ymean,
            SE.simulated = yse)

        ## biasest not reported, yet
        yest <- as.numeric(B %*% matrix((ymean - lambda), ncol = 1))
        biasest <- data.frame (estimate = 100 * (par - yest) / yest,
            SE = 100 * (par - yest) / yest)

        ## adjust var-covar matrix for g0 using delta method
        tx <- diag(np)
        dimnames (tx) <- list(pnames, pnames)
        tx['g0','g0'] <- -(1/(1 + par['g0'])^2)  ## gradient invodds(y) wrt y
        vcov <- tx %*% vcov %*% t(tx)
    }
    else {
        vcov <- matrix(nrow = np, ncol = np)
        bootstrap <- NA
    }

    dimnames(vcov) <- list(pnames, pnames)
    par['g0'] <- invodds(par['g0'])

    if (ncores > 1) {
        stopCluster(clust)
    }

    list(call = cl,
        IP = data.frame(estimate=unlist(par), SE.estimate=diag(vcov)^0.5),
        vcov = vcov,
        ip.nsim = nrow(sim),
        variance.bootstrap = bootstrap,
        proctime = as.numeric((proc.time() - ptm)[3])
    )
}
##################################################


##################################################
## population size likelihoods
negloglikM0   <- function (theta, n) {
  NN   <- exp(theta)
  ndot <- n[1]  # total captures
  Mt1  <- n[2]  # total animals
  K    <- n[3]
  LL   <- lgamma (NN+1) - lgamma(NN-Mt1+1) + ndot * log(ndot) +
          (K*NN - ndot) * log(K*NN - ndot) - K*NN*log(K*NN)
  if (is.finite(LL)) -LL
  else 1e10
}

negloglikMb   <- function (theta, u) {
  p <- invlogit(theta)
  K <- length(u)
  n  <- sum(u)
  -(n * log(p)+ sum(u *(0:(K-1))) * log (1-p) - n * log(1 - (1-p)^K))
}
##################################################

jack.est <- function (inp, deads = 0, full = F)
{

# Calculate Burnham & Overton's jackknife estimate for closed populations
#
# inp may be a vector of capture frequencies, or
#            a list comprising such a vector as its first element and
#                              the scalar number of 'deads' as its second
#
# 'deads' are assumed to be additional to the tabulated capture frequencies
# They are added to the calculated population size (cf Otis et al. 1978)
#
# MODIFIED
# 30/3/95 Fix bug when mt=5 (last jacknife selected) - length(test) s/b 5
# 30/3/95 Optional full output
# 30/3/95 Add confIDence limits and deads to 'short' output
# 3/4/95  Implement unconditional se using method of K.P.Burnham
# 3/4/95  Minor changes to full output
#

        first	<- function(vec) match(1, vec)

	jack.fill <- function(tt)
	{

	# Input:  tt is the number of capture occasions (e.g., days)
	# Output: matrix of jackknife coefficients for Burnham & Overton estimator

	# Murray Efford 30/3/95

		T1 <- tt - 1
		T2 <- tt - 2
		T3 <- tt - 3
		T4 <- tt - 4
		T5 <- tt - 5
		occ5 <- min(tt, 5)
		fcoeff <- matrix(data = 0, ncol = 5, nrow = 5)
		fcoeff[1, 1] <- T1/tt
		fcoeff[2, 1] <- (2 * tt - 3)/tt
		fcoeff[2, 2] <-  - T2^2/(tt * T1)
		fcoeff[3, 1] <- (3 * tt - 6)/tt
		fcoeff[3, 2] <-  - (3 * tt^2 - 15 * tt + 19)/(tt * T1)
		fcoeff[3, 3] <- T3^3/(tt * T1 * T2)
		fcoeff[4, 1] <- (4 * tt - 10)/tt
		fcoeff[4, 2] <-  - (6 * tt^2 - 36 * tt + 55)/(tt * T1)
		fcoeff[4, 3] <- (4 * tt^3 - 42 * tt^2 + 148 * tt - 175)/(tt * T1 * T2)
		fcoeff[4, 4] <-  - T4^4/(tt * T1 * T2 * T3)
		fcoeff[5, 1] <- (5 * tt - 15)/tt
		fcoeff[5, 2] <-  - (10 * tt^2 - 70 * tt + 125)/(tt * T1)
		fcoeff[5, 3] <- (10 * tt^3 - 120 * tt^2 + 485 * tt - 660)/(tt * T1 * T2)
		fcoeff[5, 4] <-  - (T4^5 - T5^5)/(tt * T1 * T2 * T3)
		fcoeff[5, 5] <- T5^5/(tt * T1 * T2 * T3 * T4)
		fcoeff <- fcoeff[1:occ5, 1:occ5]	# Use sub-matrix if tt < 5
		if(tt > 5) {
			fcoeff <- cbind(fcoeff, matrix(data = 0, nrow = 5, ncol = tt - 5))
		}
		fcoeff + 1	# Adding one allows: Nj<-fcoeff%*%fi
	}

	if(is.list(inp)) {
	     fi 	<- inp[[1]]
	     deads	<- inp[[2]]
	}
	else fi	<- inp

	S	<- sum(fi)
	occ	<- length(fi)
	occ5	<- min(5, occ)
	aki	<- jack.fill(occ)
	nj		<- aki %*% fi
	varnj 	<- (aki^2 %*% fi) - nj
	difnk	<- nj[2:occ5] - nj[1:(occ5 - 1)]
	b2f	<- ((aki[2:occ5,  ] - aki[1:(occ5 - 1),  ])^2) %*% fi
	test	<- (difnk/sqrt(S/(S - 1) * (b2f - difnk^2/S)))^2
	test[is.na(test)] <- 0
	test	<- rbind(test, 0)
	pk		<- 1 - pchisq((difnk/sqrt(S/(S - 1) * (b2f - difnk^2/S)))^2, 1)
	pk[difnk == 0] <- 1

	# Select jacknife on basis of test results, and interpolate estimate
	# The conditional variance calculations commented out here are superceded by the
	# unconditional calculations below

	     mt		<- first(test < 3.84) #
	     if(mt == 1) {
	          xtest	<- (nj[1] * fi[1])/(nj[1] - fi[1])
	          if(xtest > 3.84) {
	               alpha	<- (xtest - 3.84)/(xtest - test[1])
	               beta	<- 1 - alpha
	               z	<- aki[1,  ] * alpha + beta
	               N <- z %*% fi  # varN <- (z * z) %*% fi - N
	          }
	          else {
	               N	<- nj[1]
	               varN	<- varnj[1]
	          }
	     }
	     else {
	          alpha	<- (test[mt - 1] - 3.84)/(test[mt - 1] - test[mt])
	          beta	<- 1 - alpha
	          z 	<- aki[mt,  ] * alpha + aki[mt - 1,  ] * beta
	          N 	<- z %*% fi  # varN <- ((z * z) %*% fi) - N
	     }
	     k1		<- occ5 - 1
	     Pi		<- rep(1, occ5)
	     Beta 	<- pchisq(rep(3.8415, occ5), 1, pmax(test - 1, 0))
	     for(i in 2:occ5) Pi[i] <- Pi[i - 1] * (1 - Beta[i - 1])
	     Pi[1:k1]	<- Pi[1:k1] * Beta[1:k1]
	     if(occ5 < 5) Pi[occ5] <- 1 - sum(Pi[1:k1])
	     varN 	<- sum(Pi * varnj) + sum(Pi * (nj - sum(Pi * nj))^2)
             sejack <- sqrt(varN)

	# Confidence limits as in Rexstad & Burnham 1991 CAPTURE Users' Guide

	     f0		<- N - S
             if (f0>0)
                 cc1 	<- exp(1.96 * sqrt(log(1 + varN/f0^2)))
             else
                 cc1    <- NA
	     jacklcl	<- S + f0 / cc1
	     jackucl	<- S + f0 * cc1

	# Output
	     if(full) list(fi = fi, jack = N + deads, deads = deads, sejack = sejack, jacklcl
	                = jacklcl + deads, jackucl = jackucl + deads, mt = mt, aki = aki,
	               nj = data.frame(nj, se = sqrt(varnj), X2 = test, Pi = Pi, Beta = Beta))
	     else c(N + deads, sqrt(varN), jacklcl + deads, jackucl + deads, deads)
	}
##################################################

M0 <- function (ni)
## data vector ni is c(total captures, total animals, n occasions)
  { Mt1 <- ni[2]
    if (Mt1==0) c(0,NA,NA)
    else {
      theta <- optimize (f = negloglikM0, lower = log(Mt1), upper = log(1000*Mt1),
                         n = ni)$minimum
      c(exp(theta), ni[1]/exp(theta)/ni[3])
    }
  }
##################################################

Mb <- function (ui)
## data vector ui is number of new animals on each occasion
  { n <- sum(ui)
    if (n==0) c(0,NA)
    else {
      theta <- optimize (f = negloglikMb, lower = logit(0.0001), upper = logit(0.9999),
                         u = ui)$minimum
      p <- invlogit(theta)
      p <- 1 - (1-p)^length(ui)
      c(n/p, p)  ## Nhat, p
    }
  }
##################################################

Mh <- function (fi) {
## data vector fi is number of animals caught on exactly i occasions
    temp <- jack.est(fi)
    nocc <- length(fi)
    p <- sum(fi*(1:nocc)) / temp[1] / nocc
    c(temp[1], p)
}
##################################################

pfn <- function (capthist, N.estimator =  c("n", "null","zippin","jackknife")) {
    N.estimator <- tolower(N.estimator)
    N.estimator <- match.arg(N.estimator)
    ## capthist single-session only; ignoring deads
    n <- nrow(capthist)         ## number of individuals
    cts <- abind(counts(capthist, c('n','f','u')), along=1)
    ni <- cts['n', - ncol(cts)] ## drop total at end
    ui <- cts['u', - ncol(cts)] ## drop total at end
    fi <- cts['f', - ncol(cts)] ## drop total at end
    nocc <- ncol(capthist)      ## number of occasions
    estimates <- switch (N.estimator,
        n = c(n, sum(ni)/n/nocc),
        null = M0(c(sum(ni), n, nocc)),
        zippin = Mb(ui),
        jackknife = Mh(fi)
    )
    c(N=estimates[1], odds.p=odds(estimates[2]), rpsv=RPSV(capthist))
}
##################################################


# ip.secr (captdata, pfn, start=c(5,0.2,30))
