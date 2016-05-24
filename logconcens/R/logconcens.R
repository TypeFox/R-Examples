#  
#   logconcens.R
#
#   $Revision: 0.16.3 $   $Date: 2013/12/12 20:18:00 $
#
#   Code by Dominic Schuhmacher
#   Algorithm based on ideas by Lutz Duembgen, Kaspar Rufibach and Dominic Schuhmacher
#
#
#   Next steps:  -- Implement clc.fixdom and everything below in C (maybe also the repeat loop in logcon)
#
#                -- Enable domain selection based on loglikelihood comparison
#
#                -- Enable survival objects as arguments for logcon
#
#
#

# input data is n by 2 matrix x, with left and right interval endpoint in each row
# right endpoints may take value Inf;
# alternatively x may be a vector or n by 1 matrix in which case 
#

lc.control <- function(maxiter=49, move.prec=1e-5, domind1l=1, domind2r=1, force.inf=FALSE, red.thresh=NULL, check.red=TRUE, addpoints=FALSE, addeps=NULL, preweights=NULL, minw=0, show=FALSE, verbose=FALSE) {
	
  stopifnot(maxiter >= 0, move.prec > 0, domind1l >= 1, domind2r >= 1, is.logical(force.inf),
    is.null(red.thresh) | (red.thresh > 0), is.logical(check.red), is.logical(addpoints),
    is.null(addeps) | (addeps > 0), is.null(preweights) | (preweights > 0), minw >= 0, 
    is.logical(show), is.logical(verbose))
    	
return(list(maxiter=maxiter, move.prec=move.prec, domind1l=domind1l, domind2r=domind2r, 
  force.inf=force.inf, red.thresh=red.thresh, check.red=check.red, addpoints=addpoints, addeps=addeps,
  preweights=preweights, minw=minw, show=show, verbose=verbose))
}

logconcure <- function(x, p0=0, knot.prec=IQR(x[x<Inf])/75, reduce=TRUE, control=lc.control()) {
  logcon(x, adapt.p0=TRUE, p0=p0, knot.prec=knot.prec, reduce=reduce, control=control)
}

logConCens <- function(x, adapt.p0=FALSE, p0=0, knot.prec=IQR(x[x<Inf])/75, reduce=TRUE, control=lc.control()) {
  logcon(x, adapt.p0, p0, knot.prec, reduce, control)
}



#
logcon <- function(x, adapt.p0=FALSE, p0=0, knot.prec=IQR(x[x<Inf])/75, reduce=TRUE, control=lc.control()) {
#prec=.Machine$double.eps^(1/2))
  if (is.vector(x) || (is.matrix(x) && ncol(x) == 1)) {
    n <- length(x)
    x <- sort(x)
    w <- rep((1-p0)/n,n)
    # is *very* quick anyway if x is already sorted
    res <- .C("logcon_slope", sl=as.integer(FALSE), pn = as.integer(n), x = as.double(x), w = as.double(w),
                 wslr = as.double(0),  p0=as.double(p0), is_knot = as.integer(numeric(n)),
                 phi_cur = as.double(numeric(n)), phi_cur_slr = as.double(numeric(1)),
                 Fhat = as.double(numeric(n)), Fhatfin = as.double(numeric(1)), L = as.double(numeric(1)))
    res2 <- list(basedOn="exact", status=0, L=res$L, x=x, isKnot=res$is_knot, phi=res$phi_cur,
                 phislr=res$phi_cur_slr, cure=p0, Fhat=res$Fhat, Fhatfin=res$Fhatfin)
    class(res2) <- "lcdensity"
    return(res2)
  } else if (!(is.matrix(x) && ncol(x) == 2)){
    stop("argument x must be a matrix with at most two columns or a vector")
  }

  n <- dim(x)[1]
  
  wh <- (x[,1] == x[,2])
  if (length(unique(x[wh,1])) == 1) {
  	mu <- x[wh,1]
  	if (!adapt.p0 && all(x[,1] <= mu) && all(x[,2] >= mu)) {
      stop("Maximum likelihood estimator for phi does not exist if the intersection of
        all data intervals is equal to an exact observation")
  	} else if (adapt.p0 && all((x[,1] <= mu & x[,2] >= mu) | x[,2] == Inf)) {
      stop("Maximum likelihood estimator for (phi,p0) does not exist if there is an exact observation that is
      contained in all bounded intervals")  		
  	}
  }
    
  tau <- sort(unique(as.vector(x)))
  kk <- length(tau)
  k <- ifelse(is.finite(tau[kk]), kk, kk-1) # number of finite interval end points
  if ((tau[k-1]-tau[2])/knot.prec > 10000) { warning("knot.prec is very small compared to scale of data. Computations may take a long time") }
  if (kk <= 2) {
  	stop("There are only two distinct interval endpoints. If there are not exact observations at both endpoints,
  	  every log-concave function on the interval [", tau[1], ", ", tau[2], "] will do if it has the right mass.")
  }
  
  ctrl <- do.call(lc.control, as.list(control))
  maxiter <- ctrl$maxiter
  move.prec <- ctrl$move.prec
  domind1 <- ctrl$domind1l
  domind2 <- kk+1-ctrl$domind2r 
  force.inf <- ctrl$force.inf
  red.thresh <- ctrl$red.thresh
  check.red <- ctrl$check.red
  addpoints <- ctrl$addpoints
  addeps <- ctrl$addeps
  if (is.null(addeps)) {addeps=1/n^2}
  preweights <- ctrl$preweights
  if (is.null(preweights)) {preweights <- rep(1,dim(x)[1])}
  minw <- ctrl$minw
  show <- ctrl$show
  verbose <- ctrl$verbose

  if (p0 >= 1 || p0 < 0) { stop("p0 must be in [0,1)") }
  
  if (any(x[,1] > x[,2])) { stop("left interval point greater than right interval point") }
  if (any(x[,2] < tau[domind1]) || any(x[,2] < tau[domind1+1] & x[,1] < x[,2])) {
    d1 <- max(which(sapply(1:kk, function(j) {(all(x[,2] >= tau[j]) && all(x[,2] >= tau[j+1] | x[,1] == x[,2]))})))
    stop("domind1 has to be an integer from 1 to ", d1)
  }
  if (p0 == 0) {
    if (any(x[,1] > tau[domind2]) || any(x[,1] > tau[domind2-1] & x[,1] < x[,2])) {
      d2 <- min(which(sapply(1:kk, function(j) {(all(x[,1] <= tau[j]) && all(x[,1] <= tau[j-1] | x[,1] == x[,2]))})))
      stop("domind2 has to be an integer from ", d2, " to ", kk)
    }
  } else {
    if (any(is.finite(x[,2]) & x[,1] > tau[domind2]) || any(is.finite(x[,2]) & x[,1] > tau[domind2-1] & x[,1] < x[,2])) {
      d2 <- min(which(sapply(1:kk, function(j) {(all(!is.finite(x[,2]) | x[,1] <= tau[j]) && all(!is.finite(x[,2]) | x[,1] <= tau[j-1] | x[,1] == x[,2]))})))
      stop("domind2 has to be an integer from ", d2, " to ", kk)
    }
  }
  
  if (addpoints) {
  	#print(tau[1])
  	#print(tau[k])
    addtau1 <- !any((x[,1] == tau[1] & x[,2] == tau[1]))
    addtauk <- !any((x[,1] == tau[k] & x[,2] == tau[k]))
    if (addtau1) {	
  	  x <- rbind(x, c(tau[1],tau[1]))
  	  preweights <- c(preweights, addeps)
  	}
  	if (addtauk) {  
  	  x <- rbind(x, c(tau[k],tau[k]))
  	  preweights <- c(preweights, addeps)
  	}
  	n <- n + addtau1 + addtauk  
  }
  
  # subdivide[j]: does the interval [tau_j,tau_{j+1}] need subdivision
  # according to DRS11, Theorem 2.1 
  subdivide <- rep(TRUE, kk-1)
  subdivide[1:domind1] <- FALSE
  subdivide[(domind2:kk)-1] <- FALSE
  if (!any(x[,2] == tau[domind1]))  { subdivide[domind1+1] <- FALSE }
  if (!any(x[,1] == tau[domind2]))  { subdivide[domind2-2] <- FALSE }
  for (j in (domind1+1):(domind2-2)) {
    if (!any(x[,1] <= tau[j] & x[,2] >= tau[j+1])) { subdivide[j] <- FALSE }
  }
  repeat {
    fixres <- clc.fixdom(x, preweights, minw, p0, adapt.p0, reduce, red.thresh, check.red, force.inf, tau, subdivide, domind1, domind2,
                         maxiter=maxiter, knot.prec=knot.prec, move.prec=1e-5, show=show, verbose=verbose)
    p0 <- fixres$p0new                    
    if (fixres$status >= 0) break
    if (fixres$status == -1) {
      domind1 <- domind1+1
      subdivide[domind1] <- FALSE
      if (!any(x[,2] == tau[domind1]))  { subdivide[domind1+1] <- FALSE }
      message("Domain reduced on left hand side! New indices: ", domind1, " ",domind2)
    } else {
      domind2 <- domind2-1
      fixres$phislr <- -Inf  # for aesthetic reasons
      subdivide[domind2-1] <- FALSE
      if (!any(x[,1] == tau[domind2]))  { subdivide[domind2-2] <- FALSE }
      message("Domain reduced on right hand side! New indices: ", domind1, " ",domind2)
    }
    if (domind2-domind1 == 1) {
      stop("Domain reduction led to a domain limited by two consecutive interval endpoints. If there are not
        exact observations at both endpoints, every log-concave function on the interval [", tau[domind1], ", ", 
        tau[domind2], "] will do if it has the right mass.")
    }
  }

  mm <- length(fixres$tplus)
  m <- ifelse(is.finite(fixres$tplus[mm]), mm, mm-1)

  if (show) {
    dev.new()
    uplim <- tau[k] + ifelse(k == kk, 0, 0.15*(tau[k]-tau[1]))
    plot(fixres$tplus[1:m], fixres$phi,xlim=c(tau[1],uplim), xlab="tau and t", ylab="phi", type="l",col=4,lwd=2)
    if (fixres$phislr > -Inf) { lines(c(fixres$tplus[m],uplim),
                               c(fixres$phi[m], fixres$phi[m]+fixres$phislr*(uplim-fixres$tplus[m])),col=4,lwd=2) }
    abline(v=fixres$tplus[as.logical(fixres$isKnot)], lty=3, col=4)
    rug(fixres$tplus, ticksize=0.02)
    rug(tau, ticksize=0.04, lwd=1)
  }

  cure.range = NA
  phislr.range = NA
  if (mm > m & adapt.p0) {
    rint <- -exp(fixres$phi[m])/fixres$phislr
    phislrleft <- (fixres$phi[m]-fixres$phi[m-1])/(fixres$tplus[m]-fixres$tplus[m-1])
    p0max <- p0+rint
    p0min <- max(0, p0max+exp(fixres$phi[m])/phislrleft)
    phislrmax <- min(-exp(fixres$phi[m])/p0max, phislrleft)    
    phislrmin <- -Inf
    cure.range = c(p0min,p0max)
    phislr.range = c(phislrmin,phislrmax)
  } else if (kk > k & domind2==k & adapt.p0) {
    # rint <- 0
    # m is correct, since tplus[m] == tau[k] here
    phislrleft <- (fixres$phi[m]-fixres$phi[m-1])/(fixres$tplus[m]-fixres$tplus[m-1])
    p0max <- p0
    p0min <- ifelse(phislrleft < 0, max(0,p0max+exp(fixres$phi[m])/phislrleft), 0)
    phislrmax <- min(-exp(fixres$phi[m])/p0max, phislrleft)
    phislrmin <- -Inf
    cure.range = c(p0min,p0max)
    phislr.range = c(phislrmin,phislrmax)  	
  } 

  res <- list(basedOn="censored", status=fixres$status, x=x, tau=tau, domind1=domind1, domind2=domind2,
              tplus=fixres$tplus, isKnot=fixres$isKnot, phi=pmax(fixres$phi,-1e100), phislr=fixres$phislr,
              phislr.range=phislr.range, cure=p0, cure.range=cure.range,
              Fhat=fixres$Fhat, Fhatfin=fixres$Fhatfin)
  class(res) <- "lcdensity"
  return(res)
}
# 10/12/2013:
# pmax with -1e100 above is a quick an dirty fix for the case that phis at the boundary are -Inf
# This can (apparantly) only happen with (one or ?) several alligned [L_i, Inf] intervals to the right of
# everything


# x is the original data (nx2 matrix)
# tau is the ordered unique set of interval endpoints
# domind1 and domind2 are the indices (for the tau-vector) of the left end right
#    domain endpoints that shall be considered
clc.fixdom <- function(x, preweights=rep(1,dim(x)[1]), minw=0, p0, adapt.p0 = FALSE, reduce=TRUE, red.thresh=NULL, check.red=TRUE, force.inf=FALSE, tau, subdivide, domind1, domind2, maxiter=60, knot.prec=IQR(x[x<Inf])/75,
                       move.prec=1e-5, show=TRUE, verbose=FALSE) {
  if (dim(x)[2] != 2) { stop("x should be an n times 2 matrix") }
  n <- dim(x)[1]
  
  needsl <- !is.finite(tau[domind2])  # for checking if we have to use slope code
                                      # i.e. there are right endpoints = Inf and the slope
                                      # has not dropped out yet
  canReduceL <- reduce
  canReduceR <- reduce & !force.inf
  if (any(x[,2] == tau[domind1]) || any(x[,2] == tau[domind1+1] & x[,1] < x[,2])) {
    canReduceL <- FALSE
  }
  if (p0 == 0) {
    if (any(x[,1] == tau[domind2]) || any(x[,1] == tau[domind2-1] & x[,1] < x[,2])) {
      canReduceR <- FALSE
    }
  } else {
    if (any(is.finite(x[,2]) & x[,1] == tau[domind2]) || any(is.finite(x[,2]) & x[,1] == tau[domind2-1] & x[,1] < x[,2])) {
      canReduceR <- FALSE
    }
  }
  if (verbose) {
    print(paste("canReduce:", canReduceL, canReduceR))
  }

  # tplus is generated based on [tau[domind1],tau[domind2]] not on [tau[1],tau[kk]]
  # Note that this assumes that subdivide[j]=FALSE for j<=domind1-1 and j>=domind2
  # if this is not the case all hell breaks loose (catered for when reducing domains in logcon)
  tplus <- c(tau[c(domind1:domind2)], unlist(lapply(which(subdivide), subdivisor, tau=tau, eps=knot.prec)))
  tplus <- sort(tplus)  # there should be a faster way, since tau and tplus are already
                        # sorted and tau usually much shorter than tplus (also since we are computing the ranks)
  dtplus <- diff(tplus)
  mm <- length(tplus)
  m <- ifelse(needsl, mm-1, mm)
  xx <- cbind(pmax(x[,1],tplus[1]), pmin(x[,2],tplus[mm]))
  xx[,1] <- pmin(xx[,1],tplus[mm])
  xindex <- matrix(sapply(xx, function(xval) {which(tplus == xval)}), n, 2)
  # if p0 > 0 certain data intervals might lie to the right of our current domain
  # we have to store which intervals (they may show up in xx as an interval
  # of length zero (indistinguishable from exact data points) or as an interval
  # with left endpoint > right endpoint) => analogously for xindex
  xindex[(x[,1] >= tplus[mm] & x[,1] < x[,2]), ] <- 0
  # index numbers of xx values in terms of tplus
  # Note that tplus consists exactly of values from xx (floating point representations are exactly equal),
  # so no problem with the ==
  rightinf <- (x[,2] == Inf)
  # needed for deciding whether to add p0 for prob that rv in \tX_i if domain is finite
  
  # initialize phi, note: phislr is the right hand slope (=-Inf) in the non-slope case
  # isKnot and phi both have length m (NOT mm in general)
  isKnot <- rep(FALSE,m)
  if (needsl) {
    isKnot[1] <- TRUE
    phislr <- -1/(tplus[m]-tplus[1])
    phi <- -log(tplus[m]-tplus[1]) + phislr*(tplus[1:m]-tplus[1])
  } else {
    isKnot[c(1,m)] <- TRUE
    phislr <- -Inf
    phi <- rep(-log(tplus[m]-tplus[1]), m)
  }

  if (show) {
    dev.new(width=16,height=9)
    par(mfrow=c(5,10), mai=rep(0.1,4))
    plot(tplus[1:m],phi,type="l")
    rug(tau[domind1:domind2])
  }
  iter <- 1
  move <- 1
  while (iter <= maxiter && move > move.prec) {
    ww <- GetWeights(preweights, tplus, p0, xindex, rightinf, phi, phislr, needsl)
    # With the new version from 15/11/2013 of the function J10 (and J00) the next block should not be entered anymore!
    if (!all(is.finite(unlist(ww)))) {
      if(!is.finite(ww$w[1]) && canReduceL) {
        message("Cannot compute leftmost weight. Emergency domain reduction to the left. phi[1] is ", phi[1])
        return(list(status=-1, tplus=tplus, isKnot=asres$isKnot, phi=phi, phislr=phislr, p0new=p0))
      } else if (needsl && (!is.finite(ww$wslr) || !is.finite(ww$w[m])) && canReduceR) {
        message("Cannot compute rightmost weights. Emergency domain reduction to the right. phi[m] is ",
                phi[m], ", phislr is ", phislr)
        return(list(status=-2, tplus=tplus, isKnot=asres$isKnot, phi=phi, phislr=phislr, p0new=p0))
      } else if (!needsl && !is.finite(ww$w[m]) && canReduceR) {
        message("Cannot compute rightmost weight. Emergency domain reduction to the right. phi[m] is ", phi[m])
        return(list(status=-2, tplus=tplus, isKnot=asres$isKnot, phi=phi, phislr=phislr, p0new=p0))
      } else {
        print(cat("w is \n", ww$w, "\n", "phi is \n", phi, "\n"))
        warning("Cannot compute all weights, but no further domain reduction possible.")
        return(list(status=2, tplus=tplus, isKnot=asres$isKnot, phi=phi, phislr=phislr, p0new=p0))
      }
    }

    if (verbose) {
      message("Minimal weight: ", min(ww$w), "   wslr: ", ww$wslr)
    }
    if ((min(ww$w) < 1e-10) && verbose) { warning("small weights in call to logcon_slope") }
    if ((needsl & ww$wslr < 1e-10) && verbose) { warning("small slope weight in call to logcon_slope") }
    # The weights should always be positive (barring numerical problems)
    # so I skip the following; same for the block commented out below
    # wpos <- (w>0)
    # ww <- w[wpos]
    # ww <- ww/sum(ww)
    # tt <- tplus[wpos]
    ww$w[1] <- max(ww$w[1], minw)
    if (needsl) {
      ww$wslr <- max(ww$wslr, minw)
    } else {
      ww$w[m] <- max(ww$w[m], minw)
    }
    ww$w <- (1-p0)*ww$w/sum(ww$w)
    ww$wslr <- (1-p0)*ww$wslr
    # if p0 = 0 the above normalizations are not necessary (checked by manual computation)
    # otherwise yes
    asres <- .C("logcon_slope", sl=as.integer(needsl), pn = as.integer(m), x = as.double(tplus[1:m]),
                 w = as.double(ww$w),
                 wslr = as.double(ww$wslr), p0=as.double(p0), is_knot = as.integer(numeric(m)),
                 phi_cur = as.double(numeric(m)), phi_cur_slr = as.double(numeric(1)),
                 Fhat = as.double(numeric(m)), Fhatfin = as.double(numeric(1)), L = as.double(numeric(1)))
    J00phimax <- J00( pmax(phi[1:(m-1)], asres$phi_cur[1:(m-1)]), pmax(phi[2:m], asres$phi_cur[2:m]) )
    J00phimin <- J00( pmin(phi[1:(m-1)], asres$phi_cur[1:(m-1)]), pmin(phi[2:m], asres$phi_cur[2:m]) )
    move1 <- sum(dtplus[1:m-1] * (J00phimax - J00phimin))
    move <- move1 + ifelse(needsl, exp(max(phi[m],asres$phi_cur[m])) * abs(1/phislr-1/asres$phi_cur_slr), 0)
    # move <- max(c(abs(asres$phi_cur - phi), ifelse(needsl, 1e-4*abs(asres$phi_cur_slr - phislr), 0)))
                                  # 1e-4 is a fudge factor as slope movement gets to important otherwise
                                  # on the other hand this means that our slope precision is for digits
                                  # worse than precision of phi-values
                                  #
                                  # Untenstehendes bezog sich auf Summe
                                  # Scheint mir nicht so das geeignete Kriterium 25/08/11
                                  # besseres was mit Integral, aber das hier ist natürlich einfacher zu berechnen
                                  # Etwas positives hier ist natürlich, dass bei abschmierendem Rand
                                  # sicher nicht ploetzlich der move klein werden kann
    if (verbose) {
      message("Move: ", move, "   iter: ", iter)
    }
    isKnot <- asres$is_knot
    phi <- asres$phi_cur
    phislr <- asres$phi_cur_slr
    if (show) {
      plot(tplus[1:m],phi,type="l")
      rug(tau[domind1:domind2])
    }
    J00phi <- J00(phi[1:(m-1)],phi[2:m])
    # restriction of domain
    leftint <- (tplus[2]-tplus[1]) * J00(phi[1],phi[2])
    rightint <- ifelse( needsl, exp(phi[m])/(-phislr),  (tplus[m]-tplus[m-1]) * J00(phi[m-1],phi[m]))
    leftslope <- (phi[2]-phi[1]) / (tplus[2]-tplus[1])
    rightslope <- ifelse( needsl, phislr, (phi[m]-phi[m-1]) / (tplus[m]-tplus[m-1]) )
    if (is.null(red.thresh)) {
      red.thresh <- ifelse(check.red, 0.01, min(5e-4, (1e-2)/(domind2-domind1)))
    }
    if (verbose) {
      cat("left int:", leftint, "\n")
      cat("right int:", rightint, "\n")
      cat("left slope:", leftslope, "\n")
      cat("right slope:", rightslope, "\n")
      cat("left phi:", phi[1], "\n")
      cat("right phi:", phi[m], "\n")
      cat("compared to:", red.thresh, "\n")
    }
#    lefttest <- (leftslope > 1e3)
#    righttest <- (rightslope < -1e3)
    lefttest <- (leftint < red.thresh)   ###### ADD AND MOVE < XXX
    #lefttest <- righttest <- (move < 1e-3)
      # eigentlich sollte das 10^(-5)/n sein, aber das dauert ewig!!
      # evtl. leftint und rightint durch entsprechende Intervallaengen der aeussersten Intervalle teilen
    if (lefttest && verbose) { message("Leftmost integral small") }
    righttest <- (rightint < red.thresh)    ###### ADD AND MOVE < XXX
    if (righttest && verbose) { message("Rightmost integral small") }
    lefttest <- (lefttest && canReduceL)
    righttest <- (righttest && canReduceR)
    if (lefttest && righttest) {
      if (leftint <= rightint) {righttest <- FALSE} else {lefttest <- FALSE}
    }
    if (lefttest) {
      if (!check.red) {
      	return(list(status=-1, tplus=tplus, isKnot=isKnot, phi=phi, phislr=phislr, p0new=p0))
      } else {
        fstart <- rep(2,n)
        fend <- pmin(xindex[,2],m) 
        critfun <- function(i) {
          integ <- ifelse(fstart[i] < fend[i], sum(dtplus[fstart[i]:(fend[i]-1)] * J00phi[fstart[i]:(fend[i]-1)]), 0)
          crit0 <- ifelse(xindex[i,1] == 1, 1/(integ + (x[i,2] == Inf) * (p0 - ifelse(needsl, exp(phi[fend[i]])/phislr, 0))), 0)
          return(crit0)
        }
        crit <- unlist(lapply(1:n, critfun))
        # print(paste("Left",mean(crit)))
        if (mean(crit) <= 1) { 
          return(list(status=-1, tplus=tplus, isKnot=isKnot, phi=phi, phislr=phislr, p0new=p0))
        }
      }
    }
    if (righttest) {
      if (!check.red) {
      	return(list(status=-2, tplus=tplus, isKnot=isKnot, phi=phi, phislr=phislr, p0new=p0))
      } else {    	
      	# The following is not so great with adapt.p0
      	# not clear why, (means just that phislr.range may be -Inf to -Inf or -Inf to very small
      	# we catch this in plot)
        if (needsl) {
          fstart <- xindex[,1]
          fend <- rep(mm-2,n) 
          critfun <- function(i) {
          	integ <- ifelse(fstart[i] < fend[i], sum(dtplus[fstart[i]:(fend[i]-1)] * J00phi[fstart[i]:(fend[i]-1)]), 0)
          	crit0 <- ifelse(x[i,2] == Inf, 1/(integ + p0), 0)
          # included here also (x[i,2] == Inf) * p0 ) unlike Lutz's suggestion, would like to reduce
      	  # domain on the right also if p0>0, since the main argument for domain reduction is numerical stability
      	  # and this is equally valid if p0>0
          	return(crit0)
          }
          crit <- unlist(lapply(1:n, critfun))
          #print(paste("Right+",mean(crit)))
          if (mean(crit) <= 1) { 
            return(list(status=-2, tplus=tplus, isKnot=isKnot, phi=phi, phislr=phislr, p0new=p0))
          }	
        } else {
          fstart <- xindex[,1]
          fend <- rep(m-1,n) 
          critfun <- function(i) {
          	integ <- ifelse(fstart[i] < fend[i], sum(dtplus[fstart[i]:(fend[i]-1)] * J00phi[fstart[i]:(fend[i]-1)]), 0)
          	crit0 <- ifelse(xindex[i,2] >= m & xindex[i,1] < m, 1/(integ + (x[i,2] == Inf) * p0), 0)
          	return(crit0)
          }
          crit <- unlist(lapply(1:n, critfun))
          # print(paste("Right",mean(crit)))
          if (mean(crit) <= 1) { 
            return(list(status=-2, tplus=tplus, isKnot=isKnot, phi=phi, phislr=phislr, p0new=p0))
          }	          
        }
      }
    }
    if (adapt.p0) {  
        m1fun <- function(q) {
          m1s <- subint <- rep(0,n)
          #print(tplus)
          #print(x)
          #print(xx)
          #print(xindex)
          for (i in 1:n) {
          	if(x[i,2] == Inf) {
          	  # cat("\n subint \n",subint)
              fstart <- xindex[i,1]
              # The case of (automatically right-infinite) cut off intervals:
              if (fstart == 0) {
                subint[i] <- 0
              } else {
                fend <- m
                #cat(i,fstart,fend,length(dtplus),length(J00phi),"\n")
                subint[i] <- ifelse(fstart == fend, 0, sum(dtplus[fstart:(fend-1)] * J00phi[fstart:(fend-1)]))
                #cat(i,subint[i],"\n")
                if (needsl) {  
                  subint[i] <- subint[i] - exp(phi[fend])/phislr
                }
              }
              m1s[i] <- 1/(subint[i]+q)
            }
          }
          return(mean(m1s)-1)
        }
        canReduceR <- reduce & !force.inf
        if (m1fun(0) <= 0) {
          p0 <- 0;
          if (any(x[,1] == tau[domind2]) || any(x[,1] == tau[domind2-1] & x[,1] < x[,2])) {
            canReduceR <- FALSE
          }
          cat("0 is opt. m1fun(0) =",m1fun(0), "\n")	
        } else {
          solu <- uniroot(m1fun,c(0,1))
          cat("new p0:", solu$root, "\n")
          p0 <- solu$root
          if (any(is.finite(x[,2]) & x[,1] == tau[domind2]) || any(is.finite(x[,2]) & x[,1] == tau[domind2-1] & x[,1] < x[,2])) {
            canReduceR <- FALSE
          }
        }
        if (verbose) {
          print(paste("canReduceR:", canReduceR))
        }
    }
    iter <- iter + 1
  }
  status <- ifelse(move <= move.prec, 0, 1)
  return(list(status=status, tplus=tplus, isKnot=isKnot, phi=phi, phislr=phislr, Fhat=asres$Fhat, Fhatfin=asres$Fhatfin, p0new=p0))
}


# GetWeights HAS NUMERICAL PROBLEMS if phi is very small (< -600) at boundary points;
# 03/12/13 fixed
GetWeights <-
function (preweights, tplus, p0, xindex, rightinf, phi, phislr, needsl)
# isKnot not used, but would be the goal for more efficient computation
# "from kink to kink"
{
  n <- dim(xindex)[1]
  mm <- length(tplus)
  m <- ifelse(needsl, mm-1, mm)
  dtplus <- diff(tplus)
  J00phi <- J00(phi[1:(m-1)],phi[2:m])
  J01phi <- J10(phi[2:m],phi[1:(m-1)])
  J10phi <- J10(phi[1:(m-1)],phi[2:m])

  PXinTX <- rep(0,n)
  for (i in 1:n) {
    # xindex[i,1] == 0 can only happen when p0 > 0 and an interval has been cut-off
    # -Inf is not used for any calculation just to make errors more obvious
    if (xindex[i,1] == 0) {
      PXinTX[i] <- -Inf
    } else {
      fstart <- xindex[i,1]
      fend <- min(xindex[i,2], m)
      PXinTX[i] <- ifelse(fstart == fend, 0, sum(dtplus[fstart:(fend-1)] * J00phi[fstart:(fend-1)]))
      if (xindex[i,2] == m+1) {  
        PXinTX[i] <- PXinTX[i] - exp(phi[fend])/phislr
      }
      if (rightinf[i]) {
        PXinTX[i] <- PXinTX[i] + p0
      }
    }
  }
  w <- rep(0,m)  # weights for the positions in tplus except for a possible last position = Inf
                 # remember: length(tplus) is mm
  wslr <- 0      # weight for the slope
  for (j in 1:m) {
    fstart <- xindex[,1]
    fend <- pmin(xindex[,2], m)
    # subw has length n and contains the summands of (18) in DHR07
    #
    subw <- (xindex[,1] == j & xindex[,2] == j)
    subw <- subw + ifelse(j-1 >= fstart & j <= fend, dtplus[j-1]*J01phi[j-1]/PXinTX, 0) 
    subw <- subw + ifelse(j >= fstart & j+1 <= fend, dtplus[j]*J10phi[j]/PXinTX, 0)
    # evtl exp(log(a)-log(b)) statt a/b
    w[j] <- sum(preweights*subw)/sum(preweights)
  }
  if (needsl) {
    subw <- ifelse(xindex[,2] == m+1, exp(phi[m])/(-phislr * PXinTX), 0)
    w[m] <- w[m] + sum(preweights*subw)/sum(preweights)
    subw <- ifelse(xindex[,2] == m+1, exp(phi[m])/(phislr^2 * PXinTX), 0)
    wslr <- sum(preweights*subw)/sum(preweights)
  }
  return(list(w=w, wslr=wslr))
}


loglike <- function(lcd) {
  if (!inherits(lcd, "lcdensity")) {
    stop(paste("argument", sQuote(deparse(substitute(lcd))), "is not of class", sQuote("lcdensity")))
  }
  x <- lcd$x
  tplus <- lcd$tplus
  phi <- lcd$phi
  phislr <- lcd$phislr
  p0 <- ifelse(is.null(lcd$cure), 0, lcd$cure)

  n <- dim(x)[1]
  mm <- length(tplus)
  m <- ifelse(is.finite(lcd$tplus[mm]), mm, mm-1)
  dtplus <- diff(tplus)  
  J00phi <- J00(phi[1:(m-1)],phi[2:m])
  
  xx <- cbind(pmax(x[,1],tplus[1]), pmin(x[,2],tplus[mm])) 
  xx[,1] <- pmin(xx[,1],tplus[mm])
  xindex <- matrix(sapply(xx, function(xval) {which(tplus == xval)}), n, 2)
  xindex[(x[,1] >= tplus[mm] & x[,1] < x[,2]), ] <- 0
  contrib <- rep(0,n)
  for (i in 1:n) {
    fstart <- xindex[i,1]
    # The case of (automatically right-infinite) cut-off intervals:
    if (fstart == 0) {
      contrib[i] <- log(p0)
    } else if (fstart == xindex[i,2]) {
      contrib[i] <- phi[fstart]
    } else {
      fend <- min(xindex[i,2], m)
      contrib[i] <- ifelse(fstart == fend, 0, sum(dtplus[fstart:(fend-1)] * J00phi[fstart:(fend-1)]))
      if (xindex[i,2] == m+1) {  
        contrib[i] <- contrib[i] - exp(phi[fend])/phislr
      }
      if (x[i,2] == Inf) {
        contrib[i] <- contrib[i] + p0
      }
      contrib[i] <- log(contrib[i])
    }
  }
  return(mean(contrib))
}



cure.profile <- function(x, p0grid=seq(0,0.95,0.05), knot.prec=IQR(x[x<Inf])/75, reduce=TRUE, control=lc.control()) {
  N <- length(p0grid)
  llike <- rep(0,N)
  for (i in 1:N) {
    lcd <- logcon(x, adapt.p0=FALSE, p0=p0grid[i], knot.prec=knot.prec, reduce=reduce, control=control)
    llike[i] <- loglike(lcd)
    cat("p0 =", p0grid[i], "    loglike =", llike[i], "\n")
  }
  return(list(p0hat=p0grid[which.max(llike)], loglike=llike))
}


# ------------------------------------------
# functions from logcondens
# ------------------------------------------

J00 <- function (x, y, v = 1) 
{
    m <- length(x)
    z <- exp(x)
    d <- y - x
    II <- (1:m)[abs(d) > 0.005 & d < 200]
    z[II] <- z[II] * (exp(v * d[II]) - 1)/d[II]
    II <- (1:m)[abs(d) <= 0.005]
    z[II] <- z[II] * (v + d[II] * (v/2 + d[II] * (v/6 + d[II] * 
        (v/24 + d[II] * v/120))))
    II <- (1:m)[abs(d) > 0.005 & d >= 200]
    z[II] <- (exp(y[II]) - exp(x[II]))/d[II]
    return(z)
}

J10 <- function (x, y) 
{
    m <- length(x)
    z <- exp(x)
    d <- y - x
    II <- (1:m)[abs(d) > 0.01 & d < 200]
    z[II] <- z[II] * (exp(d[II]) - 1 - d[II])/(d[II]^2)
    II <- (1:m)[abs(d) <= 0.01]
    z[II] <- z[II] * (1/2 + d[II] * (1/6 + d[II] * (1/24 + d[II] * 
        (1/120 + d[II]/720))))
    II <- (1:m)[abs(d) > 0.01 & d >= 200]
    z[II] <- (exp(y[II]) - exp(x[II]) * (d[II] + 1))/(d[II]^2)
    # The case d large (usually from d > 715 on or so) is
    # the only one that poses a practical problem for us
    return(z)
}



# ------------------------------------------
# Little Helpers
# ------------------------------------------

# Create shortest arithmetic sequence from tau[j] to tau[j+1]<Inf with subsequent elements
# at most eps apart
# j has to be from 1 to k-1 (not kk-1, since tau[kk] may be =Inf)
subdivisor <- function(j,tau,eps=0.01) {
  a <- tau[j]
  b <- tau[j+1]
  N <- ceiling((b-a)/eps)+1
  return(seq(a, b, length.out=N)[-c(1,N)])
}

# compute the phi-values (linear interpolations)
# at the points generated above
phidivisor <- function(j,tau,phi,eps=0.01) {
  a <- tau[j]
  b <- tau[j+1]
  N <- ceiling((b-a)/eps)+1
  return(seq(phi[j], phi[j+1], length.out=N)[-c(1,N)])
}


# ------------------------------------------
# Constructors and print and plot methods (also general)
# ------------------------------------------

# Pretty plot of x with t-grid
plotint <- function(x, knot.prec=IQR(x[x<Inf])/75, imarks=NULL) {
  finiteonly <- is.finite(max(x[,2]))
  finmax <- max(x[is.finite(x)])
  upper <- ifelse(finiteonly, max(x[,2]), finmax + 0.3*(finmax - min(x[,1])))
  upper2 <- ifelse(finiteonly, max(x[,2]), finmax + 0.5*(finmax - min(x[,1])))
  plot(c(min(x[,1]),upper),c(0,-1),type="n",axes=FALSE,xlab="tau and t",ylab="")
  box()
  axis(1)

  n <- dim(x)[1]
  tau <- sort(unique(as.vector(x)))
  kk <- length(tau)
  k <- ifelse(is.finite(tau[kk]), kk, kk-1) # number of finite interval end points

  # subdivide[j]: does the interval [tau_j,tau_{j+1}] need subdivision
  # according to DRS11, Theorem 2.1 
  subdivide <- rep(TRUE, kk-1)
  subdivide[1] <- FALSE
  subdivide[kk-1] <- FALSE
  if (!any(x[,2] == tau[1]))  { subdivide[2] <- FALSE }
  if (!any(x[,1] == tau[kk]))  { subdivide[kk-2] <- FALSE }
  for (j in 1:(kk-1)) {
    if (!any(x[,1] <= tau[j] & x[,2] >= tau[j+1])) { subdivide[j] <- FALSE }
  }
  tplus <- unlist(lapply(which(subdivide), subdivisor, tau=tau, eps=knot.prec))
  tplus <- sort(c(tau, tplus))

  rug(tplus, ticksize=0.02)
  rug(tau, ticksize=0.04, lwd=1)
  xx <- cbind(sort(x[,1]), x[,2][order(x[,1])])
  if (!finiteonly) {
    xx[,2][xx[,2] == Inf] <- upper2
  }
  segments(xx[,1],seq(0,-1,length.out=n),xx[,2],seq(0,-1,length.out=n))
  wh <- which(xx[,1]==xx[,2])
  points(xx[wh],seq(0,-1,length.out=n)[wh])
  
  if (!is.null(imarks)) {
    points(imarks[order(x[,1])], seq(0,-1,length.out=n), pch="x", cex=max(0.2,0.9-n/300))	
  }
  
  invisible()
}


# plot method for lcdensity class
plot.lcdensity <- function(x, type = c("log-density", "density", "CDF", "survival"), sloperange=TRUE, kinklines=TRUE,
                           kinkpoints=FALSE, xlim=NULL, ylim=NULL, ...)
{
  lcd <- x
  
  type <- match.arg(type)
  
  if (lcd$basedOn == "exact")
  {
    xvec <- lcd$x
    if (is.null(xlim)) { xlim <- range(xvec) }
    stp <- diff(range(xlim))/300  # 300 is approximate number of points computed in density/CDF/survival plots
    xiknots <- xvec[as.logical(lcd$isKnot)]
    etaknots <- lcd$phi[as.logical(lcd$isKnot)]
    xi <- c(xiknots, unlist(lapply(seq(1, sum(lcd$isKnot)-1), subdivisor, tau=xiknots, eps=stp)))
    eta <- c(etaknots, unlist(lapply(seq(1, sum(lcd$isKnot)-1), phidivisor, tau=xiknots, phi=etaknots, eps=stp)))
    eta <- eta[order(xi)]
    xi <- sort(xi)
    knotlist <- match(xiknots,xi)
    n <- length(xi)
    
    if (type == "log-density") {
      if (is.null(ylim)) { ylim <- range(lcd$phi) }
      plot(xi, eta, xlim=xlim, ylim=ylim, type="l", col=4, lwd=2,
           main=paste("log-density", ifelse(is.null(lcd$cure),"",paste(" (cure param. = ", lcd$cure, ")", sep="")), sep=""),
           xlab="x", ylab="phi", ...)
      if (kinkpoints) {
        points(xiknots, etaknots, pch=19, ...)
      }
      if (kinklines) {
        abline(v=xiknots, lty=3, col=4, ...)
      }
      rug(xvec)
    }

    if (type == "density") {
      if (is.null(ylim)) { ylim <- c(0,exp(max(lcd$phi))) }
      plot(xi, exp(eta), xlim=xlim, ylim=ylim, type="l", col=4, lwd=2,
           main=paste("density", ifelse(is.null(lcd$cure),"",paste(" (cure param. = ", lcd$cure, ")", sep="")), sep=""),
           xlab="x", ylab="f", ...)
      lines(c(xi[1]-0.04*diff(xlim), xi[1]), c(0,0), col=4, lwd=2, ...)
      lines(c(xi[n], xi[n]+0.04*diff(xlim)), c(0,0), col=4, lwd=2, ...)  
      if (kinkpoints) {
        points(xiknots, exp(etaknots), pch=19, ...)
      }
      if (kinklines) {
        abline(v=xiknots, lty=3, col=4, ...)
      }
      rug(xvec)
    }
    
    # the following is for type="CDF" and "survival"
    if (is.null(ylim)) { ylim <- c(0,1) }
    subint <- diff(xi) * J00(eta[-n],eta[-1])
    Fhat <- c(0,cumsum(subint))
    if (type == "CDF") {
      plot(xi, Fhat, xlim=xlim, ylim=ylim, type="l", col=4, lwd=2,
           main=paste("Fhat", ifelse(is.null(lcd$cure),"",paste(" (cure param. = ", lcd$cure, ")", sep="")), sep=""),
           xlab="x", ylab="F", ...)
      lines(c(xi[1]-0.04*diff(xlim), xi[1]), c(0,0), col=4, lwd=2, ...)
      lines(c(xi[n], xi[n]+0.04*diff(xlim)), c(1-lcd$cure,1-lcd$cure), col=4, lwd=2, ...)  
      if (kinkpoints) {
        points(xiknots, Fhat[knotlist], pch=19, ...)
      }
      if (kinklines) {
        abline(v=xiknots, lty=3, col=4, ...)
      }
      rug(xvec)
    }
    if (type == "survival") {
      plot(xi, 1-Fhat, xlim=xlim, ylim=ylim, type="l", col=4, lwd=2,
           main=paste("Shat", ifelse(is.null(lcd$cure),"",paste(" (cure param. = ", lcd$cure, ")", sep="")), sep=""),
           xlab="x", ylab="S", ...)
      lines(c(xi[1]-0.04*diff(xlim), xi[1]), c(1,1), col=4, lwd=2, ...)
      lines(c(xi[n], xi[n]+0.04*diff(xlim)), c(lcd$cure,lcd$cure), col=4, lwd=2, ...)  
      if (kinkpoints) {
        points(xiknots, 1-Fhat[knotlist], pch=19, ...)
      }
      if (kinklines) {
        abline(v=xiknots, lty=3, col=4, ...)
      }
      rug(xvec)
    }
#
# 
  } else if (lcd$basedOn == "censored")
  {
    tplus <- lcd$tplus
    mm <- length(tplus)
    m <- ifelse(is.finite(tplus[mm]), mm, mm-1)

    if (is.null(xlim)) { xlim <- range(tplus[1:m]) + c(0,ifelse(lcd$phislr > -Inf, 0.2*diff(range(tplus[1:m])), 0)) }
    stp <- diff(range(xlim))/300  # 300 is approximate number of points computed in density/CDF/survival plots
    xiknots <- tplus[1:m][as.logical(lcd$isKnot)]
    etaknots <- lcd$phi[as.logical(lcd$isKnot)]
    if (sum(lcd$isKnot) > 1) {
      xi <- c(xiknots, unlist(lapply(seq(1, sum(lcd$isKnot)-1), subdivisor, tau=xiknots, eps=stp)))
      eta <- c(etaknots, unlist(lapply(seq(1, sum(lcd$isKnot)-1), phidivisor, tau=xiknots, phi=etaknots, eps=stp)))
      eta <- eta[order(xi)]
      xi <- sort(xi)
    } else {
      xi <- xiknots
      eta <- etaknots
    }
    knotlist <- match(xiknots,xi)
    n <- length(xiknots)
    
    # for the grey areas (since lcd$phislr.range may be == -Inf and we still want to draw something)
    phislrlow <- max(lcd$phislr.range[1], -1e8)
    phislrhigh <- max(lcd$phislr.range[2], -1e8)
    # since by numerical problems we can get -Inf -Inf for lcd$phislr.range[2]

    if (type == "log-density") {
      if (is.null(ylim)) { ylim <- range(lcd$phi) + c(ifelse(lcd$phislr > -Inf, min(-0.6*diff(range(lcd$phi)), lcd$phislr * 0.2*diff(range(tplus[1:m]))), 0), 0) }
      plot(xi, eta, xlim=xlim, ylim=ylim, type="l", col=4, lwd=2,
           main=paste("log-density", ifelse(is.null(lcd$cure),"",paste(" (cure param. = ", lcd$cure, ")", sep="")), sep=""),
           xlab="x", ylab="phi", ...)
      if (sloperange & !is.na(lcd$cure.range[1])) {
         polygon(c(tplus[m], xlim[2]+0.04*diff(xlim), xlim[2]+0.04*diff(xlim)), c(lcd$phi[m], lcd$phi[m]+phislrlow*(xlim[2]+0.04*diff(xlim)-tplus[m]), lcd$phi[m]+phislrhigh*(xlim[2]+0.04*diff(xlim)-tplus[m])), border=NA, col=grey(0.9), lwd=1, ...)
         box()      	
      }
      if (lcd$phislr > -Inf) {
        lines(c(xiknots[n], xlim[2]+0.04*diff(xlim)), c(etaknots[n], etaknots[n]+lcd$phislr*(xlim[2]+0.04*diff(xlim)-xiknots[n])), col=4, lwd=2, ...) 
      }
      if (kinkpoints) {
        points(xiknots, etaknots, pch=19, ...)
      }
      if (kinklines) {
        abline(v=xiknots, lty=3, col=4, ...)
      }
      rug(lcd$tplus, ticksize=0.02)
      rug(lcd$tau, ticksize=0.04, lwd=1)
    }

    if (type == "density") {
      if (is.null(ylim)) { ylim <- c(0,exp(max(lcd$phi))) }
      plot(xi, exp(eta), xlim=xlim, ylim=ylim, type="l", col=4, lwd=2,
           main=paste("density", ifelse(is.null(lcd$cure),"",paste(" (cure param. = ", lcd$cure, ")", sep="")), sep=""),
           xlab="x", ylab="f", ...)
      lines(c(xi[1]-0.04*diff(xlim), xi[1]), c(0,0), col=4, lwd=2, ...)
      if (sloperange & !is.na(lcd$cure.range[1])) {
        extraxi <- subdivisor(1, tau=c(tplus[m], xlim[2]+0.04*diff(xlim)), eps=stp)
        extraxi <- c(tplus[m], extraxi, xlim[2]+0.04*diff(xlim))
        extraeta1 <- phidivisor(1, tau=c(tplus[m], xlim[2]+0.04*diff(xlim)), phi=c(lcd$phi[m], lcd$phi[m]+phislrlow*(xlim[2]+0.04*diff(xlim)-tplus[m])), eps=stp)
        extraeta1 <- c(lcd$phi[m], extraeta1, lcd$phi[m]+phislrlow*(xlim[2]+0.04*diff(xlim)-tplus[m]))
        extraeta2 <- phidivisor(1, tau=c(tplus[m], xlim[2]+0.04*diff(xlim)), phi=c(lcd$phi[m], lcd$phi[m]+phislrhigh*(xlim[2]+0.04*diff(xlim)-tplus[m])), eps=stp)
        extraeta2 <- c(lcd$phi[m], extraeta2, lcd$phi[m]+phislrhigh*(xlim[2]+0.04*diff(xlim)-tplus[m]))
        extraxi <- c(extraxi,rev(extraxi))
        extraeta <- c(extraeta1,rev(extraeta2))
        polygon(extraxi, exp(extraeta), border=NA, col=grey(0.9), lwd=1, ...)
        box()      	
      }
      if (lcd$phislr > -Inf) {
        extraxi <- subdivisor(1, tau=c(xiknots[n], xlim[2]+0.04*diff(xlim)), eps=stp)
        extraxi <- c(xiknots[n], extraxi, xlim[2]+0.04*diff(xlim))
        extraeta <- phidivisor(1, tau=c(xiknots[n], xlim[2]+0.04*diff(xlim)), phi=c(etaknots[n], etaknots[n]+lcd$phislr*(xlim[2]+0.04*diff(xlim)-xiknots[n])), eps=stp)
        extraeta <- c(etaknots[n], extraeta, etaknots[n]+lcd$phislr*(xlim[2]+0.04*diff(xlim)-xiknots[n]))
        lines(extraxi, exp(extraeta), col=4, lwd=2, ...)
      } else {
        lines(c(xiknots[n], xlim[2]+0.04*diff(xlim)), c(0,0), col=4, lwd=2, ...)
      }
      if (kinkpoints) {
        points(xiknots, exp(etaknots), pch=19, ...)
      }
      if (kinklines) {
        abline(v=xiknots, lty=3, col=4, ...)
      }
      rug(lcd$tplus, ticksize=0.02)
      rug(lcd$tau, ticksize=0.04, lwd=1)
    }
    
    if (type == "CDF" | type == "survival") {
      if (is.null(ylim)) { ylim <- c(0,1) }
      subint <- diff(xi) * J00(eta[-n],eta[-1])
      Fhat <- c(0,cumsum(subint))
      if (sloperange & !is.na(lcd$cure.range[1])) {
        extraxi <- subdivisor(1, tau=c(tplus[m], xlim[2]+0.04*diff(xlim)), eps=stp)
        extraxi <- c(tplus[m], extraxi, xlim[2]+0.04*diff(xlim))	
        extraeta1 <- phidivisor(1, tau=c(tplus[m], xlim[2]+0.04*diff(xlim)), phi=c(lcd$phi[m],
          lcd$phi[m]+phislrlow*(xlim[2]+0.04*diff(xlim)-tplus[m])), eps=stp)
        extraeta1 <- c(lcd$phi[m], extraeta1, lcd$phi[m]+phislrlow*(xlim[2]+0.04*diff(xlim)-tplus[m]))
        extrasubint1 <- diff(extraxi) * J00(extraeta1[-n],extraeta1[-1])
        extraFhat1 <- tail(lcd$Fhat,1) + c(0,cumsum(extrasubint1))
        extraeta2 <- phidivisor(1, tau=c(tplus[m], xlim[2]+0.04*diff(xlim)), phi=c(lcd$phi[m],
          lcd$phi[m]+phislrhigh*(xlim[2]+0.04*diff(xlim)-tplus[m])), eps=stp)
        extraeta2 <- c(lcd$phi[m], extraeta2, lcd$phi[m]+phislrhigh*(xlim[2]+0.04*diff(xlim)-tplus[m]))
        extrasubint2 <- diff(extraxi) * J00(extraeta2[-n],extraeta2[-1])
        extraFhat2 <- tail(lcd$Fhat,1) + c(0,cumsum(extrasubint2)) 

        extraxigrey <- c(extraxi,rev(extraxi))
        extraFhatgrey <- c(extraFhat1,rev(extraFhat2))   
      }
      if (lcd$phislr > -Inf) {
        extraxi <- subdivisor(1, tau=c(xiknots[n], xlim[2]+0.04*diff(xlim)), eps=stp)
        extraxi <- c(xiknots[n], extraxi, xlim[2]+0.04*diff(xlim))
        extraeta <- phidivisor(1, tau=c(xiknots[n], xlim[2]+0.04*diff(xlim)), phi=c(etaknots[n],
          etaknots[n]+lcd$phislr*(xlim[2]+0.04*diff(xlim)-xiknots[n])), eps=stp)
        extraeta <- c(etaknots[n], extraeta, etaknots[n]+lcd$phislr*(xlim[2]+0.04*diff(xlim)-xiknots[n]))
        extrasubint <- diff(extraxi) * J00(extraeta[-n],extraeta[-1])
        extraFhat <- tail(Fhat,1) + c(0,cumsum(extrasubint))
        #print(extraxi)
        #print(extraeta)
        #print(extrasubint)
        #print(extraFhat)
      }
      #
      if (type == "CDF") {
        plot(xi, Fhat, xlim=xlim, ylim=ylim, type="l", col=4, lwd=2,
             main=paste("Fhat", ifelse(is.null(lcd$cure),"",paste(" (cure param. = ", lcd$cure, ")", sep="")), sep=""),
             xlab="x", ylab="F", ...)
        lines(c(xi[1]-0.04*diff(xlim), xi[1]), c(0,0), col=4, lwd=2, ...)
        if (sloperange & !is.na(lcd$cure.range[1])) {      
          polygon(extraxigrey, extraFhatgrey, border=NA, col=grey(0.9), lwd=1, ...)
          box()      	      
        }
        if (lcd$phislr > -Inf) {
          lines(extraxi, extraFhat, col=4, lwd=2, ...)
        } else {
          lines(c(xiknots[n], xlim[2]+0.04*diff(xlim)), c(1-lcd$cure,1-lcd$cure), col=4, lwd=2, ...)  
        }
        abline(h=c(0,1-lcd$cure), lty=2, col="gray70")
        if (kinkpoints) {
          points(xiknots, Fhat[knotlist], pch=19, ...)
        }
        if (kinklines) {
          abline(v=xiknots, lty=3, col=4, ...)
        }
        rug(lcd$tplus, ticksize=0.02)
        rug(lcd$tau, ticksize=0.04, lwd=1)
      }
      #
      if (type == "survival") {
        plot(xi, 1-Fhat, xlim=xlim, ylim=ylim, type="l", col=4, lwd=2,
             main=paste("Shat", ifelse(is.null(lcd$cure),"",paste(" (cure param. = ", lcd$cure, ")", sep="")), sep=""),
             xlab="x", ylab="S", ...)
        lines(c(xi[1]-0.04*diff(xlim), xi[1]), c(1,1), col=4, lwd=2, ...)
        if (sloperange & !is.na(lcd$cure.range[1])) {      
          polygon(extraxigrey, 1-extraFhatgrey, border=NA, col=grey(0.9), lwd=1, ...)
          box()      	      
        }
        if (lcd$phislr > -Inf) {
          lines(extraxi, 1-extraFhat, col=4, lwd=2, ...)
        } else {
          lines(c(xiknots[n], xlim[2]+0.04*diff(xlim)), c(lcd$cure,lcd$cure), col=4, lwd=2, ...)  
        }
        abline(h=c(lcd$cure,1), lty=2, col="gray70")
        if (kinkpoints) {
          points(xiknots, 1-Fhat[knotlist], pch=19, ...)
        }
        if (kinklines) {
          abline(v=xiknots, lty=3, col=4, ...)
        }
        rug(lcd$tplus, ticksize=0.02)
        rug(lcd$tau, ticksize=0.04, lwd=1)
      }
    }
  }
  invisible()
}

# print method for lcdensity class
print.lcdensity <- function(x, ...) {
  lcd <- x
  if (lcd$basedOn == "exact") {
    message("Pretty print not implemented yet")
    return(print.default(lcd))
  }
  mm <- length(lcd$tplus)
  m <- ifelse(lcd$phislr > -Inf, mm-1, mm)
  n <- dim(lcd$x)[1]
  if (n <= 20 && mm <= 50) {
    print.default(lcd)
  } else {
    cat("\nLog-concave maximum likelihood estimate for", lcd$basedOn, "data\n")
    cat("  based on", n, "data intervals/points\n\n")
    cat("Exit status:", lcd$status, "\n")
    cat("Domain: [", format(lcd$tau[lcd$domind1]), ", ", format(lcd$tau[lcd$domind2]),
        ifelse(lcd$tau[lcd$domind2] == Inf,")","]"), "\n", sep="")
    cat("In total ", sum(lcd$isKnot), " knots",
        ifelse(lcd$phislr > -Inf, paste(" and a right-hand slope of", format(lcd$phislr), "\n"), ", no slope \n"), sep="")
    if (!is.null(lcd$cure)) {
      cat("Cure parameter p0 = ", lcd$cure, "\n", sep="")
    }
#    if (lcd$force.inf) {
#      cat("Right-hand slope was enforced.\n")
#    }
  }
  invisible(lcd)
}


# summary method for lcdensity class
summary.lcdensity <- function(object, ...) {
  lcd <- object
  if (lcd$basedOn == "exact") {
    message("pretty summary not implemented yet")
    return(summary.default(lcd))
  }
  mm <- length(lcd$tplus)
  m <- ifelse(lcd$phislr > -Inf, mm-1, mm)
  n <- dim(lcd$x)[1]
  basedOn <- lcd$basedOn
  status <- lcd$status
  domind <- c(lcd$domind1, lcd$domind2)
  dom <- lcd$tau[domind]
  knots <- lcd$tplus[1:m][as.logical(lcd$isKnot)]
  phiatknots <- lcd$phi[as.logical(lcd$isKnot)]
  phislr <- lcd$phislr
  Fhatatknots <- lcd$Fhat[as.logical(lcd$isKnot)]
  if (is.na(lcd$cure.range[1])) {
    return(list(basedOn=basedOn, status=status, domind=domind, dom=dom, knots=knots, phiatknots=phiatknots,
              phislr=phislr, cure=lcd$cure, Fhatatknots=Fhatatknots))
  } else {
    return(list(basedOn=basedOn, status=status, domind=domind, dom=dom, knots=knots, phiatknots=phiatknots,
              phislr=phislr, phislr.range=lcd$phislr.range, cure=lcd$cure, cure.range=lcd$cure.range,
              Fhatatknots=Fhatatknots))
  }
}
