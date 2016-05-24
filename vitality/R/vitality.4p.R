## Four parameter simple model r s lambda beta model, no childhood hook

#' Fitting routine for the 2-process, 4-parameter vitality model (no childhood hook). 
#' 
#' Based on code by D.H. Salinger, J.J. Anderson and O. Hamel (2003).
#' "A parameter fitting routine for the vitality based survival model."
#' Ecological Modeling 166(3): 287--294.
#' 
#' @param time Vector. Time component of data: Defaults to \code{0:(1-length(sdata))}.
#' @param sdata Required. Survival or mortality data.  The default expects cumulative 
#'        survival fraction.  If providing incremental mortality fraction 
#'        instead, use option: datatype = "INC". 
#'        The default also expects the data to represent full mortality. 
#'        Otherwise, use option: rc.data = T to indicate right censored data.
#' @param rc.data Optional, boolean. Specifies Right Censored data.   If the data does not 
#'        represent full mortality, it is probably right censored.  The default 
#'        is rc.data = F.  A third option is rc.data = "TF".  Use this case to add
#'        a near-term zero survival point to data which displays nearly full 
#'        mortality ( <.01 survival at end).  If rc.data = F but the data does
#'        not show full mortality, rc.data = "TF" will be 
#'        invoked automatically. 
#' @param se Optional, boolean. Calculates the standard errors for the MLE parameters.
#'        Default is FALSE. Set equal to the initial study population to 
#'        compute standard errors. 
#' @param datatype Optional. Defaults to \code{"CUM"} for cumulative survival fraction data.
#'        Use \code{"INC"} - for incremental mortality fraction data. 
#' @param ttol Optional. Stopping criteria tolerance.  Default is 1e-6.
#'        Specify as ttol = .0001. If one of the liklihood plots (esp. for "k") does not look optimal, 
#'        try decreasing ttol.   If the program crashes, try increasing ttol.
#' @param init.params Optional. Please specify the initial param values.
#'        specify \code{init.params = c(r, s, lambda, beta)} in that order
#'        (eg. init.params = c(.1, .02, .3, 0.12)).
#' @param pplot Optional, boolean. Plots of cumulative survival for both data and fitted curves?
#'        Default \code{TRUE}. \code{FALSE} Produce no plots. A
#'        A third option:  \code{pplot = n}  (n >= 1) extends the time axis of 
#'        the fitting plots (beyond the max time in data).  For example: 
#'        \code{pplot = 1.2} extends the time axis by 20%. Note:  the incremental 
#'        mortality plot is a continuous representation of the appropriately-
#'        binned histogram of incremental mortalities.
#' @param Iplot Optional, boolean. Incremental mortality for both data and fitted curves?
#'        Default: \code{FALSE}.
#' @param Mplot Optional, boolean. Plot fitted mortality curve? Default is \code{FALSE}.
#' @param tlab Optional, character. specifies units for x-axis of plots.  Default is "days".
#' @param silent Optional, boolean. Stops all print and plot options (still get most warning and all 
#'          error messages) Default is \code{FALSE}.  A third option, \code{"verbose"} also 
#'          enables the trace setting in the ms (minimum sum) S-Plus routine.
#' @export 
#' @return vector of final MLE r, s, lambda, beta w/wo gamma and alpha parameter estimates.
#'     standard errors of MLE parameter estimates (if se = <population> is specified).
vitality.4p <- function(time = 0:(length(sdata)-1), 
                         sdata, 
                         init.params = FALSE,
                         lower = c(0, 0, 0, 0),
                         upper = c(100,50,100,50),
                         rc.data = FALSE, 
                         se = FALSE, 
                         datatype = c("CUM", "INC"), 
                         ttol = 1e-6, 
                         pplot = TRUE, 
                         Iplot = FALSE, 
                         Mplot = FALSE, 
                         tlab = "years", 
                         silent = FALSE) {
    #  --Check/prepare Data---
    datatype <- match.arg(datatype)
    if (length(time) != length(sdata)) {
        stop("time and sdata must have the same length")
    }
    in.time <- time
    dTmp <- dataPrep(time, sdata, datatype, rc.data)
    time <- dTmp$time
    sfract <- dTmp$sfract
    x1 <- dTmp$x1
    x2 <- dTmp$x2
    Ni <- dTmp$Ni
    rc.data <- dTmp$rc.data
    if(in.time[1]>0){
      time <- time[-1]
      sfract <- sfract[-1]
      x1 <- c(x1[-c(1,length(x1))], x1[1])
      x2 <- c(x2[-c(1,length(x2))], 0)
      Ni <- Ni[-1]
      rc.data <- rc.data[-1]
    }
    
    
    #  --Produce initial parameter values---
    if(length(init.params) == 1) {  
        ii <- indexFinder(sfract, 0.5)
        if (ii == -1) {
            warning("ERROR: no survival fraction data below the .5 level.\n
                    Cannot use the initial r s k u estimator.  You must supply initial r s k u estimates")
            return(-1)
        }
        else rsk <- c(1/time[ii], 0.01, 0.1, 0.1)  
    } else {                      # use user specified init params
        rsk <- init.params
    }
    if (rsk[1] == -1) {
        stop
    }
    if (silent == FALSE) {
        print(cbind(c("Initial r", "initial s", "initial lambda", "initial beta"), rsk))
    }
    
    #  --create dataframe for sa---
    dtfm <- data.frame(x1 = x1, x2 = x2, Ni = Ni)
    
    
    #  --run MLE fitting routine---   
    # --conduct Newton-Ralphoson algorithm directly --
    fit.nlm <- nlminb(start = rsk, objective = logLikelihood.4p, lower = lower,
                      upper = upper, xx1 = x1, xx2 = x2, NNi = Ni)
    
    # --save final param estimates---
    r.final <- fit.nlm$par[1]
    s.final <- abs(fit.nlm$par[2])
    lambda.final <- fit.nlm$par[3]
    beta.final <- fit.nlm$par[4]
    mlv <- fit.nlm$obj  
    if (silent == FALSE) {print(cbind(c("estimated r", "estimated s", "estimated lambda",
                                        "estimated beta", "minimum -loglikelihood value"),
                                      c(r.final, s.final, lambda.final, beta.final, mlv)))}
    
    #  == end MLE fitting == = 
    
    #  --compute standard errors---
    if (se != FALSE) {
        s.e. <- stdErr.4p(r.final, s.final, lambda.final, beta.final, x1, x2, Ni, se)  
        if (silent == FALSE){print(cbind(c("sd for r", "sd for s", "sd for lambda", "sd for beta"), s.e.))}  
        
    }
    
    #  --plotting and goodness of fit---
    if (pplot != FALSE) {
        plotting.4p(r.final, s.final, lambda.final, beta.final, mlv, time, sfract, x1, x2, Ni, pplot, Iplot, Mplot, tlab, rc.data)
    }
    # ............................................................................................
    
    # --return final param values---
    sd <- 5   #significant digits of output
    if(se != F ) {
        params <- c(r.final, s.final, lambda.final, beta.final)
        pvalue <- c(1-pnorm(r.final/s.e.[1]), 1-pnorm(s.final/s.e.[2]), 1-pnorm(lambda.final/s.e.[3]), 1-pnorm(beta.final/s.e.[4]))
        std <- c(s.e.[1], s.e.[2], s.e.[3], s.e.[4])
        out <- signif(cbind(params, std, pvalue), sd)
        return(out)
    }
    else {
        return(signif(c(r.final, s.final, lambda.final, beta.final), sd))
    }
}


#' Plotting function for 2-process vitality model. 4-param
#' 
#' None.
#' 
#' @param r.final r estimate
#' @param s.final s estimate
#' @param lambda.final lambda estimate
#' @param beta.final beta estimate
#' @param mlv TODO mlv
#' @param time time vector
#' @param sfract survival fraction
#' @param x1 Time 1
#' @param x2 Time 2
#' @param Ni Initial population
#' @param pplot Boolean. Plot cumulative survival fraction?
#' @param Iplot Boolean. Plot incremental survival?
#' @param Mplot Boolean. Plot mortality rate? 
#' @param tlab Character, label for time axis
#' @param rc.data Booolean, right-censored data?
plotting.4p <- function(r.final,
                         s.final,
                         lambda.final,
                         beta.final,
                         mlv,
                         time,
                         sfract,
                         x1, 
                         x2, 
                         Ni, 
                         pplot, 
                         Iplot, 
                         Mplot, 
                         tlab, 
                         rc.data) {
    # --plot cumulative survival---
    if (pplot != FALSE) {
        #win.graph()
        ext <- max(pplot, 1)
        par(mfrow = c(1, 1))
        len <- length(time)
        tmax <- ext * time[len]
        plot(time, sfract/sfract[1], xlab = tlab, ylab = "survival fraction",
             ylim = c(0, 1), xlim = c(min(time), tmax), col = 1)
        xxx <- seq(min(time), tmax, length = 200)
        xxx1 <- c(0, xxx[-1])
        lines(xxx, SurvFn.4p(xxx1, r.final, s.final, lambda.final, beta.final), col = 2, lwd=2)
        lines(xxx, SurvFn.in.4p(xx=xxx1, r=r.final, s=s.final), col=3, lwd=2, lty=3)
        lines(xxx, SurvFn.ex.4p(xx=xxx1, r=r.final, s=s.final, lambda = lambda.final, beta = beta.final), col=4, lwd=2, lty=2)
        title("Cumulative Survival Data and Vitality Model Fitting")
        legend(x="bottomleft", legend=c("Total", "Intrinsic", "Extrinsic"), lty=c(1, 3, 2), col=c(2,3,4), bty="n", lwd=c(2,2,2))
    } 
    
    if ( Mplot != FALSE) {
      lx <- round(sfract*100000)
      lx <- c(lx,0)
      ndx <- -diff(lx)
      lxpn <- lx[-1]
      n <- c(diff(time), 1000)
      nax <- .5*n
      nLx <- n * lxpn + ndx * nax
      mu.x <- ndx/nLx
      mu.x[length(mu.x)] <- NA
        # qx <- Ni/sfract
        # mu.x <- 2 * qx/(2 - qx)
        #win.graph()
        ext <- max(pplot, 1)
        par(mfrow = c(1, 1))
        len <- length(time)
        tmax <- ext * time[len]
        xxx <- seq(min(time), tmax, length = 200)
        xxx1 <- xxx
        mu.i <- mu.vd1.4p(xxx1, r.final, s.final)
        mu.e <- mu.vd2.4p(xxx1, r.final, lambda.final, beta.final)
        mu.t <- mu.vd.4p(xxx1, r.final, s.final, lambda.final, beta.final)
        plot(time, mu.x, xlim = c(time[1], tmax), xlab = tlab, ylab = "estimated mortality rate", log = "y",
             main = "Log Mortality Data and Vitality Model Fitting", ylim=c(min(mu.x,mu.t,na.rm=T),max(mu.x,mu.t,na.rm=T)))
        #plot(xxx, mu.t, xlim = c(0, tmax), xlab = tlab, ylab = "estimated mortality rate", log = "y",
        #    type = "l")
        lines(xxx, mu.t, col = 2, lwd=2)
        lines(xxx, mu.i, col=3, lwd=2, lty=3)
        lines(xxx, mu.e, col=4, lwd=2, lty=2)
        legend(x="bottomright", legend=c("data (approximate)", expression(mu[total]), expression(mu[i]), expression(mu[e])), lty=c(NA, 1, 3, 2), pch=c(1,NA,NA,NA), col=c(1,2,3,4), bty="n", lwd=c(1,2,2,2))
    } 
    
    # --Incremental mortality plot
    if (Iplot != FALSE) {
        #win.graph()
        par(mfrow = c(1, 1))
        
        ln <- length(Ni)-1
        x1 <- x1[1:ln]
        x2 <- x2[1:ln]
        Ni <- Ni[1:ln]
        
        ln <- length(Ni)
        #scale <- (x2-x1)[Ni == max(Ni)]
        scale <- max( (x2-x1)[Ni == max(Ni)] )
        
        ext <- max(pplot, 1)
        
        npt <- 200*ext
        xxx <- seq(x1[1], x2[ln]*ext, length = npt)
        xx1 <- xxx[1:(npt-1)]
        xx2 <- xxx[2:npt]
        sProbI <- survProbInc.4p(r.final, s.final, lambda.final, beta.final, xx1, xx2)
        
        ytop <- 1.1 * max(max(sProbI/(xx2-xx1)), Ni/(x2-x1)) * scale		
        plot((x1+x2)/2, Ni*scale/(x2-x1), ylim = c(0, ytop), xlim = c(x1[1], ext*x2[ln]),
             xlab = tlab, ylab = "incremental mortality")
        title("Probability Density Function")
        lines((xx1+xx2)/2, sProbI*scale/(xx2-xx1), col=2)
    }
    
    #return()	 
}

#' The intrinsic cumulative survival distribution function for 2-process 4-parameter
#' 
#' None.
#' 
#' @param xx vector of ages
#' @param r r value
#' @param s s value
#' @return vector of FF?
SurvFn.in.4p <- function(xx, r, s) {
  yy <- s^2*xx
  # pnorm is: cumulative prob for the Normal Dist.
  tmp1 <- sqrt(1/yy) * (1 - xx * r)    #  xx = 0 is ok.  pnorm(+-Inf) is defined
  tmp2 <- sqrt(1/yy) * (1 + xx * r)
  
  # --safeguard if exponent gets too large.---
  tmp3 <- 2*r/(s*s)
  
  if (tmp3 > 250) {   
    q <- tmp3/250 
    
    if (tmp3 > 1500) {
      q <- tmp3/500
    }
    
    valueFF <- (1.-(pnorm(-tmp1) + (exp(tmp3/q) *pnorm(-tmp2)^(1/q))^(q)))#*exp(-lambda*exp(-1/beta)/(1/beta*r)*(exp(1/beta*r*xx)-1))
    #valueFF <- (1.-(pnorm(-tmp1) + (exp(tmp3/q) *pnorm(-tmp2)^(1/q))^(q)))*exp(-a/b*(exp(b*xx)-1))  
  }		    
  else {
    valueFF <- (1.-(pnorm(-tmp1) + exp(tmp3) *pnorm(-tmp2)))#*exp(-lambda*exp(-1/beta)/(1/beta*r)*(exp(1/beta*r*xx)-1))   #1-G
    #valueFF <- (1.-(pnorm(-tmp1) + exp(tmp3) *pnorm(-tmp2)))*exp(-a/b*(exp(b*xx)-1))   #1-G
    
  }
  if ( all(is.infinite(valueFF)) ) {
    warning(message = "Inelegant exit caused by overflow in evaluation of survival function.  
            Check for right-censored data. Try other initial values.")
  }
  
  return(valueFF)	
  }

#' The cumulative survival distribution function for 2-process 4-parameter
#' 
#' None.
#' 
#' @param xx vector of ages
#' @param r r value
#' @param s s value
#' @param lambda lambda value
#' @param beta beta value
#' @return vector of FF?
SurvFn.ex.4p <- function(xx, r, s, lambda, beta) {
    yy <- s^2*xx
    # pnorm is: cumulative prob for the Normal Dist.
    tmp1 <- sqrt(1/yy) * (1 - xx * r)    #  xx = 0 is ok.  pnorm(+-Inf) is defined
    tmp2 <- sqrt(1/yy) * (1 + xx * r)
    
    # --safeguard if exponent gets too large.---
    tmp3 <- 2*r/(s*s)
    
    if (tmp3 > 250) {   
        q <- tmp3/250 
        
        if (tmp3 > 1500) {
            q <- tmp3/500
        }
        
        #valueFF <- exp(-lambda*exp(-1/beta)/(r/beta)*(exp(r*xx/beta)-1) +gamma/alpha*(exp(-alpha*xx)-1))
        valueFF <- exp(-lambda*exp(-1/beta)/(r/beta)*(exp(r*xx/beta)-1))
    } else {
      #valueFF <-exp(-lambda*exp(-1/beta)/(r/beta)*(exp(r*xx/beta)-1)+ gamma/alpha*(exp(-alpha*xx)-1))
      valueFF <-exp(-lambda*exp(-1/beta)/(r/beta)*(exp(r*xx/beta)-1)) 
        
    }
    if ( all(is.infinite(valueFF)) ) {
        warning(message = "Inelegant exit caused by overflow in evaluation of survival function.  
		Check for right-censored data. Try other initial values.")
	}
    
    return(valueFF)	
}


#' Calculates incremental survival probability for 2-process 4-parameter r, s, lambda, beta
#' 
#' None
#' 
#' @param r r value
#' @param s s value
#' @param lambda lambda value
#' @param beta beta value
#' @param xx1 xx1 vector
#' @param xx2 xx2 vector
#' @return Incremental survival probabilities.
survProbInc.4p <- function(r, s, lambda, beta, xx1, xx2){
    value.iSP <- -(SurvFn.4p(xx2, r, s, lambda, beta) - SurvFn.4p(xx1, r, s, lambda, beta))
    value.iSP[value.iSP < 1e-18] <- 1e-18   # safeguards against taking Log(0)
    value.iSP
}


#' Gives log likelihood of 2-process 4-parameter model
#' 
#' None
#' 
#' @param par vector of parameter(r, s, lambda, beta)
#' @param xx1 xx1 vector
#' @param xx2 xx2 vector
#' @param NNi survival fractions
#' @return log likelihood
logLikelihood.4p <- function(par, xx1, xx2, NNi) {
    # --calculate incremental survival probability--- (safeguraded > 1e-18 to prevent log(0))
    iSP <-  survProbInc.4p(par[1], par[2], par[3], par[4], xx1, xx2)
    loglklhd <- -NNi*log(iSP)
    return(sum(loglklhd))
}


#' Standard errors for 4-param r, s, lambda, beta
#' 
#' Note: if k <= 0, can not find std Err for k.
#' 
#' @param r r value
#' @param s s value
#' @param k k value
#' @param u u value corresponding to beta?
#' @param x1 age 1 (corresponding 1:(t-1) and 2:t)
#' @param x2 age 2
#' @param Ni survival fraction
#' @param pop initial population (total population of the study)
#' @return standard error for r, s, k, u.
stdErr.4p <- function(r, s, k, u, x1, x2, Ni, pop) {	
    LL <- function(a, b, c, d, r, s, k, u, x1, x2, Ni) {logLikelihood.4p(c(r+a, s+b, k+c, u+d), x1, x2, Ni)}
    
    #initialize hessian for storage
    hess <- matrix(0, nrow = 4, ncol = 4)
    
    #set finite difference intervals
    h <- .001
    hr <- abs(h*r)
    hs <- h*s*.1
    hk <- h*k*.1
    hu <- h*u*.1
    
    #Compute second derivitives (using 5 point)
    # LLrr
    f0 <- LL(-2*hr, 0, 0, 0, r, s, k, u, x1, x2, Ni)
    f1 <- LL(-hr, 0, 0, 0, r, s, k, u, x1, x2, Ni)
    f2 <- LL(0, 0, 0, 0, r, s, k, u, x1, x2, Ni)
    f3 <- LL(hr, 0, 0, 0, r, s, k, u, x1, x2, Ni)
    f4 <- LL(2*hr, 0, 0, 0, r, s, k, u, x1, x2, Ni)
    
    fp0 <- (-25*f0 +48*f1 -36*f2 +16*f3 -3*f4)/(12*hr)
    fp1 <- (-3*f0 -10*f1 +18*f2 -6*f3 +f4)/(12*hr)
    fp3 <- (-f0 +6*f1 -18*f2 +10*f3 +3*f4)/(12*hr)
    fp4 <- (3*f0 -16*f1 +36*f2 -48*f3 +25*f4)/(12*hr)
    
    LLrr <- (fp0 -8*fp1 +8*fp3 -fp4)/(12*hr)
    
    # LLss
    f0 <- LL(0, -2*hs, 0, 0, r, s, k, u, x1, x2, Ni)
    f1 <- LL(0, -hs, 0, 0, r, s, k, u, x1, x2, Ni)
    # f2 as above
    f3 <- LL(0, hs, 0, 0, r, s, k, u, x1, x2, Ni)
    f4 <- LL(0, 2*hs, 0, 0, r, s, k, u, x1, x2, Ni)
    
    fp0 <- (-25*f0 +48*f1 -36*f2 +16*f3 -3*f4)/(12*hs)
    fp1 <- (-3*f0 -10*f1 +18*f2 -6*f3 +f4)/(12*hs)
    fp3 <- (-f0 +6*f1 -18*f2 +10*f3 +3*f4)/(12*hs)
    fp4 <- (3*f0 -16*f1 +36*f2 -48*f3 +25*f4)/(12*hs)
    
    LLss <- (fp0 -8*fp1 +8*fp3 -fp4)/(12*hs)
    
    # LLkk
    
    f0 <- LL(0, 0, -2*hk, 0, r, s, k, u, x1, x2, Ni)
    f1 <- LL(0, 0, -hk, 0, r, s, k, u, x1, x2, Ni)
    # f2 as above
    f3 <- LL(0, 0, hk, 0, r, s, k, u, x1, x2, Ni)
    f4 <- LL(0, 0, 2*hk, 0, r, s, k, u, x1, x2, Ni)
    
    fp0 <- (-25*f0 +48*f1 -36*f2 +16*f3 -3*f4)/(12*hk)
    fp1 <- (-3*f0 -10*f1 +18*f2 -6*f3 +f4)/(12*hk)
    fp3 <- (-f0 +6*f1 -18*f2 +10*f3 +3*f4)/(12*hk)
    fp4 <- (3*f0 -16*f1 +36*f2 -48*f3 +25*f4)/(12*hk)
    
    LLkk <- (fp0 -8*fp1 +8*fp3 -fp4)/(12*hk)
    
    
    # LLuu
    f0 <- LL(0, 0, 0, -2*hu, r, s, k, u, x1, x2, Ni)
    f1 <- LL(0, 0, 0, -hu, r, s, k, u, x1, x2, Ni)
    # f2 as above
    f3 <- LL(0, 0, 0, hu, r, s, k, u, x1, x2, Ni)
    f4 <- LL(0, 0, 0, 2*hu, r, s, k, u, x1, x2, Ni)
    
    fp0 <- (-25*f0 +48*f1 -36*f2 +16*f3 -3*f4)/(12*hu)
    fp1 <- (-3*f0 -10*f1 +18*f2 -6*f3 +f4)/(12*hu)
    fp3 <- (-f0 +6*f1 -18*f2 +10*f3 +3*f4)/(12*hu)
    fp4 <- (3*f0 -16*f1 +36*f2 -48*f3 +25*f4)/(12*hu)
    
    LLuu <- (fp0 -8*fp1 +8*fp3 -fp4)/(12*hu)
    
    #-------end second derivs---
    # do mixed partials (4 points)
    # LLrs
    m1 <- LL(hr, hs, 0, 0, r, s, k, u, x1, x2, Ni)
    m2 <- LL(-hr, hs, 0, 0, r, s, k, u, x1, x2, Ni)
    m3 <- LL(-hr, -hs, 0, 0, r, s, k, u, x1, x2, Ni)
    m4 <- LL(hr, -hs, 0, 0, r, s, k, u, x1, x2, Ni)
    
    LLrs <- (m1 -m2 +m3 -m4)/(4*hr*hs)
    
    # LLru
    m1 <- LL(hr, 0, 0, hu, r, s, k, u, x1, x2, Ni)
    m2 <- LL(-hr, 0, 0, hu, r, s, k, u, x1, x2, Ni)
    m3 <- LL(-hr, 0, 0, -hu, r, s, k, u, x1, x2, Ni)
    m4 <- LL(hr, 0, 0, -hu, r, s, k, u, x1, x2, Ni)
    
    LLru <- (m1 -m2 +m3 -m4)/(4*hr*hu)
    
    # LLsu
    m1 <- LL(0, hs, 0, hu, r, s, k, u, x1, x2, Ni)
    m2 <- LL(0, -hs, 0, hu, r, s, k, u, x1, x2, Ni)
    m3 <- LL(0, -hs, 0, -hu, r, s, k, u, x1, x2, Ni)
    m4 <- LL(0, hs, 0, -hu, r, s, k, u, x1, x2, Ni)
    
    LLsu <- (m1 -m2 +m3 -m4)/(4*hu*hs)
    
    
    # LLrk
    m1 <- LL(hr, 0, hk, 0, r, s, k, u, x1, x2, Ni)
    m2 <- LL(-hr, 0, hk, 0, r, s, k, u, x1, x2, Ni)
    m3 <- LL(-hr, 0, -hk, 0, r, s, k, u, x1, x2, Ni)
    m4 <- LL(hr, 0, -hk, 0, r, s, k, u, x1, x2, Ni)
    
    LLrk <- (m1 -m2 +m3 -m4)/(4*hr*hk)
    
    # LLsk
    m1 <- LL(0, hs, hk, 0, r, s, k, u, x1, x2, Ni)
    m2 <- LL(0, -hs, hk, 0, r, s, k, u, x1, x2, Ni)
    m3 <- LL(0, -hs, -hk, 0, r, s, k, u, x1, x2, Ni)
    m4 <- LL(0, hs, -hk, 0, r, s, k, u, x1, x2, Ni)
    
    LLsk <- (m1 -m2 +m3 -m4)/(4*hs*hk)
    
    # LLku
    m1 <- LL(0, 0, hk, hu, r, s, k, u, x1, x2, Ni)
    m2 <- LL(0, 0, hk, -hu, r, s, k, u, x1, x2, Ni)
    m3 <- LL(0, 0, -hk, -hu, r, s, k, u, x1, x2, Ni)
    m4 <- LL(0, 0, -hk, hu, r, s, k, u, x1, x2, Ni)
    
    LLku <- (m1 -m2 +m3 -m4)/(4*hu*hk)
    
    
    
    diag(hess) <- c(LLrr, LLss, LLkk, LLuu)*pop
    hess[2, 1] = hess[1, 2] <- LLrs*pop
    hess[3, 1] = hess[1, 3] <- LLrk*pop
    hess[3, 2] = hess[2, 3] <- LLsk*pop
    hess[4, 1] = hess[1, 4] <- LLru*pop
    hess[4, 2] = hess[2, 4] <- LLsu*pop
    hess[4, 3] = hess[3, 4] <- LLku*pop
    
    
    #print(hess)
    hessInv <- solve(hess)
    #print(hessInv)
    
    #compute correlation matrix:
    sz <- 4
    corr <- matrix(0, nrow = sz, ncol = sz)
    for (i in 1:sz) {
        for (j in 1:sz) {
            corr[i, j] <- hessInv[i, j]/sqrt(abs(hessInv[i, i]*hessInv[j, j]))
        }
    }
    #print(corr)
    if ( abs(corr[2, 1]) > .98 ) {
        warning("WARNING: parameters r and s appear to be closely correlated for this data set.  
                s.e. may fail for these parameters.")
	}
    if (  sz == 4 && abs(corr[3, 2]) > .98 ) {
        warning("WARNING: parameters s and lambda appear to be closely correlated for this data set.  
                s.e. may fail for these parameters.")
	}
    if (  sz == 4 && abs(corr[3, 1]) > .98 ) {
        warning("WARNING: parameters r and lambda appear to be closely correlated for this data set.  
                s.e. may fail for these parameters.")
	}
    if (  sz == 4 && abs(corr[4, 2]) > .98 ) {
        warning("WARNING: parameters s and beta appear to be closely correlated for this data set.  
                s.e. may fail for these parameters.")
	}
    if (  sz == 4 && abs(corr[4, 1]) > .98 ) {
        warning("WARNING: parameters r and beta appear to be closely correlated for this data set.  
                s.e. may fail for these parameters.")
	}
    
    se <- sqrt(diag(hessInv))
    
    #  Approximate s.e. for cases where calculation of s.e. failed:
    if( sum( is.na(se) ) > 0 ) {
        seNA <- is.na(se)
        se12 <- sqrt(diag(solve(hess[c(1, 2)	, c(1, 2) ])))
        se13 <- sqrt(diag(solve(hess[c(1, 3)	, c(1, 3) ])))
        se23 <- sqrt(diag(solve(hess[c(2, 3)	, c(2, 3) ])))
        se14 <- sqrt(diag(solve(hess[c(1, 4)	, c(1, 4) ])))
        se24 <- sqrt(diag(solve(hess[c(2, 4)	, c(2, 4) ])))
        se34 <- sqrt(diag(solve(hess[c(3, 4)	, c(3, 4) ])))
        
        if(seNA[1]) {
            if(!is.na(se12[1]) ){
                se[1] = se12[1]
                warning("* s.e. for parameter r is approximate.")
            }
            else if(!is.na(se13[1])){
                se[1] = se13[1]
                warning("* s.e. for parameter r is approximate.")
            }   
            else if(!is.na(se14[1])){
                se[1] = se14[1]
                warning("* s.e. for parameter r is approximate.")
            }   
            else warning("* unable to calculate or approximate s.e. for parameter r.")
        }
        if(seNA[2]) {
            if(!is.na(se12[2]) ){
                se[2] = se12[2]
                warning("* s.e. for parameter s is approximate.")
            }
            else if(!is.na(se23[1])){
                se[2] = se23[1]
                warning("* s.e. for parameter s is approximate.")
            }   
            else if(!is.na(se24[1])){
                se[2] = se24[1]
                warning("* s.e. for parameter s is approximate.")
            }   
            else warning("* unable to calculate or approximate s.e. for parameter s.")
        }
        if(seNA[3]) {
            if(!is.na(se13[2]) ){
                se[3] = se13[2]
                warning("* s.e. for parameter lambda is approximate.")
            }
            else if(!is.na(se23[2])){
                se[3] = se23[2]
                warning("* s.e. for parameter lambda is approximate.")
            }   
            else if(!is.na(se34[1])){
                se[3] = se34[1]
                warning("* s.e. for parameter lambda is approximate.")
            }   
            else warning("* unable to calculate or approximate s.e. for parameter lambda.")
        }
        if(seNA[4]) {
            if(!is.na(se14[2]) ){
                se[4] = se14[2]
                warning("* s.e. for parameter beta is approximate.")
            }
            else if(!is.na(se24[2])){
                se[4] = se24[2]
                warning("* s.e. for parameter beta is approximate.")
            }   
            else if(!is.na(se34[1])){
                se[4] = se34[2]
                warning("* s.e. for parameter beta is approximate.")
            }   
            else warning("* unable to calculate or approximate s.e. for parameter beta.")
        }
        
    }
    
    
    #######################
    return(se)
}

SurvFn.4p <- function(xx,r,s,lambda,beta){
  yy <- s^2*xx
  # pnorm is: cumulative prob for the Normal Dist.
  tmp1 <- sqrt(1/yy) * (1 - xx * r)    #  xx = 0 is ok.  pnorm(+-Inf) is defined
  tmp2 <- sqrt(1/yy) * (1 + xx * r)
  
  # --safeguard if exponent gets too large.---
  tmp3 <- 2*r/(s*s)
  
  if (tmp3 > 250) {   
    q <- tmp3/250 
    
    if (tmp3 > 1500) {
      q <- tmp3/500
    }
    
    valueFF <-(1.-(pnorm(-tmp1) + (exp(tmp3/q) *pnorm(-tmp2)^(1/q))^(q)))*exp(-lambda*exp(-1/beta)/(r/beta)*(exp(r*xx/beta)-1)) # This requires 1/alpha
  } else {
    valueFF <-(1.-(pnorm(-tmp1) + exp(tmp3) *pnorm(-tmp2)))*exp(-lambda*exp(-1/beta)/(r/beta)*(exp(r*xx/beta)-1)) 
  }
  if ( all(is.infinite(valueFF)) ) {
    warning(message="Inelegant exit caused by overflow in evaluation of survival function. Check for right-censored data. Try other initial values.")
  }
  
  return(valueFF)  
}


#' Intrinsic cumulative survival distribution
#' 
#' None
#' 
#' @param xx vector of ages
#' @param r r value
#' @param s s value
#' @param u u value
#' @return Cumulative survival distribution
SurvFn.h.4p <- function(xx, r, s, u){
    yy <- u^2+s^2*xx
    # pnorm is: cumulative prob for the Normal Dist.
    tmp1 <- sqrt(1/yy) * (1 - xx * r)    #  xx = 0 is ok.  pnorm(+-Inf) is defined
    tmp2 <- sqrt(1/yy) * (1 + xx * r+2*u^2*r/s^2)
    
    # --safeguard if exponent gets too large.---
    tmp3 <- 2*r/(s*s)+2*u^2*r^2/s^4
    
    if (tmp3 > 250) {   
        q <- tmp3/250 
        
        if (tmp3 > 1500) {
            q <- tmp3/500
        }
        
        valueFF <- (1.-(pnorm(-tmp1) + (exp(tmp3/q) *pnorm(-tmp2)^(1/q))^(q))) 
    }		    
    else {
        valueFF <- (1.-(pnorm(-tmp1) + exp(tmp3) *pnorm(-tmp2)))  #1-G
    }
    if ( all(is.infinite(valueFF)) ) {
        warning(message = "Inelegant exit caused by overflow in evaluation of survival function.  
		Check for right-censored data. Try other initial values.")
	}
    
    return(valueFF)	
}
