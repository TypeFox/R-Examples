#  File R/stergm.EGMME.GD.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
stergm.EGMME.GD <- function(theta.form0, theta.diss0, nw, model.form, model.diss, model.mon,
                            control, MHproposal.form, MHproposal.diss, cl=NULL,
                            verbose=FALSE){

  
  ###### Set the constants and convenience variables. ######
  offsets <- c(model.form$etamap$offsettheta, model.diss$etamap$offsettheta) # which parameters are offsets?
  p.form.free <- sum(!model.form$etamap$offsettheta) # number of free formation parameters
  p.form <- length(model.form$etamap$offsettheta) # total number of formation parameters
  p.diss.free <- sum(!model.diss$etamap$offsettheta) # number of free dissolution parameters
  p.diss <- length(model.diss$etamap$offsettheta) # total number of dissolution parameters
  p.free <- p.form.free+p.diss.free  # number of free parameters (formation and dissolution)
  p <- p.form+p.diss # total number of parameters (free and offset)
  p.names<-c(paste("f.(",model.form$coef.names,")",sep=""),paste("d.(",model.diss$coef.names,")",sep=""))
  
  q <- length(model.mon$etamap$offsettheta) # number of target statistics
  q.names<-model.mon$coef.names

  # Define the function to set optimization parameters.  
  eval.optpars <- function(states, history, control, test.G,window,update.jitter){
    ## Regress statistics on parameters.
    # This uses GLS to account for serial correlation in statistics,
    # since we want p-values. First row is the intercept.

    oh <- if(window) history$oh else history$oh.all
    oh.last <- history$oh.last
    ind <- if(window) history$ind else history$ind.all
    tid <- if(window) history$tid else history$tid.all

    x<-oh[,1:p,drop=FALSE][,!offsets,drop=FALSE] # #$%^$ gls() doesn't respect I()...
    ys <- oh[,-(1:p),drop=FALSE]
    n <- nrow(ys)
    h.fits <-
      if(!is.null(cl)){
        requireNamespace('parallel')
        if(verbose) {cat("Calling lm/lmrob:\n"); print(gc())}
        out <- parallel::clusterApplyLB(cl, 1:q,
                       function(i){
                         y<-ys[,i]
                         suppressWarnings(try({
                           fit <- if(control$SA.robust) lmrob(y~x,model=FALSE)
                                  else lm(y~x,model=FALSE)
                           
                           list(coef=coef(fit), resid=resid(fit), tvals=summary(fit)$coefficients[,3])},
                                              silent=TRUE))
                         
                       })
        if(verbose) print(gc())
        out
      }else{
        lapply(1:q,
               function(i){
                 y<-ys[,i]
                 suppressWarnings(try({
                   fit <- if(control$SA.robust) lmrob(y~x,model=FALSE)
                          else lm(y~x,model=FALSE)
                   
                   list(coef=coef(fit), resid=resid(fit), tvals=summary(fit)$coefficients[,3])},
                                      silent=TRUE))
               })
      }

    bad.fits <- sapply(h.fits, inherits, "try-error") | apply(diff(ys)==0,2,all)

#    bad.fits <-     # Also, ignore fits where the statistics are too concentrated.    
#      bad.fits | (apply(ys,2,function(y){
#        freqs <- table(y)
#        sum(freqs[-which.max(freqs)])
#      })<nrow(h)/2)

    rm(x, ys); gc()
    
    if(all(bad.fits)) stop("The optimization appears to be stuck. Try better starting parameters, lower SA.init.gain, etc.")
    
    ## Grab the coefficients, t-values, and residuals.
    
    h.nfs <- h.fit <- h.pvals <- h.tvals <- matrix(NA, nrow=p.free+1,ncol=q)

    h.fit[,!bad.fits] <- sapply(h.fits[!bad.fits], "[[", "coef")[seq_len(p.free+1),]
    
    h.resid <- matrix(NA, nrow=n, ncol=q)
    h.resid[,!bad.fits] <- sapply(h.fits[!bad.fits], "[[", "resid")

    h.tvals[,!bad.fits] <- sapply(h.fits[!bad.fits], "[[", "tvals")[seq_len(p.free+1),]

    rm(h.fits)
    
    h.nfs[,!bad.fits] <- matrix(apply(h.resid[,!bad.fits,drop=FALSE],2,function(x) sum(tapply(x,list(tid),length))/sum(tapply(x,list(tid),effectiveSize))), nrow=p.free+1, ncol=sum(!bad.fits), byrow=TRUE)

    h.tvals[,!bad.fits] <- h.tvals[,!bad.fits,drop=FALSE]/sqrt(h.nfs[,!bad.fits,drop=FALSE])

    h.pvals[,!bad.fits] <- 2*pnorm(abs(h.tvals[,!bad.fits,drop=FALSE]),0,1,lower.tail=FALSE)
    
    G.pvals <- t(h.pvals[-1,,drop=FALSE])

    G.pvall <- sort(ifelse(is.na(c(G.pvals)),1,c(G.pvals)))
    
    p.max <- max(G.pvall[G.pvall<=seq_along(G.pvall)/length(G.pvall)*control$SA.phase1.max.q],0)
    
    G.signif <- G.pvals <= p.max
    
    G.signif[is.na(G.signif)] <- FALSE

    ## Compute the variances (robustly) and the statistic weights.
    v <- matrix(NA, q,q)
    v[!bad.fits,!bad.fits] <- if(control$SA.robust) covMcd(h.resid[,!bad.fits,drop=FALSE])$cov else cov(h.resid[,!bad.fits,drop=FALSE])
    v[is.na(v)] <- 0
    v <- v*(nrow(h.resid)-1)/(nrow(h.resid)-p.free-1)
    
    w <- ginv(v)
    
    ## Adjust the number of time steps between jumps.
    edge.ages <- unlist(sapply(states, function(state) state$nw%n%"time"-ergm.el.lasttoggle(state$nw)[,3]+1))
    control$SA.burnin<-control$SA.interval<- round(min(control$SA.max.interval, max(control$SA.min.interval, if(length(edge.ages)>0) control$SA.interval.mul*mean(edge.ages)))/2)
    
    
    if(is.nan(control$SA.burnin)|is.null(control$SA.burnin)|is.na(control$SA.burnin))
        control$SA.burnin <- control$SA.interval <- 10 # TODO: Kirk : check this
    
    
    if(verbose>1){
      cat("New interval:",control$SA.interval ,"\n")
    }
    
    ## Detect parameters whose effect we aren't able to reliably detect.
    ineffectual.pars <- !apply(G.signif,2,any)

    if(all(ineffectual.pars)){
      if(test.G || verbose>0) cat("None of the parameters have a detectable effect. Increasing jitter.\n" )
      control$jitter[!offsets] <- control$jitter[!offsets]*2
    }
    else if(any(ineffectual.pars)){
      if(test.G) cat("Parameters",paste.and(p.names[!offsets][ineffectual.pars]),"do not have a detectable effect. Shifting jitter to them.\n" )
      control$jitter[!offsets] <- control$jitter[!offsets] * (ineffectual.pars+1/2) / mean(control$jitter[!offsets] * (ineffectual.pars+1/2))
    }

    ## Evaluate the dstat/dpar gradient matrix.
    G <- t(h.fit[-1,,drop=FALSE])
    if(test.G)
      G[!G.signif] <- 0
    #else G <- G*(1-G.pvals)
    G[is.na(G)] <- 0

    ## Adaptively scale the estimating equations so that none of the
    ## parameters are "neglected".
    par.eff <- apply(sweep(G,1,ifelse(bad.fits,1,sqrt(diag(v))),"/"),2,function(z)mean(z^2))
    par.eff <- sqrt(par.eff)^control$SA.par.eff.pow
    par.eff <- par.eff / mean(par.eff)

    
    rownames(w)<-colnames(w)<-rownames(v)<-colnames(v)<-q.names
    
    names(par.eff)<-colnames(G.pvals)<-colnames(G)<-p.names[!offsets]
    rownames(G.pvals)<-rownames(G)<-q.names
    if(verbose>1){
      cat("Most recent parameters:\n")
      cat("Formation:\n")
      for(state in states) print(state$eta.form)
      cat("Dissolution:\n")
      for(state in states) print(state$eta.diss)
      cat("Target differences (most recent):\n")
      for(state in states) print(state$nw.diff)
      cat("Target differences (last run):\n")
      print(colMeans(oh.last[,-(1:p),drop=FALSE]))
      cat("Approximate objective function (most recent):\n")
      print(mahalanobis(oh[nrow(oh),-(1:p),drop=FALSE],0,cov=w,inverted=TRUE))
      cat("Approximate objective function (last run):\n")
      print(mahalanobis(colMeans(oh.last[,-(1:p),drop=FALSE]),0,cov=w,inverted=TRUE))
      cat("Estimaged gradient p-values:\n")
      print(G.pvals)
      cat("Estimated gradient:\n")
      print(G)
      cat("Normalized parameter effects:\n")
      print(par.eff)
      cat("Estimated covariance of statistics:\n")
      print(v)
    }

    # Plot if requested.
    if(control$SA.plot.stats && dev.interactive(TRUE)){
      requireNamespace('lattice')
      
      try({
        get.dev("gradients")
        G.scl <- sweep(G, 1, apply(G, 1, function(x) sqrt(mean(x^2))), "/")
        G.scl[is.nan(G.scl)] <- 0
        G.scl <- sweep(G.scl, 2, apply(G.scl, 2, function(x) sqrt(mean(x^2))), "/")
        G.scl[is.nan(G.scl)] <- 0
        print(.my.levelplot(G.scl,main="Scaled Gradients"))
      },silent=TRUE)
      
      suppressWarnings(try({
        get.dev("correlations")
        print(.my.levelplot(cov2cor(v),main="Correlations"))
      },silent=TRUE))
    }
    
    control$GainM <- matrix(0, nrow=p, ncol=q)
    control$GainM[!offsets,] <- t(sweep(G,2,par.eff,"/")) %*% w * control$gain
    control$GainM[!is.finite(control$GainM)] <- 0

    control$dejitter <- matrix(0, nrow=p, ncol=p)
    control$dejitter[!offsets,!offsets] <- control$GainM[!offsets,,drop=FALSE]%*%G

    rownames(control$GainM) <- rownames(control$dejitter) <- colnames(control$dejitter) <- p.names
    colnames(control$GainM) <- q.names
    
    if(verbose>1){
      cat("New deviation -> coefficient map:\n")
      print(control$GainM)
      cat("New jitter cancelation matrix:\n")
      print(control$dejitter)
    }

    if(update.jitter){
      control$jitter[!offsets] <- apply(oh[,1:p,drop=FALSE][,!offsets,drop=FALSE]-history$jitters[,!offsets,drop=FALSE],2,sd)*control$SA.phase2.jitter.mul
      names(control$jitter) <- p.names
    }
    
    if(verbose>1){
      cat("New jitter values:\n")
      print(control$jitter)
    }

    control$dev.guard <- apply(oh[,-(1:p),drop=FALSE],2,function(x) quantile(abs(x),.9)) * control$SA.guard.mul
    if(verbose>1){
      cat("New deviation guard values:\n")
      print(control$dev.guard)
    }

    control$par.guard <- apply(abs(diff(oh[,1:p,drop=FALSE],lag=control$SA.runlength*control$SA.interval-1)),2,median) * control$SA.guard.mul
    if(verbose>1){
      cat("New parameter guard values:\n")
      print(control$par.guard)
    }
    
    list(control=control,
         G=G, w=w, v=v, oh.fit=h.fit, ineffectual.pars=ineffectual.pars, bad.fits=bad.fits)
  }
  
  stergm.EGMME.SA(theta.form0, theta.diss0, nw, model.form, model.diss, model.mon,
                  control, MHproposal.form, MHproposal.diss, eval.optpars, cl=cl,
                  verbose)
}
