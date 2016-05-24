#  File R/stergm.EGMME.SA.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
stergm.EGMME.SA <- function(theta.form0, theta.diss0, nw, model.form, model.diss, model.mon,
                            control, MHproposal.form, MHproposal.diss, eval.optpars, cl=cl,
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

  if(control$SA.plot.progress && !dev.interactive(TRUE)){
    warning("Progress plot requested on a non-interactive graphics device. Ignoring.")
    control$SA.plot.progress <- FALSE # So that we don't print a warning every step.
  }
  
  if(verbose) cat("Starting optimization with with coef_F_0 = (",theta.form0, ") and coef_D_0 = (",theta.diss0,")\n" )
  ###### Define the optimization run function. ######
  
  do.optimization<-function(states, history, control){
    ind.all <- history$ind.all
    tid.all <- history$tid.all
    oh.all <- history$oh.all
    jitters.all <- history$jitters.all
    
    if(verbose) cat("Running stochastic optimization... ")
    zs <- if(!is.null(cl)){
      requireNamespace('parallel')
      # Conveniently, the first argument of stergm.EGMME.SA.Phase2.C
      # is the state of the optimization, so giving clusterApply a
      # list of states will call it for each thread's state.
      if(verbose) {cat("Calling stergm.EGMME.SA.Phase2.C:\n"); print(gc())}
      out <- parallel::clusterApply(cl, states, stergm.EGMME.SA.Phase2.C, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss, control, verbose=verbose)
      if(verbose) print(gc())
      out
    }else{
      list(stergm.EGMME.SA.Phase2.C(states[[1]], model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss, control, verbose=verbose))
    }
    if(verbose) cat("Finished. Extracting.\n")
    for(i in seq_along(states)){

      # Extend the observation index and thread id vectors.
      ind.all <- c(ind.all, sum(tid.all==i) + 1:(control$SA.runlength*control$SA.interval))
      tid.all <- c(tid.all, rep(i, control$SA.runlength*control$SA.interval))
      
      # Extract and store the history of jitters.
      jitters.all <- rbind(jitters.all,zs[[i]]$opt.history[,p+1:p,drop=FALSE])
      colnames(jitters.all) <- p.names
      zs[[i]]$opt.history <- zs[[i]]$opt.history[,-(p+1:p),drop=FALSE]
      
      # Extract and store histroy of trials.
      oh.all <- rbind(oh.all,zs[[i]]$opt.history)
      colnames(oh.all) <- c(p.names,q.names)
    }

    # Limit the size of the full history.
    # TODO: Store a thinned or summarized version of "forgotten" iterations.
    # TODO: Add an option to log it to external storage.
    if(max(ind.all)>control$SA.oh.memory){
      remember.after.ind <- max(ind.all) - control$SA.oh.memory

      # Cut the ind.all last!
      tid.all <-tid.all[ind.all>remember.after.ind]
      oh.all <- oh.all[ind.all>remember.after.ind,,drop=FALSE]
      jitters.all <- jitters.all[ind.all>remember.after.ind,,drop=FALSE]
      ind.all <- ind.all[ind.all>remember.after.ind]

      # Shift the indices so that they would start at 1.
      ind.all <- ind.all - remember.after.ind
    }
    
    history$ind.all <- ind.all
    history$tid.all <- tid.all
    history$jitters.all <- jitters.all
    history$oh.all <- oh.all
   
    min.ind.last <- max(ind.all) - control$SA.runlength*control$SA.interval + 1
    min.ind.keep <- max(ind.all) - max(
                                     max(ind.all)*control$SA.keep.oh,
                                     min(
                                       max(
                                         control$SA.runlength*control$SA.interval*control$SA.keep.min.runs*length(states),
                                         control$SA.keep.min*length(states)
                                         ),
                                       max(ind.all)
                                       )
                                     ) + 1
    
    # Extract and store subhistories of interest.
    
    history$ind <- ind <- ind.all[ind.all>=min.ind.keep]
    history$tid <- tid <-tid.all[ind.all>=min.ind.keep]
    history$ind.last <- ind.last <-ind.all[ind.all>=min.ind.last]
    history$tid.last <- tid.last <- tid.all[ind.all>=min.ind.last]
    history$oh <- oh <- oh.all[ind.all>=min.ind.keep,,drop=FALSE]
    history$oh.last <- oh.last <- oh.all[ind.all>=min.ind.last,,drop=FALSE]
    history$jitters <- jitters <- jitters.all[ind.all>=min.ind.keep,,drop=FALSE]
    history$jitters.last <- jitters.last <- jitters.all[ind.all>=min.ind.last,,drop=FALSE]      

    rm(ind.all, tid.all, jitters.all, oh.all); gc()

    # Plot if requested.
    if(control$SA.plot.progress && dev.interactive(TRUE)){
      requireNamespace('lattice')
      
      get.dev("progress.plot")
      
      thin <- (nrow(oh)-1)%/%(control$SA.max.plot.points/length(states)) + 1
      cols <- floor(sqrt(ncol(oh)))
      layout <- c(cols,ceiling(ncol(oh)/cols))
      
      suppressWarnings(print(lattice::xyplot(window(do.call(mcmc.list,by(as.data.frame(oh),INDICES=list(tid=tid),mcmc,start=min.ind.keep)), thin=thin), panel = function(...) {lattice::panel.xyplot(...);lattice::panel.abline(0, 0)}, as.table = TRUE, layout = layout, xlab=NULL)))
    }
    
    # Extract and return the "states" and the "history".
    list(states = lapply(zs, function(z) list(nw = z$newnetwork,
           nw.diff = z$nw.diff,
           eta.form = z$eta.form,
           eta.diss = z$eta.diss)
           ),
         history = history)
  }

  # This function essentially runs the chain forward with gradient etc. set to 0.
  do.dummy.run <- function(states, history, control, burnin, steps, eta.form=NULL, eta.diss=NULL){
        control$GainM <- matrix(0, nrow=p, ncol=q)
        control$dejitter <- matrix(0, nrow=p, ncol=p) # Dejitter tries to cancel the effect of jitter on the optimizer.
        control$par.guard <- control$jitter<-rep(0,p)
        control$dev.guard <- rep(0,q)
        control$dev.guard[] <- 1e10 # A huge value
        control$par.guard[!offsets] <- 1e10

        
        control$SA.keep.min.runs <- 0

        control$SA.runlength <- 1
        control$SA.burnin <- burnin
        control$SA.interval <- steps


        
        for(i in seq_along(states)){
          if(!is.null(eta.form)) states[[i]]$eta.form <- eta.form
          if(!is.null(eta.diss)) states[[i]]$eta.diss <- eta.diss
        }

        do.optimization(states,history,control)
  }
  
  interpolate.par <- function(h.fit, w=diag(1,nrow=ncol(h.fit))){
    x <- t(h.fit[-1,,drop=FALSE])
    y <- -cbind(h.fit[1,])

    c(solve(t(x)%*%w%*%x)%*%t(x)%*%w%*%y)
  }

  V.sandwich <- function(w, G, V.stat=ginv(w)){
    solve(t(G)%*%w%*%G)%*%t(G)%*%w%*%V.stat%*%w%*%G%*%solve(t(G)%*%w%*%G)
  }

  best.states <- function(states, history){
    w <- ginv(cov(history$oh.all[,-(1:p),drop=FALSE]))
    best.i <- which.min(mahalanobis((history$oh.all[-1,-(1:p),drop=FALSE]+history$oh.all[-nrow(history$oh.all),-(1:p),drop=FALSE])/2,0,w,inverted=TRUE))
    best.par <- history$oh.all[best.i,1:p][!offsets]

    lapply(states, function(state){
      if(p.form.free) state$eta.form[!model.form$etamap$offsettheta] <- best.par[seq_len(p.form.free)]
      if(p.diss.free) state$eta.diss[!model.diss$etamap$offsettheta] <- best.par[p.form.free+seq_len(p.diss.free)]
      state
    })
  }


  
  ##### Construct the initial state. ######

  history <- list()
  history$ind.all <- history$tid.all <- history$oh.all <- history$jitters.all <- state <- NULL

  for(restart in 1:control$SA.restarts){    
    nw.diff <- model.mon$nw.stats - model.mon$target.stats  # nw.diff keeps track of the difference between the current network and the target statistics.
    
    states <- replicate(if(!is.null(cl)) control$parallel else 1,
                        {
                          list(nw=nw,
                               eta.form = ergm.eta(theta.form0, model.form$etamap),
                               eta.diss = ergm.eta(theta.diss0, model.diss$etamap),
                               nw.diff  = nw.diff)
                        },
                        simplify=FALSE
                        )
    
    cat('========  Phase 1: Burn in, get initial gradient values, and find a configuration under which all targets vary. ========\n',sep="")
    
    ###### Set up and run the burn-in. ######

    control$changes <- control$collect.form <- control$collect.diss <- FALSE
    control.phase1<-control
    control.phase1$time.samplesize <- 1
    control.phase1$time.burnin <- control$SA.burnin
    control.phase1$time.interval <- 1
    
    cat("Burning in... ")
    
    zs <- if(!is.null(cl)){
      requireNamespace('parallel')
      if(verbose) {cat("Calling stergm.getMCMCsample:\n"); print(gc())}
      out <- parallel::clusterApply(cl, seq_along(states), function(i) stergm.getMCMCsample(states[[i]]$nw, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss, states[[i]]$eta.form, states[[i]]$eta.diss, control.phase1, verbose))
      if(verbose) print(gc())
      out
    }else{
      list(stergm.getMCMCsample(states[[1]]$nw, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss, states[[1]]$eta.form, states[[1]]$eta.diss, control.phase1, verbose))
    }
    
    cat("Done.\n")
    
    # Update the state with burn-in results.
    for(i in seq_along(zs)){
      states[[i]]$nw <- zs[[i]]$newnetwork
      states[[i]]$nw.diff <- states[[i]]$nw.diff + zs[[i]]$statsmatrix.mon[NROW(zs[[i]]$statsmatrix.mon),]
    }
    
    ###### Gradient estimation and getting to a decent configuration. ######
    
    ## Here, "decent configuration" is defined as one where all target
    ## statistics have some variability --- none are stuck.
    
    # Start with jitter-only.
    
    control$GainM <- matrix(0, nrow=p, ncol=q)
    control$dejitter <- matrix(0, nrow=p, ncol=p) # Dejitter tries to cancel the effect of jitter on the optimizer.
    
    control$par.guard <- control$jitter<-rep(0,p)
    control$dev.guard <- rep(0,q)
    control$jitter[!offsets] <- control$SA.phase1.jitter
    control$dev.guard[] <- 1e10 # A huge value
    control$par.guard[!offsets] <- control$SA.phase1.jitter * 4
    
    ## Adjust the number of time steps between jumps using burn-in.
    edge.ages <- unlist(sapply(states, function(state) state$nw%n%"time"-ergm.el.lasttoggle(state$nw)[,3]+1))
    control$SA.burnin<-control$SA.interval<- round(min(control$SA.max.interval, max(control$SA.min.interval, if(length(edge.ages)>0) control$SA.interval.mul*mean(edge.ages)))/2)
  
  if(is.nan(control$SA.burnin)|is.null(control$SA.burnin)|is.na(control$SA.burnin))
    control$SA.burnin <- control$SA.interval <- 10 # TODO: Kirk : check this
  
    if(verbose>1){
      cat("New interval:",control$SA.interval ,"\n")
    }  
    
    for(try in 1:control$SA.phase1.tries){
      if(verbose) cat('======== Attempt ',try,' ========\n',sep="") else cat('Attempt',try,':\n')
      for(run in 1:control$SA.phase1.minruns){
        tmp <- if(control$SA.restart.on.err) try(do.optimization(states, history, control), silent=!verbose) else do.optimization(states, history, control)
        if(inherits(tmp, "try-error") || all(apply(tmp$history$oh.last[,-(1:p),drop=FALSE],2,var)<sqrt(.Machine$double.eps))){
          cat("Something went very wrong. Restarting with smaller gain.\n")
          control$SA.init.gain <- control$SA.init.gain * control$SA.gain.decay
          do.restart <- TRUE
          break
        }else{
          states <- tmp$states
          history <- tmp$history
          do.restart <- FALSE
        }
      }

      states <- best.states(states, history)

      control$gain <- control$SA.init.gain
      out <- if(control$SA.restart.on.err) try(eval.optpars(states, history, control, TRUE,restart>1,FALSE), silent=!verbose) else eval.optpars(states, history, control, TRUE,restart>1,FALSE)
      if(inherits(out, "try-error") || all(apply(history$oh.last[,-(1:p),drop=FALSE],2,var)<sqrt(.Machine$double.eps))){
        cat("Something went very wrong. Restarting with smaller gain.\n")
        control$SA.init.gain <- control$SA.init.gain * control$SA.gain.decay
        do.restart <- TRUE
        break
      }else do.restart <- FALSE
      control <- out$control

      if(all(!out$ineffectual.pars) && all(!out$bad.fits)){
        cat("All parameters have some effect and all statistics are moving. Proceeding to Phase 2.\n")
        break
      }
      if(try==control$SA.phase1.tries) stop("The optimizer was unable to find a reasonable configuration: one or more statistics are still stuck after multiple tries, and one or more parameters do not appear to have any robust effect.")
    }
    if(do.restart) next
    
    ###### Main optimization run. ######
    
    cat('========  Phase 2: Find and refine the estimate. ========\n',sep="")
    
    for(subphase in 1:control$SA.phase2.levels.max){
      if(verbose) cat('======== Subphase ',subphase,' ========\n',sep="") else cat('Subphase 2.',subphase,' ',sep="")
      
      control$gain <- control$SA.init.gain*control$SA.gain.decay^(subphase-1)
      stepdown.count <- control$SA.stepdown.ct
      
      for(regain in 1:control$SA.phase2.repeats){
        tmp <- if(control$SA.restart.on.err) try(do.optimization(states, history, control), silent=!verbose) else do.optimization(states, history, control)
        if(inherits(tmp, "try-error") || all(apply(tmp$history$oh.last[,-(1:p),drop=FALSE],2,var)<sqrt(.Machine$double.eps))){
          cat("Something went very wrong. Restarting with smaller gain.\n")
          control$SA.init.gain <- control$SA.init.gain * control$SA.gain.decay
          do.restart <- TRUE
          break
        }else{
          states <- tmp$states
          history <- tmp$history
          do.restart <- FALSE
        }
        
        if(verbose){
          cat("New parameters:\n")
          cat("Formation:\n")
          for(state in states) print(state$eta.form)
          cat("Dissolution:\n")
          for(state in states) print(state$eta.diss)
        }

        ## Get updated gain and other values
        out <- if(control$SA.restart.on.err) try(eval.optpars(states, history, control, FALSE,TRUE,TRUE), silent=!verbose) else eval.optpars(states, history, control, FALSE,TRUE,TRUE)
        if(inherits(out, "try-error") || all(apply(history$oh.last[,-(1:p),drop=FALSE],2,var)<sqrt(.Machine$double.eps))){
          cat("Something went very wrong. Restarting with smaller gain.\n")
          control$SA.init.gain <- control$SA.init.gain * control$SA.gain.decay
          do.restart <- TRUE
          break
        }else do.restart <- FALSE
        control <- out$control

        ## Run two tests here:
        ## 1) Are the estimating equations, on average, 0?
        ## 2) Does there appear to be a trend in their sum of squares?

        # Calculate the approximate estimating equation values:
        ys <- history$oh[,-(1:p),drop=FALSE]%*%out$w%*%out$G

        # This test is fast, so no need to thin.
        p.val.1 <- try(approx.hotelling.diff.test(ys)$p.value)
        if(is.na(p.val.1)) p.val.1 <- 0

        # Thin the data to keep from bogging down.
        x <- unique(round(seq(from=1,to=NROW(history$oh),length.out=control$SA.stepdown.maxn)))
        y <- sqrt(mahalanobis(ys,0,ginv(cov(ys)),inverted=TRUE)[x])
        i <- history$ind[x]
        t <- history$tid[x]

        for(thread in unique(t)) i[t==thread] <- rank(i[t==thread])

        fit.2 <- try(summary(gls(y~x,correlation=corAR1(form=~i|t)))$tTable[2,c(1,4)])

         if(!inherits(p.val.1, "try-error") && !inherits(fit.2, "try-error")){
          p.val.2 <- fit.2[2]/2 # We are interested in one-sided decline here.
          est.2 <- fit.2[1]
          if(est.2>0) p.val.2 <- 1-p.val.2 # If it's actually getting worse one-sided p-value is thus.
          
          if(verbose){
            cat("Estimating equations = 0 p-value:",p.val.1,", trending:", p.val.2, ".\n")
          }

          p.vals <- c(p.val.1,p.val.2)

          fisher.pval <- function(p.vals){
            p.vals <- unlist(p.vals)
            df  <- 2*length(p.vals)
            pchisq(-2*sum(log(p.vals)), df, lower.tail=FALSE)
          }
          if(fisher.pval(p.vals)>control$SA.stepdown.p){
            stepdown.count <- stepdown.count - 1
            if(stepdown.count<=0){
              if(verbose) cat("Estimating equations do not significantly differ from 0 and neither they nor the parameters exhibit a significant trend. Reducing gain.\n")
              else cat("\\")
              stepdown.count <- control$SA.stepdown.ct
              if(!verbose) cat("\n")
              break
            }else{
              if(verbose) cat("Estimating equations do not significantly differ from 0 and do not exhibit a significant trend. ",stepdown.count,"/",control$SA.stepdown.ct," to go.\n")
              else cat("\\")
            }
          }else{
            stepdown.count <- control$SA.stepdown.ct
            if(verbose) cat("Estimating equations significantly differ from 0 or exhibit a significant trend. Resetting counter.\n")
            else cat("/")
          }
        }else{
          if(verbose) cat("Problem testing estimating equations. Continuing with current gain.\n")
          else cat("!/")

        }
        if(do.restart) break
      }

      if(do.restart) break      

      
      #### Decide whether to stop.

      ## Run through minimal number of phases.
      if(subphase<control$SA.phase2.levels.min) next
      
      ## Run three tests here:
      ## 1) Is the stochastic approximation estimate of the GMM sufficiently precise?
      ## 2) Is there strong evidence of nonlinearity?
      ## 3) Does there appear to be a trend in parameter values?

      # Test 1:
      
      par.se <- {
        h <- history$oh[,1:p,drop=FALSE][,!offsets,drop=FALSE]
        apply(h,2,function(x) sd(x)/sqrt(effectiveSize(x)))
      }

      sandwich.se <- sqrt(diag(V.sandwich(out$w,out$G,out$v)))

      if(verbose){
        cat("Approximate standard error of the estimate:\n")
        print(sandwich.se)        
        cat("Approximate standard error of window means:\n")
        print(par.se)
        cat("par. var. / (std. var. + par. var.):\n")
        print(par.se^2/(sandwich.se^2+par.se^2))
      }
      
      # Test 2:

      ys <- history$oh[, -(1:p), drop=FALSE]
      xs <- history$oh[, 1:p, drop=FALSE][,!offsets,drop=FALSE]

      nlin.totest <- rep(c(FALSE,TRUE),c(p.free+1, p.free*2))
      
      nlin.fit <- lm(ys~xs+I(xs^2)+I(xs^3))
      nlin.nfs <- matrix(apply(cbind(resid(nlin.fit)),2,function(x) sum(tapply(x,list(history$tid),length))/sum(tapply(x,list(history$tid),effectiveSize))), nrow=1+p.free*3, ncol=q, byrow=TRUE)
      nlin.coef <- cbind(coef(nlin.fit))
      nlin.vcov <- vcov(nlin.fit)
      
      drop <- apply(is.na(nlin.coef),1,any)
      nlin.coef <- nlin.coef[nlin.totest & !drop,,drop=FALSE]
      nlin.nfs <- nlin.nfs[nlin.totest & !drop,,drop=FALSE]
      nlin.vcov <- t(nlin.vcov[nlin.totest[!drop],nlin.totest[!drop],drop=FALSE]*sqrt(c(nlin.nfs)))*sqrt(c(nlin.nfs))
      chi2 <- mahalanobis(c(nlin.coef),0,ginv(nlin.vcov),inverted=TRUE)

      p.val.1 <- pchisq(chi2, length(nlin.coef), lower.tail=FALSE)
      if(verbose) cat("Local nonlinearity p-value:",p.val.1,"\n")

      if(subphase == control$SA.phase2.levels.max){
        if(verbose) cat("Maximum number of gain levels exceeded. Stopping.")
      }else{
        if(all(par.se^2/(sandwich.se^2+par.se^2) > control$SA.phase2.max.mc.se)){
          if(verbose) cat("EGMME does not appear to be estimated with sufficient prescision. Continuing.\n")
          next
        }
        
        if(p.val.1 < control$SA.stop.p){
          if(verbose) cat("There is evidence of local nonlinearity. Continuing.\n")
          next
        }
      }

      ###### If we've gotten this far, proceed to Part 3. ######

      ## Refine the estimate.
      if(inherits(state, "try-error")) stop("Something went wrong too many times. Try better starting values or reducing control$SA.init.gain.") 
      
      eta.form <- states[[1]]$eta.form
      eta.diss <- states[[1]]$eta.diss
      
      eta.free <- switch(control$SA.refine,
                         mean = colMeans(history$oh[,1:p,drop=FALSE][,!offsets,drop=FALSE]),
                         linear = interpolate.par(out$oh.fit,out$w),
                         none = if(is.null(cl)) c(eta.form,eta.diss)[!offsets] else stop("No interpolation does not make sense with multithreaded fitting."))
      if(p.form.free) eta.form[!model.form$etamap$offsettheta] <- eta.free[seq_len(p.form.free)]
      if(p.diss.free) eta.diss[!model.diss$etamap$offsettheta] <- eta.free[p.form.free+seq_len(p.diss.free)]
      
      if(verbose){
        cat("Refining the estimate using the", control$SA.refine,"method. New estimate:\n")
        cat("Formation:\n")
        print(eta.form)
        cat("Dissolution:\n")
        print(eta.diss)
      }
      
      ## Sample to estimate standard error and test if we've actually arrived.
      G <- out$G
      w <- out$w
      V.par<-matrix(NA,p,p)
      if(control$SA.se){
        cat('========  Phase 3: Simulate from the fit and estimate standard errors. ========\n',sep="")
        
        if(verbose)cat("Evaluating target statistics at the estimate.\n")

        # This is to avoid the latest sample "pushing out" everything else.
        control$SA.keep.min <- round(nrow(history$oh)/length(states))+control$SA.phase3.samplesize.runs*control$SA.runlength*control$SA.interval
        tmp <- do.dummy.run(states, history, control, control$SA.burnin, control$SA.phase3.samplesize.runs*control$SA.runlength*control$SA.interval, eta.form=eta.form, eta.diss=eta.diss)
        if(verbose)cat("Finished.\n")
        states <- tmp$states
        history <- tmp$history
        rm(tmp); gc()
        
        min.ind.se <- max(history$ind.all)-control$SA.phase3.samplesize.runs*control$SA.runlength*control$SA.interval + 1
        sm.mon <- history$oh.all[history$ind.all>=min.ind.se,-(1:p),drop=FALSE]
        sm.tid <- history$tid.all[history$ind.all>=min.ind.se]

        V.stat<-cov(sm.mon)
        ee <- sm.mon %*% ginv(V.stat) %*% G
        p.val.1 <- approx.hotelling.diff.test(ee)$p.value

        if(verbose) cat("Estimating equation = 0 p-value:",p.val.1,"\n")

        if(subphase == control$SA.phase2.levels.max){
          if(verbose) cat("Maximum number of gain levels exceeded. Stopping.")
        }else{
          if(p.val.1 < control$SA.stop.p){
            if(verbose) cat("Simulated values of estimating equations are not centered around 0. Continuing.\n")
            next
          }else{
            if(verbose) cat("Simulated values of estimating equations are centered around 0. Stopping.\n")          
          }
        }
        
        V.par[!offsets,!offsets]<-V.sandwich(w,G,V.stat)
        sm.mons <- lapply(unique(sm.tid), function(t) sm.mon[sm.tid==t,,drop=FALSE])
        sm.mon <- do.call(mcmc.list, lapply(sm.mons, mcmc, start=control$SA.burnin))
      }else{
        V.par[!offsets,!offsets]<-V.sandwich(w,G,out$v)
        sm.mon <- NULL
      }

      mc.se <- rep(NA, p)
      mc.se[!offsets] <- par.se
      mc.se.form <- mc.se[1:p.form]
      mc.se.diss <- mc.se[-(1:p.form)]

      break
    }
    if(!do.restart) break # If we've gotten this far, no restart conditions have been triggered, so we are good.
  }
  
    
  list(newnetwork=if(control$SA.se) zs[[1]]$newnetwork else states[[1]]$nw,
       newnetworks=if(control$SA.se) lapply(zs,"[[","newnetwork") else lapply(states,"[[","nw"),
       init.form=theta.form0,
       init.diss=theta.diss0,
       covar=V.par,
       covar.form=V.par[seq_len(p.form),seq_len(p.form),drop=FALSE],
       covar.diss=V.par[p.form+seq_len(p.diss),p.form+seq_len(p.diss),drop=FALSE],
       mc.se = mc.se,
       mc.se.form = mc.se.form,
       mc.se.diss = mc.se.diss,
       eta.form=eta.form,
       eta.diss=eta.diss,
       opt.history=history$oh.all,
       sample=sm.mon,
       network=nw)
}

stergm.EGMME.SA.Phase2.C <- function(state, model.form, model.diss, model.mon,
                             MHproposal.form, MHproposal.diss, control, verbose) {
  Clist.form <- ergm.Cprepare(state$nw, model.form)
  Clist.diss <- ergm.Cprepare(state$nw, model.diss)
  Clist.mon <- ergm.Cprepare(state$nw, model.mon)
  maxedges <- max(control$MCMC.init.maxedges, Clist.mon$nedges)
  maxchanges <- max(control$MCMC.init.maxchanges, Clist.mon$nedges)

  eta.form <- ergm:::.deinf(state$eta.form)
  eta.diss <- ergm:::.deinf(state$eta.diss)
  
  repeat{
    z <- .C("MCMCDynSArun_wrapper",
            # Observed/starting network. 
            as.integer(Clist.form$tails), as.integer(Clist.form$heads),
            time = as.integer(min(Clist.form$time,Clist.diss$time,Clist.diss$time)),
            lasttoggle = as.integer(NVL(Clist.form$lasttoggle,Clist.diss$lasttoggle,Clist.mon$lasttoggle,0)),  
            as.integer(Clist.form$nedges),
            as.integer(Clist.form$n),
            as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
            # Formation terms and proposals. 
            as.integer(Clist.form$nterms), as.character(Clist.form$fnamestring), as.character(Clist.form$snamestring),
            as.character(MHproposal.form$name), as.character(MHproposal.form$pkgname),
            as.double(Clist.form$inputs),
            # Dissolution terms and proposals. 
            as.integer(Clist.diss$nterms), as.character(Clist.diss$fnamestring), as.character(Clist.diss$snamestring),
            as.character(MHproposal.diss$name), as.character(MHproposal.diss$pkgname),
            as.double(Clist.diss$inputs),
            # Parameter fitting.
            eta=as.double(c(eta.form, eta.diss)),
            as.integer(Clist.mon$nterms), as.character(Clist.mon$fnamestring), as.character(Clist.mon$snamestring),
            as.double(Clist.mon$inputs), 
            nw.diff=as.double(state$nw.diff),
            as.integer(control$SA.runlength),
            as.double(control$GainM),
            as.double(control$jitter), as.double(control$dejitter), # Add a little bit of noise to parameter guesses.
            as.double(control$dev.guard),
            as.double(control$par.guard),
            # Degree bounds.
            as.integer(MHproposal.form$arguments$constraints$bd$attribs), 
            as.integer(MHproposal.form$arguments$constraints$bd$maxout), as.integer(MHproposal.form$arguments$constraints$bd$maxin),
            as.integer(MHproposal.form$arguments$constraints$bd$minout), as.integer(MHproposal.form$arguments$constraints$bd$minin),
            as.integer(MHproposal.form$arguments$constraints$bd$condAllDegExact), as.integer(length(MHproposal.form$arguments$constraints$bd$attribs)), 
            # MCMC settings.              
            as.integer(control$SA.burnin),
            as.integer(control$SA.interval),
            as.integer(control$MCMC.burnin.min), as.integer(control$MCMC.burnin.max), as.double(control$MCMC.burnin.pval), as.double(control$MCMC.burnin.add),
            # Space for output.
            as.integer(maxedges),
            as.integer(maxchanges),
            newnwtails = integer(maxedges), newnwheads = integer(maxedges), 
            opt.history=double(((Clist.form$nstats+Clist.diss$nstats)*2+Clist.mon$nstats)*control$SA.runlength*control$SA.interval),
            # Verbosity.
            as.integer(max(verbose-1,0)),
            status = integer(1), # 0 = OK, MCMCDyn_TOO_MANY_EDGES = 1, MCMCDyn_MH_FAILED = 2, MCMCDyn_TOO_MANY_CHANGES = 3
            PACKAGE="tergm")
    if(z$status==0) break;
    if(z$status==1){
      maxedges <- 5*maxedges
      if(verbose>0) message("Too many edges encountered in the simulation. Increasing capacity to ", maxedges)
    }
    if(z$status==3){
      maxchanges <- 5*maxchanges
      if(verbose>0) message("Too many changes elapsed in the simulation. Increasing capacity to ", maxchanges)
    }
  }

  eta.form <- z$eta[seq_len(Clist.form$nstats)]
  names(eta.form) <- model.form$coef.names
  eta.diss <- z$eta[-seq_len(Clist.form$nstats)]
  names(eta.diss) <- model.diss$coef.names

  newnetwork<-newnw.extract(state$nw,z)
  newnetwork %n% "time" <- z$time
  newnetwork %n% "lasttoggle" <- z$lasttoggle
  
  list(nw.diff=z$nw.diff,
       newnetwork=newnetwork,
       eta.form=eta.form,
       eta.diss=eta.diss,
       opt.history=matrix(z$opt.history,ncol=(Clist.form$nstats+Clist.diss$nstats)*2+Clist.mon$nstats,byrow=TRUE))
}
