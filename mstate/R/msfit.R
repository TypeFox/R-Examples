`msfit` <-
  function(object, newdata, variance=TRUE, vartype=c("aalen","greenwood"), trans)
{
### Modeled after survfit.coxph from the survival library by Therneau;
### indeed much of the code is just copied. Survfit is also called
### for vartype="greenwood".
### Stripped some of the features (cluster, individual) and limited the
### choice of other (type, vartype).
### Input:
###     object: a coxph object
###     newdata: a dataset with K lines for K transitions, as in survfit;
###         newdata should contain a column "strata", specifying to which
###         strata in the coxph object the K lines of newdata correspond.
###         Argument newdata may only be omitted when the formula ends
###         with ~strata(trans)
###     variance: if TRUE (default), estimated (co)variance of cumulative
###         hazards is returned
###     vartype: aalen (default), greenwood is only supported in absence
###         of newdata
###     trans: the transition matrix of the multi-state model
    if(!is.null((object$call)$weights) || !is.null(object$weights))
        stop("msfit cannot (yet) compute the result for a weighted model")
    Terms <- terms(object)
    strat <- attr(Terms, "specials")$strata
    cluster <- attr(Terms, "specials")$cluster
    if (length(cluster)) stop("cluster terms are not supported")
    if (!is.null(attr(object$terms, "specials")$tt))
      stop("msfit cannot yet process coxph models with a tt term")

    resp <-  attr(Terms, "variables")[attr(Terms, "response")]
    nvar <- length(object$coefficients)
    score <- exp(object$linear.predictors) # exp(beta^T Z_i)    
    vartype <- match.arg(vartype)
    if (is.na(vartype)) stop("Invalid variance type specified")

    # Recreate a copy of the data (formerly using coxph.getdata, now with model.frame)
    # Two purposes. One to find largest time point to be used in end result.
    # Second later for call to C.
    has.strata <- !is.null(attr(object$terms, 'specials')$strata)
    if (is.null(object$y) || is.null(object[['x']]) ||
      !is.null(object$call$weights) || 
      (has.strata && is.null(object$strata)) ||
      !is.null(attr(object$terms, 'offset'))) {
      
      mf <- model.frame(object)
    }
    else mf <- NULL  #useful for if statements later
    if (is.null(mf)) y <- object[['y']]
    else {
      y <- model.response(mf)
      y2 <- object[['y']]
      if (!is.null(y2) && any(as.matrix(y2) != as.matrix(y)))
        stop("Could not reconstruct the y vector")
    }
    
    if (is.null(object[['x']])) x <- model.matrix(object, data=mf)
    else x <- object[['x']]
    
    n <- nrow(y)
    if (n != object$n[1] || nrow(x) !=n) 
      stop("Failed to reconstruct the original data set")

    # Here we find largest time point
    type <- attr(y, 'type')
    if (type=='counting') lasty <- max(y[,2]) # start, stop, status
    else if (type=='right') lasty <- max(y[,1]) # time, status
    else stop("Cannot handle \"", type, "\" type survival data")

    if (is.null(mf)) offset <- 0
    else {
      offset <- model.offset(mf)
      if (is.null(offset)) offset <- 0
    }
    
    Terms <- object$terms
    if (!has.strata)  strata <- rep(0L,n)
    else {
      stangle <- untangle.specials(Terms, 'strata') # used multiple times
      strata <- object$strata #try this first
      if (is.null(strata)){
        if (length(stangle$vars)==1) strata <- mf[[stangle$vars]]
        else strata <- strata(mf[, stangle$vars], shortlabel=TRUE)
      }
    }
    if (has.strata) {
      temp <- attr(Terms, "specials")$strata
      factors <- attr(Terms, "factors")[temp,]
      strata.interaction <- any(t(factors)*attr(Terms, "order") >1)
      if (strata.interaction) stop("Interaction terms with strata not supported")
    }
    
    if (vartype=="greenwood") {
      if (missing(trans))
        stop("argument trans missing; needed for vartype=\"greenwood\"")
      labels <- attr(Terms, "term.labels")
      if (length(labels)!=1) stop("Invalid formula for greenwood, ~strata(trans) needed, no covariates allowed")
      if (attr(Terms, "term.labels") != "strata(trans)") # This will happen if length=1 and different from strata(trans)
        stop("Invalid formula for greenwood, ~strata(trans) needed, no covariates allowed")            
      
      sf0 <- summary(survfit(object))
      norisk <- sf0$n.risk
      noevent <- sf0$n.event
      sf0 <- data.frame(time=sf0$time,Haz=-log(sf0$surv),norisk=norisk,
                        noevent=noevent, trans=as.numeric(sf0$strata))
      allt <- sort(unique(c(sf0$time,lasty)))
      nt <- length(allt)
      K <- nrow(to.trans2(trans))
      ### Set up empty structures for Haz and varHaz
      Haz <- data.frame(time=rep(allt,K),Haz=NA,trans=rep(1:K,rep(nt,K)))
      if (variance) {
        tr12 <- data.frame(trans1=rep(1:K,rep(K,K)),trans2=rep(1:K,K))
        tr12 <- tr12[tr12$trans1<=tr12$trans2,]
        varHaz <- data.frame(time=rep(allt,K*(K+1)/2),varHaz=0,
                             trans1=rep(tr12$trans1,rep(nt,K*(K+1)/2)),
                             trans2=rep(tr12$trans2,rep(nt,K*(K+1)/2)))
      }
      ### Easiest to loop over starting states
      S <- nrow(trans)
      for (s in 1:S) {
        trs <- trans[s,]
        trs <- trs[!is.na(trs)]
        ntrs <- length(trs)
        if (ntrs > 0) {
          for (i in 1:ntrs) {
            trans1 <- trs[i]
            sf1 <- sf0[sf0$trans==trans1,]
            Haz$Haz[(trans1-1)*nt + match(sf1$time,allt)] <- sf1$Haz
            Haz$Haz[(trans1-1)*nt + 1:nt] <- NAfix(Haz$Haz[(trans1-1)*nt + 1:nt],subst=0)
            if (variance) {
              varHaz1 <- cumsum((sf1$norisk-sf1$noevent)*sf1$noevent/sf1$norisk^3)
              varHaz11 <- varHaz[varHaz$trans1==trans1 & varHaz$trans2==trans1,]
              varHaz11$varHaz <- NA
              varHaz11$varHaz[match(sf1$time,allt)] <- varHaz1
              varHaz11$varHaz <- NAfix(varHaz11$varHaz,subst=0)
              varHaz[varHaz$trans1==trans1 & varHaz$trans2==trans1,] <- varHaz11
              if (i<ntrs) {
                for (j in ((i+1):ntrs)) {
                  trans2 <- trs[j]
                  sf2 <- sf0[sf0$trans==trans2,]
                  jointt <- intersect(sf1$time,sf2$time)
                  if (length(jointt)>0) {
                    varHazij <- rep(NA,length(jointt))
                    ik <- match(jointt,sf1$time)
                    jk <- match(jointt,sf2$time)
                    varHazij <- cumsum(-sf1$noevent[ik]*sf2$noevent[jk]/sf1$norisk[ik]^3)
                    ## sf1$norisk[ik] should be equal to sf2$norisk[jk]
                    varHaz12 <- varHaz[varHaz$trans1==trans1 & varHaz$trans2==trans2,]
                    varHaz12$varHaz <- NA
                    varHaz12$varHaz[match(jointt,allt)] <- varHazij
                    varHaz12$varHaz <- NAfix(varHaz12$varHaz,subst=0)
                    varHaz[varHaz$trans1==trans1 & varHaz$trans2==trans2,] <- varHaz12
                  }
                }
              }
            }
          }
        }
      }
    }
    else { # aalen
      labels <- attr(Terms, "term.labels")
      if (length(labels)==1) {
        if (labels == "strata(trans)") { # no covariates, no call to C
          sf0 <- summary(survfit(object))
          norisk <- sf0$n.risk
          noevent <- sf0$n.event
          sf0 <- data.frame(time=sf0$time,Haz=-log(sf0$surv),norisk=norisk,
                            noevent=noevent,var=sf0$std.err^2/(sf0$surv)^2,trans=as.numeric(sf0$strata))
          allt <- sort(unique(c(sf0$time,lasty)))
          nt <- length(allt)
          K <- max(sf0$trans)
          ### Set up empty structures for Haz and varHaz
          Haz <- data.frame(time=rep(allt,K),Haz=NA,trans=rep(1:K,rep(nt,K)))
          if (variance) {
            tr12 <- data.frame(trans1=rep(1:K,rep(K,K)),trans2=rep(1:K,K))
            tr12 <- tr12[tr12$trans1<=tr12$trans2,]
            varHaz <- data.frame(time=rep(allt,K*(K+1)/2),varHaz=0,
                                 trans1=rep(tr12$trans1,rep(nt,K*(K+1)/2)),
                                 trans2=rep(tr12$trans2,rep(nt,K*(K+1)/2)))
          }
          for (k in 1:K) {
            # no need for inner loop, since estimated hazards are uncorrelated
            sfk <- sf0[sf0$trans==k,]
            wht <- match(sfk$time,allt)
            Hazk <- Haz[Haz$trans==k,]
            Hazk$Haz[wht] <- sfk$Haz
            Hazk$Haz <- NAfix(Hazk$Haz,subst=0)
            Haz[Haz$trans==k,] <- Hazk
            if (variance) {
              varHazkk <- varHaz[varHaz$trans1==k & varHaz$trans2==k,]
              varHazkk$varHaz <- NA
              varHazkk$varHaz[wht] <- sfk$var
              varHazkk$varHaz <- NAfix(varHazkk$varHaz,subst=0)
              varHaz[varHaz$trans1==k & varHaz$trans2==k,] <- varHazkk
            }
          }
        }
      }
      else {
        method <- object$method
        if (method=="breslow") method <- 1
        else if (method=="efron") method <- 2
        else stop("Only \"efron\" and \"breslow\" methods for ties supported")
        # Copy of the data was already recreated. Finish the job for second purpose.
        # Get the sort index for the data, and add a column to y if
        #  necessary to make it of the "counting process" type
        #  (only 1 C routine to handle both cases).
        #
        type <- attr(y, 'type')
        if (type=='counting') {
          if (has.strata) ord <- order(mf[,strat], y[,2], -y[,3])
          else            ord <- order(y[,2], -y[,3]) 
        }
        else if (type=='right') {
          if (has.strata) ord <- order(mf[,strat], y[,1], -y[,2])
          else            ord <- order(y[,1], -y[,2]) 
          miny <- min(y[,1])
          if (miny < 0) y <- cbind(2*miny -1, y)
          else          y <- cbind(-1, y)
        }
        else stop("Cannot handle \"", type, "\" type survival data")
        
        # Previously, in call to coxph.getdata(), x=variance was used as argument,
        # meaning that data$x would be non-null iff variance=TRUE, hence the previous
        # if (!is.null(data$x)) x <- data$x[ord,], is changed into: if (variance) ...
        if (variance) x <- x[ord,]
        else          x <- 0
        
        # The strata part works differently from survfit; each line in
        # newdata is associated with a single stratum, explicitly given in
        # newdata
        if (has.strata) newstrat <- (as.numeric(mf[,strat]))[ord]
        else newstrat <- rep(1,n)
        newstrat <- cumsum(table(newstrat)) # !!! need to check what happens without strata
        H <- length(newstrat)

        subterms <- function(tt, i) {
          dataClasses <- attr(tt, "dataClasses")
          predvars <- attr(tt, "predvars")
          oldnames <-  dimnames(attr(tt, 'factors'))[[1]]
          tt <- tt[i]
          index <- match(dimnames(attr(tt, 'factors'))[[1]], oldnames)
          if (length(index) >0) {
            if (!is.null(predvars)) 
              attr(tt, "predvars") <- predvars[c(1, index+1)]
            if (!is.null(dataClasses))
              attr(tt, "dataClasses") <- dataClasses[index]
          }
          tt
        }
        
        if (has.strata) {
          temp <- untangle.specials(Terms, 'strata')
          if (length(temp$vars)) 
            Terms <- subterms(Terms, -temp$terms)
        }

        # Get the second, "new" data set
        Terms2 <- delete.response(Terms)
        if (has.strata) {
          if (length(attr(Terms2, 'specials')$strata))
            Terms2 <- subterms(Terms2, -attr(Terms2, 'specials')$strata)
          if (!is.null(object$xlevels)) {
            myxlev <- object$xlevels[match(attr(Terms2, "term.labels"),
                                           names(object$xlevels), nomatch=0)]
            if (length(myxlev)==0) myxlev <- NULL
          }
          else myxlev <- NULL

          mf2 <- model.frame(Terms2, data=newdata, xlev=myxlev)

        }
        else mf2 <- model.frame(Terms2, data=newdata, xlev=object$xlevels)

        offset2 <- 0   #offset variable for the new data set
        if (!missing(newdata)) {
          offset2 <- model.offset(mf2)
          if (length(offset2) >0) offset2 <- offset2 - mean(offset)
          else offset2 <- 0
          x2 <- model.matrix(Terms2, mf2)[,-1, drop=FALSE]  #no intercept
        }
        else stop("newdata missing")
        if (has.strata & is.null(newdata$strata)) stop("no \"strata\" column present in newdata")
        n2 <- nrow(x2)
        
        # Compute risk scores for the new subjects
        coef <- ifelse(is.na(object$coefficients), 0, object$coefficients)
        newrisk <- exp(c(x2 %*% coef) + offset2 - sum(coef*object$means))
        
        dimnames(y) <- NULL # I only use part of Y, so names become invalid
        storage.mode(y) <- 'double'
        ndead <- sum(y[,3])
        untimes <- sort(unique(y[,2][y[,3]==1]))
        nt <- length(untimes)
        
        # n2,nvar are called K,p in 'agmssurv'
        surv <- .C('agmssurv', ### modeled after 'agsurv2' from survival, changes indicated
                   sn = as.integer(n),
                   sp = as.integer(nvar),
                   svar = as.integer(variance),
                   smethod = as.integer(method),
                   sH = as.integer(H),
                   sK = as.integer(n2),
                   snt = as.integer(nt), ## NEW: length of unt
                   y = y[ord,],
                   score = as.double(score[ord]),
                   xmat = as.double(x),
                   varcov = as.double(object$var),
                   strata = as.integer(c(0,newstrat)), ### CHANGED: format to contain 0, end of stratum 1, end of stratum 2 etc in y
                   kstrata = as.integer(newdata$strata), ### NEW: kstrata[k] = stratum in newdata[k,]
                   unt = as.double(untimes), ### NEW: unique event times from all strata
                   newx = as.double(x2),
                   newrisk = as.double(newrisk),
                   Haz = double(nt*n2), ### NEW: OUTPUT, what used to be surv, dimension has changed
                   varHaz = double(nt*n2*(n2+1)/2), ### NEW: OUTPUT, what used to be varhaz, dimension has changed, also contains covariances
                   d = double(3*nvar), ### working vectors d
                   work = double(nt*n2*(nvar+1)) ### NEW: another working vector
        )
        
        # Restructure Haz and varHaz as data frames
        Haz <- data.frame(time=rep(untimes,n2),Haz=surv$Haz,trans=rep(1:n2,rep(nt,n2)))
        varHaz <- as.vector(t(matrix(surv$varHaz,ncol=nt))) # change row and column order
        
        hlp <- matrix(c(rep(1:n2,rep(n2,n2)),rep(1:n2,n2)),n2^2,2)
        hlp <- hlp[hlp[,1]<=hlp[,2],]
        varHaz <- data.frame(time=rep(untimes,n2*(n2+1)/2),varHaz=varHaz,
                             trans1=rep(hlp[,1],rep(nt,n2*(n2+1)/2)),
                             trans2=rep(hlp[,2],rep(nt,n2*(n2+1)/2)))
        
        # Squeeze in the lasty, unless lasty is event time
        if (lasty > max(untimes)) {
          Hmat <- matrix(Haz$Haz,nrow=nt)
          Hmat <- rbind(Hmat,Hmat[nt,])
          vHmat <- matrix(varHaz$varHaz,nrow=nt)
          vHmat <- rbind(vHmat,vHmat[nt,])
          untimes <- c(untimes,lasty)
          nt <- nt+1
          Haz <- data.frame(time=rep(untimes,n2),Haz=as.vector(Hmat),trans=rep(1:n2,rep(nt,n2)))
          varHaz <- data.frame(time=rep(untimes,n2*(n2+1)/2),varHaz=as.vector(vHmat),
                               trans1=rep(hlp[,1],rep(nt,n2*(n2+1)/2)),
                               trans2=rep(hlp[,2],rep(nt,n2*(n2+1)/2)))
        }
      }
    }
    if (variance) res <- list(Haz=Haz,varHaz=varHaz,trans=trans)
    else res <- list(Haz=Haz,trans=trans)
    class(res) <- "msfit"
    return(res)
  }
