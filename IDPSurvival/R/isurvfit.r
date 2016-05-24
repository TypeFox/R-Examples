## File IDPSurvival/R/isurvfit.r
##
## IDPSurvival package for R (http://www.R-project.org)
##
## Copyright (C) 2014 Dalle Molle Institute for Artificial Intelligence.
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

isurvfit <- function(formula, data, s=0.5, weights, subset, display=TRUE,
                     conf.type=c('exact',  'approx', 'none'), nsamples=10000,
                     conf.int= .95) {
  Call <- match.call()
  Call[[1]] <- as.name('isurvfit')  #make nicer printout for the user
  # create a copy of the call that has only the arguments we want,
  #  and use it to call model.frame()
  indx <- match(c('formula', 'data', 'weights', 'subset'), names(Call), nomatch=0)
  # It's very hard to get the next error message other than malice
  #  eg survfit(wt=Surv(time, status) ~1) 
  if (indx[1]==0) stop("a formula argument is required")
  
  conf.type <- match.arg(conf.type)
  if (conf.type=='none') {se.fit=FALSE
  }  else {se.fit=TRUE}
   
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  m <- eval.parent(temp)
  
  Terms <- terms(formula)
#   ord <- attr(Terms, 'order')
#   if (length(ord) & any(ord !=1))
#     stop("Interaction terms are not valid for this function")
  
  n <- nrow(m)
  Y <- model.extract(m, 'response')
  if (!is.Surv(Y)) stop("Response must be a survival object")
  
  casewt <- model.extract(m, "weights")
  if (is.null(casewt)) casewt <- rep(1,n)
  
#   if (!is.null(attr(Terms, 'offset'))) warning("Offset term ignored")
  
#   temp <- untangle.specials(Terms, "cluster")
#   if (length(temp$vars)>0) {
#     if (length(temp$vars) > 1) stop("can not have two cluster terms")
#     Terms <- Terms[-temp$terms]
#   }
  
  ll <- attr(Terms, 'term.labels')
  if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
  else X <- strata(m[ll])
  
  if (!is.Surv(Y)) stop("y must be a Surv object")
  
  if (attr(Y, 'type') == "right" )
     temp <- .isurvfit.default(X,Y,s, casewt,se.fit,conf.type,nsamples,conf.int)
  else {
    stop("Can only handle right censored data")
  }
  class(temp) <- 'isurvfit'
  
  if (display) {
    plot(temp,se.fit=se.fit)
  }
  temp$call <- Call
  return(temp)
}

.isurvfit.default <- function(x, y, s, casewt=rep(1,length(x)),se.fit,
                       conf.type, nsamples=10000, conf.int= .95) {
  if (is.logical(conf.int)) {
    # A common error is for users to use "conf.int = FALSE"
    #  it's illegal, but allow it
    if (!conf.int) conf.type <- "none"
    conf.int <- .95
  }
  
  if (!is.Surv(y)) stop("y must be a Surv object")
  if (!is.factor(x)) stop("x must be a factor")
  if (attr(y, 'type') != 'right')
    stop("Can only handle right censored data")
  xlev <- levels(x)   # Will supply names for the curves
  x <- as.numeric(x)  # keep only the levels
  
  n.used <- as.vector(table(x))    # This is for the printout
  nstrat <- length(n.used)
  
  
  # Each of the necessary output objects is originally a list with one
  #  element per strata.  
  time   <- vector('list', nstrat)
  n.risk <- vector('list', nstrat)
  survUP <- vector('list', nstrat)
  survLOW <- vector('list', nstrat)
  survLOW0 <- vector('list', nstrat)
  n.cens <- vector('list', nstrat)
  n.event<- vector('list', nstrat)
  upper<- vector('list', nstrat)
  lower<- vector('list', nstrat)
  lower0<- vector('list', nstrat)
  varhazLOW0 <- vector('list', nstrat)
  strata <- integer(nstrat)
  if (se.fit) {
    varhazUP <- vector('list', nstrat)
    varhazLOW <- vector('list', nstrat)
  }
  
  uniquex <- sort(unique(x))
  qval = floor(nsamples*(1-conf.int)/2)
  
  for (i in 1:nstrat) {
    who <- (x== uniquex[i])

      temp <- tapply(casewt[who], 
                     list(factor(y[who,1]), 
                          factor(y[who,2], levels=0:1)), sum)
      temp <- ifelse(is.na(temp), 0, temp)
      
      # The two lines below do not always give the same answer
      #  When two times differ by only the machine precision, unique
      #  will give more values, and then "time" will be a different
      #  length than the other components
      time[[i]] <- type.convert(dimnames(temp)[[1]], as.is=TRUE,
                                dec=getOption("OutDec"))
      #	    time[[i]] <- sort(unique(y[who,1]))   # old version of line above
      ntemp  <- (dim(temp))[1]
      nevent <- as.vector(temp[,2])
      ncens  <- as.vector(temp[,1])
      nrisk  <- rev(cumsum(rev(temp %*% c(1,1))))
      ndead  <- as.vector(table(y[who,1], factor(y[who,2], 
                                                 levels=0:1)) [,2])
    
    strata[i] <- ntemp
    trisk <- ifelse(nrisk==0, 1, nrisk) #avoid 0/0 cases
    
    survUP[[i]] <- cumprod((s+trisk-nevent)/(s+trisk))
    survLOW[[i]] <- (trisk-nevent)/(s+trisk-nevent)* cumprod((s+trisk-nevent)/(s+trisk))
    survLOW0[[i]] <- trisk[1]/(s+trisk[1])
    
    if (conf.type=='exact') {
      a=t(cbind(nevent, ncens))
      a=as.vector(a)
      Hk <- rdirichlet(nsamples,cbind(t(a),s))
      Hrisk <- 1-t(apply(Hk, 1, cumsum))
      Hrisk[which(Hrisk<0)] <- 0 
      Hs <- .repmat(matrix(Hk[,2*ntemp+1],nsamples,1),1,ntemp)
      a = seq(1,2*ntemp,by=2)
      Hrisk = Hrisk[,a]
      Hevent = Hk[,a]
      termL <-(Hrisk-Hs)/(Hrisk)
      termL[which(Hrisk==0)]<-1
      termL[which(Hrisk-Hs<=0)]<-0
      term  <- (Hrisk)/(Hrisk+Hevent)
      term[which(Hrisk+Hevent ==0)] <- 1
      term[which(Hrisk == 0)] <- 0
      if (ntemp==1) {
        SsampUP <- cumprod(term)
        SsampLOW <- termL*cumprod(term)
        varhazUP[[i]] <- var(SsampUP)
        varhazLOW[[i]] <- var(SsampLOW)
        SsampLOW0 <- (Hrisk+Hevent)/(Hrisk+Hevent+Hs)
      } else {
        SsampUP <- t(apply(term, 1, cumprod))
        SsampLOW <- termL*t(apply(term, 1, cumprod))
        varhazUP[[i]] <- apply(SsampUP, 2, var)
        varhazLOW[[i]] <- apply(SsampLOW, 2, var)
        SsampLOW0 <- (Hrisk[,1]+Hevent[,1])/(Hrisk[,1]+Hevent[,1]+Hs[,1])
      }
      

      
      if (qval==0) {
        upper[[i]]<- rep(1,ntemp)
        lower[[i]]<- rep(0,ntemp)
        lower0[[i]]<- 0
      } else {
        if (ntemp==1) {
          SsampUP <- sort(SsampUP, decreasing = TRUE, method='quick')
          upper[[i]]<-SsampUP[qval]
          SsampLOW <- sort(SsampLOW, decreasing = FALSE,method='quick')
          lower[[i]]<-SsampLOW[qval] 
        } else {
          SsampUP <- apply(SsampUP, 2, sort,decreasing = TRUE, method='quick')
          upper[[i]]<-SsampUP[qval,]
          SsampLOW <- apply(SsampLOW, 2, sort,decreasing = FALSE ,method='quick')
          lower[[i]]<-SsampLOW[qval,]
        }
        SsampLOW0 <- sort(SsampLOW0,decreasing = FALSE,method='quick')
        lower0[[i]]<-SsampLOW[qval]
      }
    }

    if (conf.type=='approx') {
      varhazUP[[i]] <- c(cumprod((s+trisk-nevent)/(s+trisk))*cumprod((1+s+trisk-nevent)/(1+s+trisk)))
      varhazLOW[[i]] <- (trisk-nevent)/(s+trisk-nevent)*(trisk-nevent+1)/(s+trisk-nevent+1)*
                            cumprod((s+trisk-nevent)/(s+trisk))*cumprod((1+s+trisk-nevent)/(1+s+trisk))
      varhazLOW0[[i]] <- trisk[1]/(trisk[1]+s)*(trisk[1]+1)/(trisk[1]+s+1)
    }
    n.event[[i]] <- nevent
    n.cens[[i]]  <- ncens
    n.risk[[i]]  <- nrisk
  }
  
    temp <- list(n=n.used,
                 time = unlist(time),
                 n.risk = unlist(n.risk),
                 n.event= unlist(n.event),
                 n.censor = unlist(n.cens),
                 survUP = unlist(survUP),
                 survLOW = unlist(survLOW),
                 survLOW0 = unlist(survLOW0))
  
  if (se.fit) {
  
    zval <- qnorm(1- (1-conf.int)/2, 0,1)
    
    if (conf.type=='exact') {
      std.up <- sqrt(unlist(varhazUP))
      std.low <- sqrt(unlist(varhazLOW))
      temp$std.up <- std.up
      temp$std.low <- std.low
      temp$upper<-unlist(upper)
      temp$lower<-unlist(lower)
      temp$lower0<-unlist(lower0)
      temp$conf.type<-'exact'
      temp$conf.int<-conf.int
    }

    if (conf.type=='approx') {
      std.up <- sqrt(unlist(varhazUP)-unlist(survUP)^2)
      std.low <- sqrt(unlist(varhazLOW)-unlist(survLOW)^2)
      std.low0 <- sqrt(unlist(varhazLOW0)-unlist(survLOW0)^2)
      temp$std.up <- std.up
      temp$std.low <- std.low
      temp1 <- temp$survUP + zval* std.up 
      temp2 <- temp$survLOW - zval* std.low 
      temp3 <- temp$survLOW0 - zval* std.low0
      temp <- c(temp, list(upper=pmin(temp1,1), lower=pmax(temp2,0), lower0=pmax(temp3,0),
                           conf.type='approx', conf.int=conf.int))
    }
    
  }
  
  names(strata) <- xlev[sort(unique(x))]
  temp$strata <- strata
  
  
  temp
}


plot.isurvfit <- function(x,se.fit=TRUE,...) {
  ns <- as.numeric(x$strata)
  start <- 1
  for (i in 1:length(ns)) {
    keep <- start:(start+ns[i]-1)
    step <- ceiling(ns[i]/15)
    s.dots <- seq(1,ns[i],step)
    keep.dots <- keep[s.dots]
    s.dots.shift <- seq(ceiling(step/2),ns[i],step)
    keep.dots.shift <- keep[s.dots.shift]
    time <- c(0, x$time[keep])
    if (i==1) {
       do.call("plot",list(time, c(1,x$survUP[keep]), type = "s", lwd=2, col = i, xlab="t", ylab="S(t)",
           ylim=c(0, 1),...))
   } else{
      lines(time, c(1,x$survUP[keep]), type = "s", lwd=2, col = i, xlab="t", ylab="S(t)",
            ylim=c(0, 1))
    }
    lines(time[s.dots.shift+1], x$survUP[keep.dots.shift], type = "p", lwd=1, col = i,pch=i)
    lines(time, c(x$survLOW0[i],x$survLOW[keep]), type = "s", lwd=1, col = i)
    lines(time[s.dots+1], c(x$survLOW[keep.dots]), type = "p", lwd=1, col = i,pch=i)
    if (se.fit & !is.null(x$upper)) {
      lines(time, c(1, x$upper[keep]), type = "s", lwd=1, col = i, lty = 2)
      lines(time[s.dots+1], c(x$upper[keep.dots]), type = "p", lwd=1, col = i,pch=i)
      lines(time, c(x$lower0[i], x$lower[keep]), type = "s", lwd=1, col = i, lty = 2)
      lines(time[s.dots.shift+1], x$lower[keep.dots.shift], type = "p", lwd=1, col = i,pch=i)
    }
    start <- start+ns[i]
  }
}

print.isurvfit <- function(x,...) {
  ns <- as.numeric(x$strata)
  start <- 1
  n.events <- integer(length(ns))
  n.censored <- integer(length(ns))
  for (i in 1:length(ns)) {
    keep <- start:(start+ns[i]-1)
    n.events[i] <- sum(x$n.event[keep])
    n.censored[i] <- sum(x$n.censor[keep])
    start <- start+ns[i]
  }
  group.name <- names(x$strata)
  n.records <- x$strata
  df = data.frame(n.records, n.events, n.censored)
  print(df)
}

"[.isurvfit" <- function(x, ..., drop=TRUE) {
  nmatch <- function(indx, target) { 
    # This function lets R worry about character, negative, or logical subscripts
    #  It always returns a set of positive integer indices
    temp <- 1:length(target)
    names(temp) <- target
    temp[indx]
  }
  
  if (missing(..1)) i<- NULL  else i <- ..1
  if (missing(..2)) j<- NULL  else j <- ..2
  if (is.null(i) && is.null(j)) return (x) #no subscripts present!
  if (!is.matrix(x$surv) && !is.null(j))
    stop("survfit object does not have 2 dimensions")
  
  if (is.null(x$strata)) {
    if (is.matrix(x$surv)) {
      if (is.null(j) && !is.null(i)) j <- i #special case noted above
      x$surv <- x$surv[,j,drop=drop]
      if (!is.null(x$std.err)) x$std.err <- x$std.err[,j,drop=drop]
      if (!is.null(x$upper)) x$upper <- x$upper[,j,drop=drop]
      if (!is.null(x$lower)) x$lower <- x$lower[,j,drop=drop]
      if (!is.null(x$cumhaz)) x$cumhaz <- x$cumhaz[,j,drop=drop]
    }
    else warning("survfit object has only a single survival curve")
  }
  else {
    if (is.null(i)) keep <- seq(along.with=x$time)
    else {
      indx <- nmatch(i, names(x$strata)) #strata to keep
      if (any(is.na(indx))) 
        stop(paste("strata", 
                   paste(i[is.na(indx)], collapse=' '),
                   'not matched'))
      
      # Now, indx may not be in order: some can use curve[3:2] to reorder
      #  The list/unlist construct will reorder the data
      temp <- rep(1:length(x$strata), x$strata)
      keep <- unlist(lapply(indx, function(x) which(temp==x)))
      
      if (length(indx) <=1 && drop) x$strata <- NULL
      else               x$strata  <- x$strata[i]
      
      x$n       <- x$n[indx]
      x$time    <- x$time[keep]
      x$n.risk  <- x$n.risk[keep]
      x$n.event <- x$n.event[keep]
      x$n.censor<- x$n.censor[keep]
      if (!is.null(x$enter)) x$enter <- x$enter[keep]
    }
    if (is.matrix(x$surv)) {
      # If the curve has been selected by strata and keep has only
      #  one row, we don't want to lose the second subscript too
      if (!is.null(i) && (is.null(j) ||length(j) >1)) drop <- FALSE
      if (is.null(j)) {
        x$surv <- x$surv[keep,,drop=drop]
        if (!is.null(x$std.err)) 
          x$std.err <- x$std.err[keep,,drop=drop]
        if (!is.null(x$upper)) x$upper <-x$upper[keep,,drop=drop]
        if (!is.null(x$lower)) x$lower <-x$lower[keep,,drop=drop]
        if (!is.null(x$cumhaz)) x$cumhaz <-x$cumhaz[keep,,drop=drop]
      }
      else {
        x$surv <- x$surv[keep,j, drop=drop]
        if (!is.null(x$std.err)) 
          x$std.err <- x$std.err[keep,j, drop=drop]
        if (!is.null(x$upper)) x$upper <- x$upper[keep,j, drop=drop]
        if (!is.null(x$lower)) x$lower <- x$lower[keep,j, drop=drop]
        if (!is.null(x$cumhaz)) x$cumhaz <- x$cumhaz[keep,j, drop=drop]
      }
    }
    else {
      x$surv <- x$surv[keep]
      if (!is.null(x$std.err)) x$std.err <- x$std.err[keep]
      if (!is.null(x$upper)) x$upper <- x$upper[keep]
      if (!is.null(x$lower)) x$lower <- x$lower[keep]
      if (!is.null(x$cumhaz)) x$cumhaz <- x$cumhaz[keep]
    }
  }
  x
}

