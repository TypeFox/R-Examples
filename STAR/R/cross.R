crossTight <- function(tMax=1,
                       h=0.001,
                       a=0.3,
                       b=2.35,
                       withBounds=TRUE,
                       logScale=FALSE) {

  n <- round(tMax/h)
  tMax <- n*h

  if (!withBounds && !logScale) {
    result <- .C(crossTight_sym,
                 as.double(tMax),
                 as.integer(n),
                 as.double(a),
                 as.double(b),
                 G=double(n)
                 )$G
    result <- list(G=result)
  }
  if (!withBounds && logScale) {
    result <- .C(crossTightlog,
                 as.double(tMax),
                 as.integer(n),
                 as.double(a),
                 as.double(b),
                 G=double(n)
                 )$G
    result <- list(G=result)
  }
  if (withBounds && !logScale) {
    result <- .C(crossTightWithB,
                 as.double(tMax),
                 as.integer(n),
                 as.double(a),
                 as.double(b),
                 G=double(n),
                 Gu=double(n),
                 Gl=double(n)
                 )[c("Gu","G","Gl")]
  }
  if (withBounds && logScale) {
    result <- .C(crossTightWithBlog,
                 as.double(tMax),
                 as.integer(n),
                 as.double(a),
                 as.double(b),
                 G=double(n),
                 Gu=double(n),
                 Gl=double(n)
                 )[c("Gu","G","Gl")]
  }

  result$time <- (1:n)*h
  result$g <- diff(c(0,result$G))/h
  result$mids <- (1:n)*h-0.5*h
  result$h <- h
  result$call <- match.call()
  class(result) <- "FirstPassageTime"
  result
  
}

mkTightBMtargetFct <- function(ci=0.95,
                               tMax=1,
                               h=0.001,
                               logScale=FALSE) {

  ## Check ci
  if (!(0<ci && ci <1)) stop("Expect 0 < ci < 1")
  trueTarget <- (1-ci)/2
  n <- round(tMax/h)
  function(logp) {
    p <- exp(logp)
    (trueTarget - crossTight(tMax,h,p[1],p[2],FALSE,logScale)$G[n])^2
  }
  
}

crossGeneral <- function(tMax=1,
                         h=0.001,
                         cFct,
                         cprimeFct,
                         bFct,
                         withBounds=FALSE,
                         Lplus
                         ) {

  ## Check arg. cFct. It should be a function or a numeric
  if (!is.function(cFct)) {
    if (!is.numeric(cFct)) {
      stop("cFct should be a function or a numeric.")
    } else {
      cCst <- cFct[1]
      cFct <- function(t) rep(cCst,length(t))
    } ## End of conditional on !is.numeric(cFct) 
  } ## End of conditional on !is.function(cFct) 

  n <- round(tMax/h)
  tMax <- n*h
  tt <- (1:n)*h
  
  ## Check arg. bFct.
  if (!missing(bFct)) {
    if (!is.function(bFct)) stop("bFct should be a function.")
    if (withBounds) {
      if (any(bFct(tt) < 0)) stop("bFct should be >= 0.")
      if (missing(Lplus) && !missing(cprimeFct)) {
        ## Do checks on bFct
        if (!is.function(cprimeFct)) stop("cprimeFct should be a function.")
        testFct <- function(u,t) 2*cprimeFct(u) - (cFct(t)-cFct(u))/(t-u)
        tests <- sapply(tt[-1],
                        function(t) {
                          b.t <- bFct(t)
                          u <- tt[tt<t]
                          tu <- testFct(u,t)
                          ub <- min(tu)
                          lb <- max(0,max(tu))
                          c(test3p9 = (b.t > ub),
                            test3p12 = (lb > b.t)
                            )
                        }
                        )
        tests <- apply(tests,1,any)
        if (tests[1]) {
          Lplus <- TRUE
        } else {
          if (tests[2]) {
            Lplus <- FALSE
          } else {
            stop("bFct not suitable for bounds calculation.")
          }
        } ## End of conditional of tests[1] 
      } ## End of conditional on !missing(cprimeFct)
    } ## End of conditional on withBounds
  } else {
    if (withBounds && missing(Lplus)) stop("Lplus required")
  } ## End of conditional on !missing(bFct)

  if (missing(bFct)) {
    fFct <- function(t) 2*pnorm(-cFct(t)/sqrt(t))
    kFct <- function(t,u) {
      numerator <- numeric(length(u))
      denominator <- numerator
      ratio <- numerator
      result <- numerator
      largeU <- u >= t
      numerator[!largeU] <- cFct(u[!largeU])-cFct(t)
      denominator[!largeU] <- sqrt(t-u[!largeU])
      ratio[!largeU] <- numerator[!largeU]/denominator[!largeU]
      ratio[!largeU][is.nan(ratio[!largeU])] <- 0
      result[!largeU] <- 2*pnorm(ratio[!largeU])
      result[u==t] <- 1
      result
    }
  } else {
    fFct <- function(t) {
      term1 <- pnorm(-cFct(t)/sqrt(t))
      term2 <- exp(-2*bFct(t)*(cFct(t)-t*bFct(t)))*
        pnorm((-cFct(t)+2*t*bFct(t))/sqrt(t))
      term1+term2
    }
    kFct <- function(t,u) {
      numerator1 <- numeric(length(u))
      numerator2 <- numerator1
      denominator <- numerator1
      ratio1 <- numerator1
      ratio2 <- numerator1
      term1 <- numerator1
      term2 <- numerator1
      largeU <- u >= t
      numerator1[!largeU] <- cFct(u[!largeU])-cFct(t)
      denominator[!largeU] <- sqrt(t-u[!largeU])
      ratio1[!largeU] <- numerator1[!largeU]/denominator[!largeU]
      ## ratio1[!largeU][is.nan(ratio1[!largeU])] <- 0
      term1[!largeU] <- pnorm(ratio1[!largeU])
      term1[u == t] <- 0.5
      numerator2[!largeU] <- cFct(u[!largeU])-cFct(t)+2*(t-u[!largeU])*bFct(t)
      ratio2[!largeU] <- numerator2[!largeU]/denominator[!largeU]
      ## ratio2[!largeU][is.nan(ratio2[!largeU])] <- 0
      term2[!largeU] <- exp(-2*bFct(t)*
                            (cFct(t)-cFct(u[!largeU])-
                             (t-u[!largeU])*bFct(t)))*pnorm(ratio2[!largeU])
      term2[u == t] <- 0.5
      term1+term2
    }
  } ## End of conditional on missing(bFct)

  ## We can now implement the mid-point method of 
  ## Loader & Deely (1987) pp 99-100
  
  mt <- tt-0.5*h
  B <- fFct(tt)
  A <- t(sapply(tt,function(x) kFct(x,mt)))
  x <- forwardsolve(A,B)
  x[x<0] <- 0
  
  result <- list(time=tt,
                 mids=mt,
                 G=cumsum(x),
                 g=x/h,
                 h=h
                 )

  class(result) <- "FirstPassageTime"

  if (withBounds) {
    n <- length(x)
    Gu <- numeric(n)
    Gl <- numeric(n)  
    if (Lplus) {
      ## First case for bounds considered by Loader and Deely
      ## Corresponding to Eq. 3.6 and 3.7
      f <- fFct(h)
      Gu[1] <- f/kFct(h,0)
      Gl[1] <- f

      for (i in 2:n) {
        f <- fFct(i*h)
        sumU <- f
        sumL <- f
        kjm <- kFct(i*h,0)
        kj <- kFct(i*h,h)
        for (j in 1:(i-1)) {
          if (j == i-1) {
            kjp <- 1
          } else {
            kjp <- kFct(i*h,(j+1)*h)
          }
          sumU <- sumU + (kj-kjm)*Gu[j]
          sumL <- sumL + (kjp-kj)*Gl[j]
          kjm <- kj
          kj <- kjp
        } ## End of the loop on j
        Gu[i] <- sumU/kjm
        Gl[i] <- sumL
      } ## End of the loop on i
    } else {
      ## Second case for bounds considered by Loader and Deely
      ## Corresponding to Eq. 3.10 and 3.11
      f <- fFct(h)
      Gu[1] <- f
      kjm <- kFct(h,0)
      kj <- 1
      Gl[1] <- f-Gu[1]*(kjm-kj)
      for (i in 2:n) {
        f <- fFct(i*h)
        sumU <- f
        sumL <- f
        kjm <- kFct(i*h,0)
        kj <- kFct(i*h,h)
        for (j in 1:(i-1)) {
          if (j == i-1) {
            kjp <- 1
          } else {
            kjp <- kFct(i*h,(j+1)*h)
          }
          sumU <- sumU - (kj-kjp)*Gl[j]
          sumL <- sumL - (kjm-kj)*Gu[j]
          kjm <- kj
          kj <- kjp
        } ## End of the loop on j
        Gu[i] <- sumU
        Gl[i] <- sumL
      } ## End of the loop on i
    } ## End of conditional on Lplus
    result$Gu <- Gu
    result$Gl <- Gl
  } ## End of conditional on withBounds
  
  result
  
}

plot.FirstPassageTime <- function(x,y,
                                  which=c("Distribution","density"),
                                  xlab,ylab,...) {
  if (missing(xlab)) xlab <- "Time (s)"
  if (missing(ylab)) {
    ylab <- switch(which[1],
                   Distribution="Probability",
                   density="Density (1/s)"
                   )
  }
  switch(which[1],
         Distribution=plot(x$time,x$G,xlab=xlab,ylab=ylab,type="l",...),
         density=plot(x$mids,x$g,xlab=xlab,ylab=ylab,type="l",...)
         )
  
}

lines.FirstPassageTime <- function(x,
                                   which=c("Distribution","density"),
                                   ...) {
  switch(which[1],
         Distribution=lines(x$time,x$G,...),
         density=lines(x$mids,x$g,...)
         )
  
}

summary.FirstPassageTime <- function(object,digits,...) {
  n <- length(object$time)
  tMax <- object$time[n]
  Gmax <- object$G[n]
  if ("Gu" %in% names(object)) {
    GuMax <- object$Gu[n]
    GlMax <- object$Gl[n]
    range <- GuMax - GlMax
    i <- 1
    rr <- range*10^i
    while (rr < 1) {
      i <- i+1
      rr <- range*10^i
    }
    cat(paste("Prob. of first passage before " ,tMax,": ",round(Gmax,digits=i),sep=""))
    cat(paste(" (bounds: [",round(GlMax,digits=i),",",round(GuMax,digits=i),"])",sep=""))
  } else {
    if (missing(digits)) digits <- 4
    cat(paste("Prob. of first passage before " ,tMax,": ",round(Gmax,digits=digits),sep=""))
  }
  cat("\n")
  cat(paste("Integration time step used: ",object$h,".\n",sep=""))
  
}

print.FirstPassageTime <- function(x,...) {

  n <- length(x$time)
  tMax <- x$time[n]
  Gmax <- x$G[n]
  result <- Gmax
  if ("Gu" %in% names(x)) {
    GuMax <- x$Gu[n]
    GlMax <- x$Gl[n]
    attr(result,"range") <- c(GlMax,GuMax)
  }

  print(result)
  
}

