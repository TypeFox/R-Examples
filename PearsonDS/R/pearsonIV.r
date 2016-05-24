# to do: permit vector-valued parameters!
# to do: improve log.p, lower.tail!
dpearsonIV <-
function(x,m,nu,location,scale,params,log=FALSE) {
  if (!missing(params)) { m <- params[[1]]; nu <- params[[2]]; 
                          location <- params[[3]]; scale <- params[[4]] }
  if (max(length(m),length(nu),length(location),length(scale))>1)
    stop("vector-valued parameters not (yet) allowed")                        
  if ((scale>0)&&(m>0.5)) {
#  stopifnot((scale>0)&&(m>0.5))
    .hasGSL <- suppressWarnings(suppressPackageStartupMessages(require(gsl)))
    logspace <- TRUE                                                              # make this an input parameter?
    if (logspace) {
      if (.hasGSL) {
        k <- -0.5*log(pi)-log(scale)-lgamma(m-0.5)+2*
              Re(gsl::lngamma_complex(m+nu/2*1i))-lgamma(m)
      } else {
        k <- .Call("logPearsonIVnorm",m,nu,scale,package="PearsonDS")
      }  
      if (log) {
        k-m*log(1+((x-location)/scale)^2)-nu*atan((x-location)/scale)
      } else {
        exp(k-m*log(1+((x-location)/scale)^2)-nu*atan((x-location)/scale))
      }  
    } else {
      if (.hasGSL) {
        k <- exp(-0.5*log(pi)-log(scale)-lgamma(m-0.5)+2*                         # improve?
              Re(gsl::lngamma_complex(m+nu/2*1i))-lgamma(m))
      } else {
        k <- .Call("PearsonIVnorm",m,nu,scale,package="PearsonDS")
      }  
      if (log) {
        log(k)-m*log(1+((x-location)/scale)^2)-nu*atan((x-location)/scale)
      } else {
        k*(1+((x-location)/scale)^2)^(-m)*exp(-nu*atan((x-location)/scale))
      }  
    }
  } else {
    tmp <- rep(NaN,length(x))
    tmp[is.na(x)&(!is.nan(x))] <- NA
    tmp
  }  
}

rpearsonIV <-
function(n,m,nu,location,scale,params) {
  if (!missing(params)) { m <- params[[1]]; nu <- params[[2]]; 
                          location <- params[[3]]; scale <- params[[4]] }
#  stopifnot((scale>0)&&(m>0.5))
  .hasGSL <- suppressWarnings(suppressPackageStartupMessages(require(gsl)))
  logspace <- TRUE                                                              # make this an input parameter?
  if (max(length(m),length(nu),length(location),length(scale))>1)
    stop("vector-valued parameters not (yet) allowed")                        
  if (m<1) stop("simulation for m<1 currently not available")
  if (logspace) {
    if (.hasGSL) {
      k <- -0.5*log(pi)-#log(scale)-                                            # bug, reported by Dave DeMers on 2013-11-19, corrected in v. 0.97, 2013-11-20
            lgamma(m-0.5)+2*Re(gsl::lngamma_complex(m+nu/2*1i))-lgamma(m)
    } else {
      k <- .Call("logPearsonIVnorm",m,nu,1.0,package="PearsonDS")               # bug, reported by Dave DeMers on 2013-11-19, corrected in v. 0.97, 2013-11-20
    }  
    .Call("rPearsonIVlogK",n,m,nu,scale,location,k,package="PearsonDS")
  } else {
    if (.hasGSL) {
      k <- -0.5*log(pi)-#log(scale)-                                            # bug, reported by Dave DeMers on 2013-11-19, corrected in v. 0.97, 2013-11-20
            lgamma(m-0.5)+2*Re(gsl::lngamma_complex(m+nu/2*1i))-lgamma(m)
    } else {
      k <- .Call("logPearsonIVnorm",m,nu,1.0,package="PearsonDS")               # bug, reported by Dave DeMers on 2013-11-19, corrected in v. 0.97, 2013-11-20
    }  
    .Call("rPearsonIVk",n,m,nu,scale,location,exp(k),package="PearsonDS")
  }  
}

F21  <- function(aa,bb,cc,zz,tol=1e-8,minit,maxit,DEBUG=FALSE) {
  if (missing(minit)) minit <- max(pmax(0,-Re(c(aa,bb,cc))))+5
  if (missing(maxit)) maxit <- max(1000000,10*minit)
  if (DEBUG) cat("Minit:",minit,", Maxit:",maxit,"\n")
  if (tol<.Machine$double.eps) {
    tol <- .Machine$double.eps
    warning(paste("tol too small, using",tol,"instead"))
  }  
  suffix <- ""
  if ((aa==1)&&is.double(bb)) suffix <- "a1bR" else
    if ((aa==1)&&is.double(cc)) suffix <- "a1cR" else
      if (is.double(aa)) suffix <- "aR"
  Dfun <- paste0("F21D",suffix)
  DDfun <- paste0("F21DD",suffix)
  QDfun <- paste0("F21QD",suffix)
  Dres <- .Call(Dfun,aa,bb,cc,zz,minit,maxit,package="PearsonDS")
  res  <- Dres$value
  ind  <- (Dres$rel*tol)<.Machine$double.eps
  ind[is.na(ind)] <- TRUE
  if (DEBUG) cat(sum(ind)," summands inexact (double)\n")
  rtolDD <- 0; rtolQD <- 0
  if (sum(!ind)>0) rtolD <- .Machine$double.eps/min(Dres$rel[!ind]) else rtolD <- 0
  if (sum(ind)>0) {
    DDres <- .Call(DDfun,aa,bb,cc,zz[ind],minit,maxit,package="PearsonDS")
    res[ind] <- DDres$value
    ind2 <- (DDres$rel*tol)<2.465190e-32
    ind2[is.na(ind2)] <- TRUE
    if (DEBUG) cat(sum(ind2)," summands inexact (double-double)\n")
    if (sum(!ind2)>0) rtolDD <- 2.465190e-32/min(DDres$rel[!ind2]) else rtolDD <- 0
    if (sum(ind2)>0) {
      QDres   <- .Call(QDfun,aa,bb,cc,zz[ind[ind2]],minit,maxit,package="PearsonDS")
      rtolQD  <- 3.038582e-64/min(QDres$rel)
      res[ind[ind2]] <- QDres$value
    }  
  }
  rtol <- max(.Machine$double.eps,rtolD,rtolDD,rtolQD)
  if (DEBUG) cat("actual tolerance: ",rtol,"\n")
  if (rtol>tol) warning(paste("actual tolerance",rtol,"bigger than target tolerance",tol))
  res
}

F21_3_5 <- function(aa,bb,cc,zz,...) {                                          # transformation 15.3.5 of 
  (1-zz)^(-bb)*                                                                 # Abramowitz&Stegun
  F21(bb,cc-aa,cc,zz/(zz-1),...)
}          

.ppearsonIVf21 <-
function(q,m,nu,location,scale,params,lower.tail=TRUE,log.p=FALSE,tol=1e-8,...) {
  if (!missing(params)) { m <- params[[1]]; nu <- params[[2]]; 
                          location <- params[[3]]; scale <- params[[4]] }
#  stopifnot((scale>0)&&(m>0.5))
  if (max(length(m),length(nu),length(location),length(scale))>1)
    stop("vector-valued parameters not (yet) allowed")
  if (TRUE) {
    x     <- (q-location)/scale
    ind1  <- x< -sqrt(3)
    ind2  <- x> sqrt(3)
    ind3  <- (x<0)&(!ind1)
    ind4  <- !(ind1|ind2|ind3)
    res   <- numeric(length(x))

    res[ind1] <- dpearsonIV(q[ind1],params=c(m,nu,location,scale)) *
                 scale/(2*m-1) * (1i - x[ind1]) *
                 F21(1,m+nu/2*1i,2*m,2/(1-1i*x[ind1]),tol,...)

    res[ind2] <- 1-dpearsonIV(-q[ind2],params=c(m,-nu,-location,scale)) *
                 scale/(2*m-1) * (1i + x[ind2]) *
                 F21(1,m-nu/2*1i,2*m,2/(1+1i*x[ind2]),tol,...)
    if (abs(nu)<(4-2*sqrt(3))*m) {
      res[ind3] <- (1-exp(-(nu+1i*2*m)*pi))^(-1) -
                   (1i*scale*dpearsonIV(q[ind3],params=c(m,nu,location,scale)))/
                   (1i*nu-2*m+2) * (1+x[ind3]^2) *
                   F21_3_5(1,2-2*m,2-m+1i*nu/2,(1+1i*x[ind3])/2,
                                       tol,...)
      res[ind4] <- 1 - (1-exp(-(-nu+1i*2*m)*pi))^(-1) +
                   (1i*scale*dpearsonIV(-q[ind4],params=c(m,-nu,-location,scale)
                   ))/(1i*(-nu)-2*m+2) * (1+(-x[ind4])^2) *
                   F21_3_5(1,2-2*m,2-m-1i*nu/2,(1-1i*x[ind4])/2,
                                       tol,...)
    } else {
      res[ind3] <- (1-exp(-(nu+1i*2*m)*pi))^(-1) -
                   (1i*scale*dpearsonIV(q[ind3],params=c(m,nu,location,scale)))/
                   (1i*nu-2*m+2) * (1+x[ind3]^2) *
                   F21(1,2-2*m,2-m+1i*nu/2,(1+1i*x[ind3])/2,tol,...)
      res[ind4] <- 1 - (1-exp(-(-nu+1i*2*m)*pi))^(-1) +
                   (1i*scale*dpearsonIV(-q[ind4],params=c(m,-nu,-location,scale)
                   ))/(1i*(-nu)-2*m+2) * (1+(-x[ind4])^2) *
                   F21(1,2-2*m,2-m-1i*nu/2,(1-1i*x[ind4])/2,tol,...)
    }
    res <- Re(res)
  }
}

.ppearsonIVint <-
function(q,m,nu,location,scale,params,lower.tail=TRUE,log.p=FALSE,tol=1e-8,...) {
  if (!missing(params)) { m <- params[[1]]; nu <- params[[2]]; 
                          location <- params[[3]]; scale <- params[[4]] }
#  stopifnot((scale>0)&&(m>0.5))
  if (max(length(m),length(nu),length(location),length(scale))>1)
    stop("vector-valued parameters not (yet) allowed")
  .hasGSL <- suppressWarnings(suppressPackageStartupMessages(require(gsl)))
  modus <- location - scale*nu/(2*m)
  res   <- numeric(length(q))
  ind   <- q>modus
  if (.hasGSL) {
    k <- -0.5*log(pi)-log(scale)-lgamma(m-0.5)+2*
          Re(gsl::lngamma_complex(m+nu/2*1i))-lgamma(m)
  } else {
    k <- .Call("logPearsonIVnorm",m,nu,scale,package="PearsonDS")
  }
  intfn <- function(x) exp(k-m*log(1+((x-location)/scale)^2)-nu*
                       atan((x-location)/scale))
  if (sum(ind)>0) res[ind] <-
    sapply(q[ind],function(x) 1-integrate(function(z) intfn(z),
                                lower=x,upper=+Inf,rel.tol=tol)$value)
  if (sum(!ind)>0) res[!ind] <-
    sapply(q[!ind],function(x) integrate(function(z) intfn(z),
                               lower=-Inf,upper=x,rel.tol=tol)$value)
  if (!lower.tail) res <- 1-res
  if (log.p)       res <- log(res)
  res
}

.ppearsonIVold <-
function(q,m,nu,location,scale,params,DEBUG=FALSE,tol=1e-8) {
  if (!missing(params)) { m <- params[[1]]; nu <- params[[2]]; 
                          location <- params[[3]]; scale <- params[[4]] }
#  stopifnot((scale>0)&&(m>0.5))
  x   <- (q-location)/scale
  phi  <- pi/2 + atan(x)
  u    <- 1-m+nu/2*1i
  norm <- exp(gsl::lngamma_complex(u+1)-gsl::lngamma_complex(u+2*m-1)+lgamma(2*m-1))
  Re((exp(-nu*phi)*exp((2-2*m)*phi*1i)*
     (F21(2-2*m,u,u+1,exp(phi*2i),DEBUG=DEBUG,tol=tol)/
     norm)-1)/
     (exp(-nu*pi)*exp((2-2*m)*pi*1i)-1))
}

ppearsonIV <- 
function(q,m,nu,location,scale,params,lower.tail=TRUE,log.p=FALSE,tol=1e-8,...) {                                              
  if (!missing(params)) { m <- params[[1]]; nu <- params[[2]]; 
                          location <- params[[3]]; scale <- params[[4]] }
  if (missing(params)) params <- c(m,nu,location,scale)
#  stopifnot((scale>0)&&(m>0.5))
  if (max(length(m),length(nu),length(location),length(scale))>1)
    stop("vector-valued parameters not (yet) allowed")      
  if ((scale>0)&&(m>0.5)) {
    .hasGSL <- suppressWarnings(suppressPackageStartupMessages(require(gsl)))
    res <- numeric(length(q))
    res[is.na(q)] <- q[is.na(q)]                    
    if (!all(is.na(q))) {
      if (isTRUE(any(m<8,m>156,!.hasGSL))) {
        res[!is.na(q)] <- .ppearsonIVint(q[!is.na(q)],params=params,tol=tol,...)
      } else if (isTRUE(all(m>28+0.72*nu,m> -150+3.75*nu,m<5+1.85*nu))) {
        res[!is.na(q)] <- .ppearsonIVold(q[!is.na(q)],params=params,tol=tol,...)
      } else if (isTRUE(m>28+0.72*nu)) {
        res[!is.na(q)] <- .ppearsonIVint(q[!is.na(q)],params=params,tol=tol,...)
      } else {
        res[!is.na(q)] <- .ppearsonIVf21(q[!is.na(q)],params=params,tol=tol,...)
      }
    }  
    if (!lower.tail) res <- 1-res
    if (log.p)       res <- log(res)
    res
  } else {
    tmp <- rep(NaN,length(q))
    tmp[is.na(q)&(!is.nan(q))] <- NA
    tmp
  }  
}

qpearsonIV <-
function(p,m,nu,location,scale,params,lower.tail=TRUE,log.p=FALSE,tol=1e-8,...) {
  if (!missing(params)) { m <- params[[1]]; nu <- params[[2]]; 
                          location <- params[[3]]; scale <- params[[4]] }
#  stopifnot((scale>0)&&(m>0.5))
  if (max(length(m),length(nu),length(location),length(scale))>1)
    stop("vector-valued parameters not (yet) allowed")          
  if ((scale>0)&&(m>0.5)) {
    if (log.p)       p <- exp(p)
    if (!lower.tail) p <- 1-p                
    ind   <- (0<p)&(p<1)&(!is.na(p))
    modus <- location - scale*nu/(2*m)
    n     <- length(p)
    res   <- numeric(n)
    res[is.na(p)] <- p[is.na(p)]
    res[p<0]      <- NaN
    res[p>1]      <- NaN
    res[p==0]     <- -Inf
    res[p==1]     <- Inf
    for (i in seq_len(n)[ind]) {
      xold  <- modus
      numit <- 0
      repeat {
        numit <- numit + 1
        dp    <- dpearsonIV(xold,m,nu,location,scale)
        ttol  <- max(tol*1e-2*dp,.Machine$double.eps)
        xnew  <- xold-(ppearsonIV(xold,m,nu,location,scale,tol=ttol,...)-p[i])/dp
        if ((abs(xnew-xold)<tol)||(numit>30)) break                
        xold <- xnew
      }
      res[i] <- xnew
    }
    res
  } else {
    tmp <- rep(NaN,length(p))
    tmp[is.na(p)&(!is.nan(p))] <- NA
    tmp
  }  
}
