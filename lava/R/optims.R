###{{{ nlminb

nlminb2 <- function(start,objective,gradient,hessian,...) {
  nlminbcontrols <- c("eval.max","iter.max","trace","abs.tol","rel.tol","x.tol","step.min")
  dots <- list(...)
  control <- list(...)$control
  control <- control[names(control)%in%nlminbcontrols]
  dots$control <- control
  if (length(dots$trace)>0 && dots$trace>0) cat("\n")
  mypar <- c(list(start=start,objective=objective,gradient=gradient,hessian=hessian),dots)
  mypar["debug"] <- NULL
  do.call("nlminb", mypar)
}

nlminb1 <- function(start,objective,gradient,hessian,...) {
  nlminb2(start,objective,gradient=gradient,hessian=NULL,...)
}
nlminb0 <- function(start,objective,gradient,hessian,...) {
  nlminb2(start,objective,gradient=NULL,hessian=NULL,...)
}

###}}} nlminb

###{{{ estfun

estfun <- function(start,objective,gradient,hessian,NR=FALSE,...) {
  myobj <- function(x,...) {
    S <- gradient(x,...)
    crossprod(S)[1]
  }
  if (!missing(hessian) && !is.null(hessian)) {
    mygrad <- function(x) {
      H <- hessian(x)
      S <- gradient(x)
      2*S%*%H
    }
  } else {
      hessian <- function(x) numDeriv::jacobian(gradient,x,method=lava.options()$Dmethod)
      mygrad <- function(x) {
          H <- hessian(x)
          S <- gradient(x)
          2*S%*%H
      }
  }
  if (NR) {
      op <- lava::NR(start,gradient=gradient,hessian=hessian,...)
  } else {
      op <- nlminb2(start,myobj,mygrad,hessian=NULL,...)
  }
  return(op)
}

estfun0 <- function(...,hessian=NULL) estfun(...,hessian=hessian)

###}}}

###{{{ Newton-Raphson/Scoring

##' @export
NR <- function(start,objective,gradient,hessian,debug=FALSE,control,...) {
    control0 <- list(trace=0,gamma=1,lambda=0,ngamma=0,gamma2=0,backtrace=TRUE,
                     iter.max=200,tol=1e-9,stabil=FALSE,epsilon=1e-9)
    if (!missing(control)) {
        control0[names(control)] <- control
    }

    if (control0$trace>0)
        cat("\nIter=0;\t\n",
            "\tp=", paste0(formatC(start), collapse=" "),"\n")

    gradFun = !missing(gradient)
    if (!gradFun & missing(hessian)) {
        hessian <- function(p) {
            ff <- objective(p)
            res <- attributes(ff)$hessian
            attributes(res)$grad <- as.vector(attributes(ff)$grad)
            return(res)
        }
    }
    oneiter <- function(p.orig,Dprev,return.mat=FALSE) {
        if (is.null(hessian)) {
            cat(".")
            I <- -numDeriv::jacobian(gradient,p.orig,method=lava.options()$Dmethod)
        } else {
            I <- -hessian(p.orig)
        }
        D <- attributes(I)$grad
        if (is.null(D)) {
            D <- gradient(p.orig)
        }
        if (return.mat) return(list(D=D,I=I))
        if (control0$stabil) {
            if (control0$lambda!=0) {
                if (control0$lambda<0) {
                    sigma <- (t(D)%*%(D))[1]
                } else {
                    sigma <- control0$lambda
                }
                sigma <- min(sigma,10)
                I <- I+control0$gamma2*sigma*diag(nrow=nrow(I))
            } else {
                sigma <- ((D)%*%t(D))
                I <- I+control0$gamma2*(sigma)
            }
        }

        svdI <- svd(I); svdI$d0 <- numeric(length(svdI$d));
        svdI$d0[abs(svdI$d)>control0$epsilon] <-
            1/svdI$d[abs(svdI$d)>control0$epsilon]
        iI <- with(svdI,  (v)%*%diag(d0,nrow=length(d0))%*%t(u))
        Delta = control0$gamma*iI%*%D

        Lambda <- 1
        if (control0$backtrace) {
            mD0 <- mean(Dprev^2)
            mD <- mean(D^2)
            while (mD>=mD0) {
                p <- p.orig + Lambda*Delta
                if (gradFun) {
                    D = gradient(p)
                } else {
                    DI <- oneiter(p,return.mat=TRUE)
                    D = DI$D
                }
                mD = mean(D^2)
                if (is.nan(mD)) mD=mD0
                Lambda <- Lambda/2
            if (Lambda<1e-12) break;
            }
        }
        p <- p.orig + Lambda*Delta
        Lambda <- Lambda^.5
        return(list(p=p,D=D,iI=iI))
    }

    count <- count2 <- 0
    thetacur <- start
    gammacount <- 0
    Dprev <- rep(Inf,length(start))
    for (jj in seq_len(control0$iter.max)) {
        gammacount <- gammacount+1
        count <-  count+1
        count2 <- count2+1
        oldpar <- thetacur
        newpar <- oneiter(thetacur,Dprev)
        Dprev <- newpar$D
        thetacur <- newpar$p
        if (!is.null(control0$ngamma) && control0$ngamma>0) {
            if (control0$ngamma<=gammacount) {
                control0$gamma <- sqrt(control0$gamma)
                gammacount <- 0
            }
        }
        if (count2==control0$trace) {
            cat("Iter=",count, ";\n\tD=", paste0(formatC(newpar$D), collapse=" "),"\n")
            cat("\tp=", paste0(formatC(thetacur), collapse=" "),"\n")
            count2 <- 0
        }
        if (mean(newpar$D^2)<control0$tol) break;
    }
    res <- list(par=as.vector(thetacur), iterations=count, method="NR",
                gradient=newpar$D, iH=newpar$iI)
    return(res)
}

###}}} Newton Raphson/Scoring
