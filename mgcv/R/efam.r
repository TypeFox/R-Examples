## (c) Simon N. Wood (ocat, tw, nb, ziP) & Natalya Pya (scat, beta), 
## 2013-2015. Released under GPL2.

## extended families for mgcv, standard components. 
## family - name of family character string
## link - name of link character string
## linkfun - the link function
## linkinv - the inverse link function
## mu.eta - d mu/d eta function (derivative of inverse link wrt eta)
## note: for standard links this information is supplemented using 
##       function fix.family.link.extended.family with functions 
##       gkg where k is 2,3 or 4 giving the kth derivative of the 
##       link over the first derivative of the link to the power k.
##       for non standard links these functions muct be supplied.
## dev.resids - function computing deviance residuals.
## Dd - function returning derivatives of deviance residuals w.r.t. mu and theta. 
## aic - function computing twice - log likelihood for 2df to be added to.
## initialize - expression to be evaluated in gam.fit4 and initial.spg 
##              to initialize mu or eta.
## preinitialize - optional expression evaluated in estimate.gam to 
##                 e.g. initialize theta parameters (see e.g. ocat)
## postproc - expression to evaluate in estimate.gam after fitting (see e.g. betar)
## ls - function to evaluated log saturated likelihood and derivatives w.r.t.
##      phi and theta for use in RE/ML optimization. If deviance used is just -2 log 
##      lik. can njust return zeroes. 
## validmu, valideta - functions used to test whether mu/eta are valid.      
## n.theta - number of theta parameters.
## no.r.sq - optional TRUE/FALSE indicating whether r^2 can be computed for family
## ini.theta - function for initializing theta.
## putTheta, getTheta - functions for storing and retriving theta values in function 
##                      environment.
## rd - optional function for simulating response data from fitted model.
## residuals - optional function for computing residuals.
## predict - optional function for predicting from model, called by predict.gam.
## family$data - optional list storing any family specific data for use, e.g. in predict
##               function.



## extended family object for ordered categorical

ocat <- function(theta=NULL,link="identity",R=NULL) {
## extended family object for ordered categorical model.
## one of theta and R must be supplied. length(theta) == R-2.
## weights are ignored. #! is stuff removed under re-definition of ls as 0
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("identity")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
            linktemp <- stats$name
    } else stop(linktemp, " link not available for ordered categorical family; available links are \"identity\"")
  }
  if (is.null(theta)&&is.null(R)) stop("Must supply theta or R to ocat")
  if (!is.null(theta)) R <- length(theta) + 2 ## number of catergories
  ## Theta <-  NULL;
  n.theta <- R-2
  ## NOTE: data based initialization is in preinitialize...
  if (!is.null(theta)&&sum(theta==0)==0) {
    if (sum(theta<0)) iniTheta <- log(abs(theta)) ## initial theta supplied
    else { 
      iniTheta <- log(theta) ## fixed theta supplied
      n.theta <- 0
    }
  } else iniTheta <- rep(-1,length=R-2) ## inital log theta value

  env <- new.env(parent = .GlobalEnv)
  assign(".Theta", iniTheta, envir = env)

  putTheta <-function(theta) assign(".Theta", theta,envir=environment(sys.function()))

  getTheta <-function(trans=FALSE) { 
    theta <- get(".Theta")
    if (trans) { ## transform to (finite) cut points...
      R = length(theta)+2
      alpha <- rep(0,R-1) ## the thresholds
      alpha[1] <- -1
      if (R > 2) { 
        ind <- 2:(R-1)
        alpha[ind] <- alpha[1] + cumsum(exp(theta))
      } 
      theta <- alpha
    }
    theta
  }
  postproc <- expression({
    object$family$family <- 
    paste("Ordered Categorical(",paste(round(object$family$getTheta(TRUE),2),collapse=","),")",sep="")
  })

  validmu <- function(mu) all(is.finite(mu))

  dev.resids <- function(y, mu, wt,theta=NULL) {
    #F <- function(x) {
    #  h <- ind <- x > 0; h[ind] <- 1/(exp(-x[ind]) + 1)
    #  x <- exp(x[!ind]); h[!ind] <- (x/(1+x))
    #  h 
    #}
    Fdiff <- function(a,b) {
      ## cancellation resistent F(b)-F(a), b>a
      h <- rep(1,length(b)); h[b>0] <- -1; eb <- exp(b*h)
      h <- h*0+1; h[a>0] <- -1; ea <- exp(a*h)
      ind <- b<0;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- bi/(1+bi) - ai/(1+ai)
      ind1 <- a>0;bi <- eb[ind1];ai <- ea[ind1]
      h[ind1] <- (ai-bi)/((ai+1)*(bi+1))
      ind <- !ind & !ind1;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- (1-ai*bi)/((bi+1)*(ai+1))
      h
    }
    if (is.null(theta)) theta <- get(".Theta")
    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    al1 <- alpha[y+1];al0 = alpha[y] ## cut points above and below y
    ## Compute sign for deviance residuals, based on sign of interval
    ## midpoint for each datum minus the fitted value of the latent 
    ## variable. This makes first and last categories like 0s and 1s in
    ## logistic regression....
    s <- sign((al1 + al0)/2-mu) ## sign for deviance residuals
    al1mu <- al1-mu;al0mu <- al0-mu
    ## f1 <- F(al1mu);f0 <- F(al0mu);
    ##f <- pmax(f1 - f0,.Machine$double.eps)
    f <- Fdiff(al0mu,al1mu)
    ##a1 <- f1^2 - f1;a0 <- f0^2 - f0; a <- a1 -a0
#!    al1al0 <- (al1-al0)/2;al0al1 <- (al0-al1)/2
    ##g1 <- F(al1al0);g0 <- F(al0al1)
    ##A <- pmax(g1 - g0,.Machine$double.eps)
#!    A <- Fdiff(al0al1,al1al0)
    rsd <- -2*log(f) #! 2*(log(A)-log(f))
    attr(rsd,"sign") <- s
    rsd
  } ## end of dev.resids

  Dd <- function(y, mu, theta, wt=NULL, level=0) {
  ## derivatives of the deviance...
   # F <- function(x) { ## e^(x)/(1+e^x) without overflow
   #   h <- ind <- x > 0; h[ind] <- 1/(exp(-x[ind]) + 1)
   #   x <- exp(x[!ind]); h[!ind] <- (x/(1+x))
   #   h 
   # }
    Fdiff <- function(a,b) {
      ## cancellation resistent F(b)-F(a), b>a 
      h <- rep(1,length(b)); h[b>0] <- -1; eb <- exp(b*h)
      h <- h*0+1; h[a>0] <- -1; ea <- exp(a*h)
      ind <- b<0;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- bi/(1+bi) - ai/(1+ai)
      ind1 <- a>0;bi <- eb[ind1];ai <- ea[ind1]
      h[ind1] <- (ai-bi)/((ai+1)*(bi+1))
      ind <- !ind & !ind1;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- (1-ai*bi)/((bi+1)*(ai+1))
      h
    }
    abcd <- function(x,level=-1) {
      bj <- cj <- dj <- NULL
      ## compute f_j^2 - f_j without cancellation error
      ## x <- 10;F(x)^2-F(x);abcd(x)$aj
      h <- rep(1,length(x)); h[x>0] <- -1; ex <- exp(x*h)
      ex1 <- ex+1;ex1k <- ex1^2
      aj <- -ex/ex1k
      if (level>=0) {
        ## compute f_j - 3 f_j^2 + 2f_j^3 without cancellation error
        ##  x <- 10;F(x)-3*F(x)^2+2*F(x)^3;abcd(x,0)$bj
        ex1k <- ex1k*ex1;ex2 <- ex^2
        bj <- h*(ex-ex^2)/ex1k
        if (level>0) {
          ## compute -f_j + 7 f_j^2 - 12 f_j^3 + 6 f_j^4
          ## x <- 10;-F(x)+7*F(x)^2-12*F(x)^3+6*F(x)^4;abcd(x,1)$cj
          ex1k <- ex1k*ex1;ex3 <- ex2*ex
          cj <- (-ex3 + 4*ex2 - ex)/ex1k
          if (level>1) {    
            ## compute d_j
            ## x <- 10;F(x)-15*F(x)^2+50*F(x)^3-60*F(x)^4+24*F(x)^5;abcd(x,2)$dj
            ex1k <- ex1k*ex1;ex4 <- ex3*ex
            dj <- h * (-ex4 + 11*ex3 - 11*ex2 + ex)/ex1k
          }
        }
      }        
      list(aj=aj,bj=bj,cj=cj,dj=dj)
    }
 
    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    al1 <- alpha[y+1];al0 = alpha[y]
    al1mu <- al1-mu;al0mu <- al0 - mu
    ##f1 <- F(al1mu);f0 <- F(al0mu);
    ##f <- pmax(f1 - f0,.Machine$double.eps)
    f <- pmax(Fdiff(al0mu,al1mu),.Machine$double.xmin)
    r1 <- abcd(al1mu,level); a1 <- r1$aj
    r0 <- abcd(al0mu,level); a0 <- r0$aj
    ## a1 <- f1^2 - f1;a0 <- f0^2 - f0; 
    a <- a1 - a0
    
    #!al1al0 <- (al1-al0)/2; #! al0al1 <- (al0-al1)/2
    ##g1 <- F(al1al0);g0 <- F(al0al1)
    ##A <- pmax(g1 - g0,.Machine$double.eps)
#!    A <- Fdiff(al0al1,al1al0)
    if (level>=0) {
      ## b1 <- f1 - 3 * f1^2 + 2 * f1^3;b0 <- f0 - 3 * f0^2 + 2 * f0^3  
      b1 <- r1$bj;b0 <- r0$bj
      b <- b1 - b0
    }
    if (level>0) {
      ##c1 <- -f1 + 7 * f1^2 - 12* f1^3 + 6 * f1^4
      ##c0 <- -f0 + 7 * f0^2 - 12* f0^3 + 6 * f0^4
      c1 <- r1$cj;c0 <- r0$cj
      c <- c1 - c0
#!      R1 <- abcd(al1al0,level-2) 
#!      R0 <- abcd(al0al1,level-2)
      ## B <- g1^2 - g1 + g0^2 - g0
#!     B <- R1$aj + R0$aj
    }
    if (level>1) {
      ##d1 <- f1 - 15 * f1^2 + 50 * f1^3 - 60 * f1^4 + 24 * f1^5
      ##d0 <- f0 - 15 * f0^2 + 50 * f0^3 - 60 * f0^4 + 24 * f0^5
      d1 <- r1$dj;d0 <- r0$dj
      d <- d1 - d0
      ##C <- 2 * g1^3 - 3 * g1^2 + g1 - 2 * g0^3 + 3 * g0^2 - g0
#!      C <- R1$bj - R0$bj
    }

    oo <- list(D=NULL,Dmu=NULL,Dmu2=NULL,Dth=NULL,Dmuth=NULL,
             Dmu2th=NULL)
    n <- length(y)
    ## deviance...
    oo$D <- -2 * log(f) #! 2*(log(A)-log(f))
    if (level >= 0) { ## get derivatives used in coefficient estimation
      oo$Dmu <- -2 * a / f
      a2 <- a^2
      oo$EDmu2 <- oo$Dmu2 <- 2 * (a2/f -  b)/f
    }
    if (R<3) level <- 0 ## no free parameters

    if (level > 0) { ## get first derivative related stuff
      f2 <- f^2;a3 <- a2*a
      oo$Dmu3 <- 2*(- c - 2 * a3/f2 + 3 * a * b/f)/f
      Dmua0 <- 2 * (a0 * a/f -  b0)/f
      Dmua1 <- -2 * (a1 * a /f - b1)/f
      Dmu2a0 <- -2* (c0 + (a0*(2*a2/f - b)- 2*b0*a  )/f)/f
      Dmu2a1 <- 2*(c1  + (2*(a1*a2/f - b1*a) - a1*b)/f)/f
      Da0 <- -2*a0/f #! + B/A; 
      Da1 <-  2 * a1/f #! - B/A 
      ## now transform to derivatives w.r.t. theta...
      oo$Dmu2th <- oo$Dmuth <- oo$Dth <- matrix(0,n,R-2)
      for (k in 1:(R-2)) { 
        etk <- exp(theta[k])
        ind <- y == k+1
        oo$Dth[ind,k] <- Da1[ind]*etk
        oo$Dmuth[ind,k] <- Dmua1[ind]*etk
        oo$Dmu2th[ind,k] <- Dmu2a1[ind]*etk
    
        if (R>k+2) { 
          ind <- y > k+1 & y < R
          oo$Dth[ind,k] <- (Da1[ind]+Da0[ind])*etk
          oo$Dmuth[ind,k] <- (Dmua1[ind]+Dmua0[ind])*etk
          oo$Dmu2th[ind,k] <- (Dmu2a1[ind]+Dmu2a0[ind])*etk
       
        }
        ind <- y == R
        oo$Dth[ind,k] <- Da0[ind]*etk
        oo$Dmuth[ind,k] <- Dmua0[ind]*etk
        oo$Dmu2th[ind,k] <- Dmu2a0[ind]*etk 
      } 
    }  
    if (level >1) { ## and the second derivative components 
      oo$Dmu4 <- 2*((3*b^2 + 4*a*c)/f + a2*(6*a2/f - 12*b)/f2 - d)/f
      Dmu3a0 <-  2 * ((a0*c  + 3*c0*a + 3*b0*b)/f - d0  + 
                      6*a*(a0*a2/f - b0*a - a0*b)/f2 )/f
      Dmu3a1 <- 2 * (d1 - (a1*c + 3*(c1*a + b1*b))/f
                     + 6*a*(b1*a - a1*a2/f  + a1*b)/f2)/f

      Dmua0a0 <- 2*(c0 + (2*a0*(b0 - a0*a/f) - b0*a)/f )/f
      Dmua1a1 <- 2*( (b1*a + 2*a1*(b1 - a1*a/f))/f - c1)/f
      Dmua0a1 <- 2*(a0*(2*a1*a/f - b1) - b0*a1)/f2 

      Dmu2a0a0 <- 2*(d0 + (b0*(2*b0 - b)  + 2*c0*(a0 - a))/f +
                     2*(b0*a2 + a0*(3*a0*a2/f  - 4*b0*a - a0*b))/f2)/f

      Dmu2a1a1 <-  2*( (2*c1*(a + a1) + b1*(2*b1 + b))/f
                + 2*(a1*(3*a1*a2/f  - a1*b) - b1*a*(a + 4*a1))/f2 - d1)/f

      Dmu2a0a1 <- 0 ## (8*a0*b1*a/f^3 + 8*b0*a1*a/f^3 - 12*a0*a1*a/f^4 - 4*b0*b1/f^2 +
                    ## 4*a0*a1*b/f^3 - 2*c0*a1/f^2 - 2*c1*a0/f^2)

      Da0a0 <- 2 * (b0 + a0^2/f)/f #! + .5 * (C - B^2/A)/A
      Da1a1 <- -2* (b1 - a1^2/f)/f #! + .5 * (C - B^2/A)/A
      Da0a1 <- -2*a0*a1/f2 #! - .5 * (C - B^2/A)/A

      ## now transform to derivatives w.r.t. theta...
      n2d <- (R-2)*(R-1)/2
      oo$Dmu3th <- matrix(0,n,R-2)
      oo$Dmu2th2 <- oo$Dmuth2 <- oo$Dth2 <- matrix(0,n,n2d)
      i <- 0
      for (j in 1:(R-2)) for (k in j:(R-2)) { 
        i <- i + 1 ## the second deriv col
        ind <- y >= j ## rest are zero
        ar1.k <- ar.k <- rep(exp(theta[k]),n)
        ar.k[y==R | y <= k] <- 0; ar1.k[y<k+2] <- 0
        ar.j <- ar1.j <- rep(exp(theta[j]),n)
        ar.j[y==R | y <= j] <- 0; ar1.j[y<j+2] <- 0
        ar.kj <- ar1.kj <- rep(0,n)
        if (k==j) {
          ar.kj[y>k&y<R] <- exp(theta[k])
          ar1.kj[y>k+1] <- exp(theta[k])
          oo$Dmu3th[ind,k] <- Dmu3a1[ind]*ar.k[ind]  + Dmu3a0[ind]*ar1.k[ind]
        }
        oo$Dth2[,i] <- Da1a1*ar.k*ar.j + Da0a1*ar.k*ar1.j + Da1 * ar.kj +
                     Da0a0*ar1.k*ar1.j + Da0a1*ar1.k*ar.j + Da0 * ar1.kj
        oo$Dmuth2[,i] <- Dmua1a1*ar.k*ar.j + Dmua0a1*ar.k*ar1.j + Dmua1 * ar.kj +
                     Dmua0a0*ar1.k*ar1.j + Dmua0a1*ar1.k*ar.j + Dmua0 * ar1.kj
        oo$Dmu2th2[,i] <- Dmu2a1a1*ar.k*ar.j + Dmu2a0a1*ar.k*ar1.j + Dmu2a1 * ar.kj +
                     Dmu2a0a0*ar1.k*ar1.j + Dmu2a0a1*ar1.k*ar.j + Dmu2a0 * ar1.kj
      } 
    }
    oo
  } ## end of Dd
 
  aic <- function(y, mu, theta=NULL, wt, dev) {
  
   Fdiff <- function(a,b) {
      ## cancellation resistent F(b)-F(a), b>a 
      h <- rep(1,length(b)); h[b>0] <- -1; eb <- exp(b*h)
      h <- h*0+1; h[a>0] <- -1; ea <- exp(a*h)
      ind <- b<0;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- bi/(1+bi) - ai/(1+ai)
      ind1 <- a>0;bi <- eb[ind1];ai <- ea[ind1]
      h[ind1] <- (ai-bi)/((ai+1)*(bi+1))
      ind <- !ind & !ind1;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- (1-ai*bi)/((bi+1)*(ai+1))
      h
    }
    if (is.null(theta)) theta <- get(".Theta")
    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    al1 <- alpha[y+1];al0 = alpha[y]
    ##f1 <- F(al1-mu);f0 <- F(al0-mu);f <- f1 - f0
    f <- Fdiff(al0-mu,al1-mu)
    -2*sum(log(f))
  } ## end aic

  ls <- function(y,w,n,theta,scale) {
    ## the log saturated likelihood function. 
    #! actually only first line used since re-def as 0
    return(list(ls=0,lsth1=rep(0,R-2),lsth2=matrix(0,R-2,R-2)))
    F <- function(x) {
      h <- ind <- x > 0; h[ind] <- 1/(exp(-x[ind]) + 1)
      x <- exp(x[!ind]); h[!ind] <- (x/(1+x))
      h 
    } 
    Fdiff <- function(a,b) {
      ## cancellation resistent F(b)-F(a), b>a
      h <- rep(1,length(b)); h[b>0] <- -1; eb <- exp(b*h)
      h <- h*0+1; h[a>0] <- -1; ea <- exp(a*h)
      ind <- b<0;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- bi/(1+bi) - ai/(1+ai)
      ind1 <- a>0;bi <- eb[ind1];ai <- ea[ind1]
      h[ind1] <- (ai-bi)/((ai+1)*(bi+1))
      ind <- !ind & !ind1;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- (1-ai*bi)/((bi+1)*(ai+1))
      h
    }
    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    al1 <- alpha[y+1];al0 = alpha[y]
    g1 <- F((al1-al0)/2);g0 <- F((al0-al1)/2)
    ##A <- pmax(g1 - g0,.Machine$double.eps)
    A <- Fdiff((al0-al1)/2,(al1-al0)/2)
    ls <- sum(log(A))
    B <- g1^2 - g1 + g0^2 - g0 
    C <- 2 * g1^3 - 3 * g1^2 + g1 - 2 * g0^3 + 3 * g0^2 - g0
    Da0 <- .5 * B/A ; Da1 <- -0.5 *B/A 
    Da0a0 <- .25 * C/A - .25 * B^2/A^2
    Da1a1 <- .25 * C/A - .25 * B^2/A^2
    Da0a1 <- - .25 * C/A + .25 * B^2/A^2
    i <- 0 
    n2d <- (R-2)*(R-1)/2
    n <- length(y)
    Dth <- matrix(0,n,R-2)
    Dth2 <- matrix(0,n,n2d)
    for (j in 1:(R-2)) for (k in j:(R-2)) { 
      i <- i + 1 ## the second deriv col
      ind <- y >= j ## rest are zero
      ar1.k <- ar.k <- rep(exp(theta[k]),n)
      ar.k[y==R | y <= k] <- 0; ar1.k[y<k+2] <- 0
      ar.j <- ar1.j <- rep(exp(theta[j]),n)
      ar.j[y==R | y <= j] <- 0; ar1.j[y<j+2] <- 0
      ar.kj <- ar1.kj <- rep(0,n)
      if (k==j) {
        ar.kj[y>k&y<R] <- exp(theta[k])
        ar1.kj[y>k+1] <- exp(theta[k])
        Dth[ind,k] <- Da1[ind]*ar.k[ind]  + Da0[ind]*ar1.k[ind]
      }
      Dth2[,i] <- Da1a1*ar.k*ar.j + Da0a1*ar.k*ar1.j + Da1 * ar.kj +
                  Da0a0*ar1.k*ar1.j + Da0a1*ar1.k*ar.j + Da0 * ar1.kj
    } 
    lsth2=colSums(Dth2)
    if (R>2) {
      ls2 <- matrix(0,R-2,R-2);ii <- 0
      for (i in 1:(R-2)) for (j in i:(R-2)) { 
        ii <- ii + 1 
        ls2[i,j] <- ls2[j,i] <- lsth2[ii]
      } 
    }
    list(ls=ls,lsth1=colSums(Dth),lsth2=ls2)
   
  } ## end of ls
  
  ## initialization is interesting -- needs to be with reference to initial cut-points
  ## so that mu puts each obs in correct category initially...
 
  preinitialize <- expression({
    ocat.ini <- function(R,y) {
    ## initialize theta from raw counts in each category
      if (R<3) return
      y <- c(1:R,y) ## make sure there is *something* in each class
      p <- cumsum(tabulate(y[is.finite(y)])/length(y[is.finite(y)]))
      eta <- if (p[1]==0) 5 else -1 - log(p[1]/(1-p[1])) ## mean of latent
      theta <- rep(-1,R-1) 
      for (i in 2:(R-1)) theta[i] <- log(p[i]/(1-p[i])) + eta 
      theta <- diff(theta)
      theta[theta <= 0.01] <- 0.01
      theta <- log(theta)
    }
    R3 <- length(G$family$getTheta())+2
    if (R3>2&&G$family$n.theta>0) { 
      Theta <- ocat.ini(R3,G$y)
      G$family$putTheta(Theta)
    } 
  })

  initialize <- expression({ 
    R <- length(family$getTheta())+2 ## don't use n.theta as it's used to signal fixed theta
    if (any(y < 1)||any(y>R)) stop("values out of range")
    n <- rep(1, nobs)
    ## get the cut points... 
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -2;alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(family$getTheta()))
    } 
    alpha[R+1] <- alpha[R] + 1
    mustart <- (alpha[y+1] + alpha[y])/2
  })

  residuals <- function(object,type=c("deviance","working","response")) {
    if (type == "working") { 
      res <- object$residuals 
    } else if (type == "response") {
       theta <- object$family$getTheta()
       mu <- object$linear.predictors
       R = length(theta)+2
       alpha <- rep(0,R+1) ## the thresholds
       alpha[1] <- -Inf;alpha[R+1] <- Inf
       alpha[2] <- -1
       if (R > 2) { 
         ind <- 3:R
         alpha[ind] <- alpha[2] + cumsum(exp(theta))
       } 
       fv <- mu*NA
       for (i in 1:(R+1)) {
         ind <- mu>alpha[i] & mu<=alpha[i+1]
         fv[ind] <- i
       } 
       res <- object$y - fv
    } else if (type == "deviance") { 
      y <- object$y
      mu <- object$fitted.values
      wts <- object$prior.weights
      res <- object$family$dev.resids(y,mu,wts)
      s <- attr(res,"sign")
      if (is.null(s)) s <- sign(y-mu)
      res <- as.numeric(sqrt(pmax(res,0)) * s) 
      
    }
    res
  } ## residuals

 
  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                beta=NULL,off=NULL,Vb=NULL) {
  ## optional function to give predicted values - idea is that 
  ## predict.gam(...,type="response") will use this, and that
  ## either eta will be provided, or {X, beta, off, Vb}. family$data
  ## contains any family specific extra information. 
    ocat.prob <- function(theta,lp,se=NULL) {
    ## compute probabilities for each class in ocat model
    ## theta is finite cut points, lp is linear predictor, se 
    ## is standard error on lp...
      R <- length(theta) 
      dp <- prob <- matrix(0,length(lp),R+2)
      prob[,R+2] <- 1
      for (i in 1:R) {
        x <- theta[i] - lp
        ind <- x > 0
        prob[ind,i+1] <- 1/(1+exp(-x[ind]))
        ex <- exp(x[!ind])
        prob[!ind,i+1] <- ex/(1+ex)
        dp[,i+1] <- prob[,i+1]*(prob[,i+1]-1)
      }
      prob <- t(diff(t(prob)))
      dp <- t(diff(t(dp))) ## dprob/deta
      if (!is.null(se)) se <- as.numeric(se)*abs(dp)
      list(prob,se)
    } ## ocat.prob

    theta <- family$getTheta(TRUE)
    if (is.null(eta)) { ## return probabilities
      mu <- X%*%beta + off 
      se <- if (se) sqrt(pmax(0,rowSums((X%*%Vb)*X))) else NULL
      ##theta <- cumsum(c(-1,exp(theta)))
      p <- ocat.prob(theta,mu,se)
      if (is.null(se)) return(p) else { ## approx se on prob also returned
        names(p) <- c("fit","se.fit")
        return(p)
      } 
    } else { ## return category implied by eta (i.e mean of latent)
      R = length(theta)+2
      alpha <- rep(0,R) ## the thresholds
      alpha[1] <- -Inf;alpha[R] <- Inf
      fv <- eta*NA
      for (i in 1:(R+1)) {
        ind <- eta>alpha[i] & eta<=alpha[i+1]
        fv[ind] <- i
      } 
      return(fv)
    }
  } ## predict

  rd <- function(mu,wt,scale) {
    ## simulate data given fitted latent variable in mu 
    theta <- get(".Theta") 
    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    ## ... cut points computed, now simulate latent variable, u
    y <- u <- runif(length(mu))
    u <- mu + log(u/(1-u)) 
    ## and allocate categories according to u and cut points...
    for (i in 1:R) {
      y[u > alpha[i]&u <= alpha[i+1]] <- i
    }
    y
  }

  environment(dev.resids) <- environment(aic) <- environment(putTheta) <-
  environment(getTheta) <- environment(rd) <- environment(predict) <- env
  structure(list(family = "Ordered Categorical", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,postproc=postproc,
        preinitialize = preinitialize, ls=ls,rd=rd,residuals=residuals,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta,
        ini.theta = iniTheta,putTheta=putTheta,predict=predict,step = 1,
        getTheta=getTheta,no.r.sq=TRUE), class = c("extended.family","family"))
} ## end of ocat

#######################
## negative binomial...
#######################

nb <- function (theta = NULL, link = "log") { 
## Extended family object for negative binomial, to allow direct estimation of theta
## as part of REML optimization. Currently the template for extended family objects.
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("log", "identity", "sqrt")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"")
  }
  ## Theta <-  NULL;
  n.theta <- 1
  if (!is.null(theta)&&theta!=0) {
      if (theta>0) { 
        iniTheta <- log(theta) ## fixed theta supplied
        n.theta <- 0 ## signal that there are no theta parameters to estimate
      } else iniTheta <- log(-theta) ## initial theta supplied
  } else iniTheta <- 0 ## inital log theta value
    
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", iniTheta, envir = env)
    getTheta <- function(trans=FALSE) if (trans) exp(get(".Theta")) else get(".Theta") # get(".Theta")
    putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

    variance <- function(mu) mu + mu^2/exp(get(".Theta"))

    validmu <- function(mu) all(mu > 0)

    dev.resids <- function(y, mu, wt,theta=NULL) {
      if (is.null(theta)) theta <- get(".Theta")
      theta <- exp(theta) ## note log theta supplied
      2 * wt * (y * log(pmax(1, y)/mu) - 
        (y + theta) * log((y + theta)/(mu + theta))) 
    }
    
    Dd <- function(y, mu, theta, wt, level=0) {
    ## derivatives of the deviance...
      ##ltheta <- theta
      theta <- exp(theta)
      yth <- y + theta
      muth <- mu + theta
      r <- list()
      ## get the quantities needed for IRLS. 
      ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
      ## Dmu is deriv w.r.t. mu once, etc...
      r$Dmu <- 2 * wt * (yth/muth - y/mu)
      r$Dmu2 <- -2 * wt * (yth/muth^2 - y/mu^2)
      r$EDmu2 <- 2 * wt * (1/mu - 1/muth) ## exact (or estimated) expected weight
      if (level>0) { ## quantities needed for first derivatives
        r$Dth <- -2 * wt * theta * (log(yth/muth) + (1 - yth/muth) ) 
        r$Dmuth <- 2 * wt * theta * (1 - yth/muth)/muth
        r$Dmu3 <- 4 * wt * (yth/muth^3 - y/mu^3)
        r$Dmu2th <- 2 * wt * theta * (2*yth/muth - 1)/muth^2
      } 
      if (level>1) { ## whole damn lot
        r$Dmu4 <- 2 * wt * (6*y/mu^4 - 6*yth/muth^4)
        r$Dth2 <- -2 * wt * theta * (log(yth/muth) +
                     theta*yth/muth^2 - yth/muth - 2*theta/muth + 1 +
                     theta /yth)
        r$Dmuth2 <- 2 * wt * theta * (2*theta*yth/muth^2 - yth/muth - 2*theta/muth + 1)/muth
        r$Dmu2th2 <- 2 * wt * theta * (- 6*yth*theta/muth^2 + 2*yth/muth + 4*theta/muth - 1) /muth^2
        r$Dmu3th <- 4 * wt * theta * (1 - 3*yth/muth)/muth^3
      }
      r
    }

    aic <- function(y, mu, theta=NULL, wt, dev) {
        if (is.null(theta)) theta <- get(".Theta")
        Theta <- exp(theta)
        term <- (y + Theta) * log(mu + Theta) - y * log(mu) +
            lgamma(y + 1) - Theta * log(Theta) + lgamma(Theta) -
            lgamma(Theta + y)
        2 * sum(term * wt)
    }
    
    ls <- function(y,w,n,theta,scale) {
       ## the log saturated likelihood function.
       Theta <- exp(theta)
       ylogy <- y;ind <- y>0;ylogy[ind] <- y[ind]*log(y[ind])
       term <- (y + Theta) * log(y + Theta) - ylogy +
            lgamma(y + 1) - Theta * log(Theta) + lgamma(Theta) -
            lgamma(Theta + y)
       ls <- -sum(term*w)
       ## first derivative wrt theta...
       yth <- y+Theta
       lyth <- log(yth)
       psi0.yth <- digamma(yth) 
       psi0.th <- digamma(Theta)
       term <- Theta * (lyth - psi0.yth + psi0.th-theta)
       lsth <- -sum(term*w)
       ## second deriv wrt theta...
       psi1.yth <- trigamma(yth) 
       psi1.th <- trigamma(Theta)
       term <- Theta * (lyth - Theta*psi1.yth - psi0.yth + Theta/yth + Theta * psi1.th + psi0.th - theta -1)        
       lsth2 <- -sum(term*w)
       list(ls=ls, ## saturated log likelihood
            lsth1=lsth, ## first deriv vector w.r.t theta - last element relates to scale, if free
            lsth2=lsth2) ## Hessian w.r.t. theta, last row/col relates to scale, if free
    }

    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
        n <- rep(1, nobs)
        mustart <- y + (y == 0)/6
    })
  
    postproc <- expression({
      object$family$family <- 
      paste("Negative Binomial(",round(object$family$getTheta(TRUE),3),")",sep="")
    })

    rd <- function(mu,wt,scale) {
      Theta <- exp(get(".Theta"))
      rnbinom(mu,size=Theta,mu=mu)
    }

    qf <- function(p,mu,wt,scale) {
      Theta <- exp(get(".Theta"))
      qnbinom(p,size=Theta,mu=mu)
    }
 

     environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
     environment(rd)<- environment(qf)<- environment(variance) <- environment(putTheta) <- env
    structure(list(family = "negative binomial", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,variance=variance,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,postproc=postproc,ls=ls,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta, 
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,rd=rd,qf=qf),
        class = c("extended.family","family"))
} ## nb


## Tweedie....


tw <- function (theta = NULL, link = "log",a=1.01,b=1.99) { 
## Extended family object for Tweedie, to allow direct estimation of p
## as part of REML optimization. 
## p = (a+b*exp(theta))/(1+exp(theta)), i.e. a < p < b
## NOTE: The Tweedie density computation at low phi, low p is susceptible
##       to cancellation error, which seems unavoidable. Furthermore 
##       there are known problems with spurious maxima in the likelihood 
##       w.r.t. p when the data are rounded. 
  
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("log", "identity", "sqrt","inverse")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else  stop(gettextf("link \"%s\" not available for Tweedie family.", 
                linktemp, collapse = ""), domain = NA)
  }
  ## Theta <-  NULL;
  n.theta <- 1
  if (!is.null(theta)&&theta!=0) {
      if (abs(theta)<=a||abs(theta)>=b) stop("Tweedie p must be in interval (a,b)")
      if (theta>0) { ## fixed theta supplied
        iniTheta <- log((theta-a)/(b-theta)) 
        n.theta <- 0 ## so no theta to estimate
      } else iniTheta <- log((-theta-a)/(b+theta)) ## initial theta supplied
  } else iniTheta <- 0 ## inital log theta value
    
  env <- new.env(parent = .GlobalEnv)
  assign(".Theta", iniTheta, envir = env) 
  assign(".a",a, envir = env);assign(".b",b, envir = env)
  getTheta <- function(trans=FALSE) { 
  ## trans transforms to the original scale...
    th <- get(".Theta")
    a <- get(".a");b <- get(".b")
    if (trans) th <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1)
    th
  }
  putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

  validmu <- function(mu) all(mu > 0)
  
  variance <- function(mu) { 
    th <- get(".Theta");a <- get(".a");b <- get(".b")
    p <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1)
    mu^p
  }
  
  dev.resids <- function(y, mu, wt,theta=NULL) {
    if (is.null(theta)) theta <- get(".Theta")
    a <- get(".a");b <- get(".b")
    p <- if (theta>0) (b+a*exp(-theta))/(1+exp(-theta)) else (b*exp(theta)+a)/(exp(theta)+1)
    y1 <- y + (y == 0)
    theta <- if (p == 1) log(y1/mu) else (y1^(1 - p) - mu^(1 - p))/(1 - p)
    kappa <- if (p == 2) log(y1/mu) else (y^(2 - p) - mu^(2 - p))/(2 - p)
    2 * (y * theta - kappa) * wt
  }
    
  Dd <- function(y, mu, theta, wt, level=0) {
  ## derivatives of the deviance...
    a <- get(".a");b <- get(".b")
    th <- theta
    p <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1)
    dpth1 <- if (th>0) exp(-th)*(b-a)/(1+exp(-th))^2 else exp(th)*(b-a)/(exp(th)+1)^2
    dpth2 <- if (th>0) ((a-b)*exp(-th)+(b-a)*exp(-2*th))/(exp(-th)+1)^3 else
                   ((a-b)*exp(2*th)+(b-a)*exp(th))/(exp(th)+1)^3
    mu1p <- mu^(1-p)
    mup <- mu^p
    r <- list()
    ## get the quantities needed for IRLS. 
    ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
    ## Dmu is deriv w.r.t. mu once, etc...
    ymupi <- y/mup
    r$Dmu <- 2*wt*(mu1p - ymupi)
    r$Dmu2 <- 2*wt*(mu^(-1-p)*p*y + (1-p)/mup)
    r$EDmu2 <- (2*wt)/mup ## expected Dmu2 (weight)
   
    if (level>0) { ## quantities needed for first derivatives
        i1p <- 1/(1-p)
        y1 <- y + (y==0)
        ##ylogy <- y*log(y1)
        logmu <- log(mu)
        mu2p <- mu * mu1p
        r$Dth <- 2 * wt * ( (y^(2-p)*log(y1) - mu2p*logmu)/(2-p) + 
                            (y*mu1p*logmu - y^(2-p)*log(y1))/(1-p) -
                            (y^(2-p) - mu2p)/(2-p)^2 + 
                            (y^(2-p) - y*mu1p)*i1p^2) *dpth1       

        r$Dmuth <- 2 * wt * logmu * (ymupi - mu1p)*dpth1
        mup1 <-  mu^(-p-1)
        r$Dmu3 <- -2 * wt * mup1*p*(y/mu*(p+1) + 1-p)    
        r$Dmu2th <- 2 * wt  * (mup1*y*(1-p*logmu)-(logmu*(1-p)+1)/mup )*dpth1
      } 
      if (level>1) { ## whole damn lot
        mup2 <- mup1/mu
        r$Dmu4 <- 2 * wt * mup2*p*(p+1)*(y*(p+2)/mu + 1 - p)
        y2plogy <- y^(2-p)*log(y1);y2plog2y <- y2plogy*log(y1)
       
        r$Dth2 <- 2 * wt * (((mu2p*logmu^2-y2plog2y)/(2-p) + (y2plog2y - y*mu1p*logmu^2)/(1-p) +
                             2*(y2plogy-mu2p*logmu)/(2-p)^2 + 2*(y*mu1p*logmu-y2plogy)/(1-p)^2
                             + 2 * (mu2p - y^(2-p))/(2-p)^3+2*(y^(2-p)-y*mu^(1-p))/(1-p)^3)*dpth1^2) +
                              r$Dth*dpth2/dpth1
        
        r$Dmuth2 <- 2 * wt * ((mu1p * logmu^2 - logmu^2*ymupi)*dpth1^2) + r$Dmuth*dpth2/dpth1

        r$Dmu2th2 <- 2 * wt * ( (mup1 * logmu*y*(logmu*p - 2) + logmu/mup*(logmu*(1-p) + 2))*dpth1^2) +
                                r$Dmu2th * dpth2/dpth1

        r$Dmu3th <- 2 * wt * mup1*(y/mu*(logmu*(1+p)*p-p -p-1) +logmu*(1-p)*p + p - 1 + p)*dpth1
      }
      r
    }

    aic <- function(y, mu, theta=NULL, wt, dev) {
        if (is.null(theta)) theta <- get(".Theta")  
        a <- get(".a");b <- get(".b")
        p <- if (theta>0) (b+a*exp(-theta))/(1+exp(-theta)) else (b*exp(theta)+a)/(exp(theta)+1)
        scale <- dev/sum(wt)
        -2 * sum(ldTweedie(y, mu, p = p, phi = scale)[, 1] * 
            wt) + 2
    }

    ls <- function(y, w, n, theta, scale) {
        ## evaluate saturated log likelihood + derivs w.r.t. working params and log(scale)
        a <- get(".a");b <- get(".b")
        LS <- colSums(w * ldTweedie(y, y, rho=log(scale), theta=theta,a=a,b=b))
        lsth1 <- c(LS[4],LS[2])
        lsth2 <- matrix(c(LS[5],LS[6],LS[6],LS[3]),2,2)
        list(ls=LS[1],lsth1=lsth1,lsth2=lsth2)
    }

 
    initialize <- expression({
        n <- rep(1, nobs)
        mustart <- y + (y == 0)*.1
    })
    postproc <- expression({
      object$family$family <- 
      paste("Tweedie(p=",round(object$family$getTheta(TRUE),3),")",sep="")
    })
  
    rd <- function(mu,wt,scale) {
     th <- get(".Theta") 
     a <- get(".a");b <- get(".b")
     p <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1)

     if (p == 2) 
            rgamma(mu, shape = 1/scale, scale = mu * scale)
     else
            rTweedie(mu, p = p, phi = scale)
    }

     environment(Dd) <- environment(ls) <-
     environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
      environment(rd) <- environment(variance) <- environment(putTheta) <- env
    structure(list(family = "Tweedie", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,variance=variance,rd=rd,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,postproc=postproc,ls=ls,
        validmu = validmu, valideta = stats$valideta,canonical="none",n.theta=n.theta, 
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,scale = -1),
        class = c("extended.family","family"))
} ## tw

## beta regression

betar <- function (theta = NULL, link = "logit",eps=.Machine$double.eps*100) { 
## Extended family object for beta regression
## length(theta)=1; log theta supplied
## This serves as a prototype for working with -2logLik
## as deviance, and only dealing with saturated likelihood 
## at the end.
## Written by Natalya Pya. 'saturated.ll' by Simon Wood 
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("logit", "probit", "cloglog", "cauchit")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for beta regression; available links are  \"logit\", \"probit\", \"cloglog\" and \"cauchit\"")
    }
   
    n.theta <- 1
    if (!is.null(theta)&&theta!=0) {
       if (theta>0) {
           iniTheta <- log(theta) ## fixed theta supplied
           n.theta <- 0 ## signal that there are no theta parameters to estimate
       } else iniTheta <- log(-theta) ## initial theta supplied
    } else iniTheta <- 0 ##  inital log theta value
    
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", iniTheta, envir = env)
    assign(".betarEps",eps, envir = env)
    getTheta <- function(trans=FALSE) if (trans) exp(get(".Theta")) else get(".Theta") 
    putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

    variance <- function(mu) { 
        th <- get(".Theta")
        mu*(1 - mu)/(1+exp(th))
    }
  
    validmu <- function(mu) all(mu > 0 & mu < 1)

    dev.resids <- function(y, mu, wt,theta=NULL) {
    ## '-2*loglik' instead of deviance in REML/ML expression
      if (is.null(theta)) theta <- get(".Theta")
      theta <- exp(theta) ## note log theta supplied
      muth <- mu*theta
      ## yth <- y*theta
      2* wt * (-lgamma(theta) +lgamma(muth) + lgamma(theta - muth) - muth*(log(y)-log1p(-y)) - 
                theta*log1p(-y) + log(y) + log1p(-y)) 
    }
    
    Dd <- function(y, mu, theta, wt, level=0) {
    ## derivatives of the -2*loglik...
      ## ltheta <- theta
      theta <- exp(theta)
      onemu <- 1 - mu; ## oney <- 1 - y
      muth <- mu*theta; ## yth <- y*theta
      onemuth <- onemu*theta  ## (1-mu)*theta
      psi0.th <- digamma(theta)
      psi1.th <- trigamma(theta)
      psi0.muth <- digamma(muth) 
      psi0.onemuth <- digamma(onemuth)
      psi1.muth <- trigamma(muth)
      psi1.onemuth <- trigamma(onemuth)
      psi2.muth <- psigamma(muth,2)
      psi2.onemuth <- psigamma(onemuth,2)
      psi3.muth <- psigamma(muth,3)
      psi3.onemuth <- psigamma(onemuth,3)
      log.yoney <- log(y)-log1p(-y)
      r <- list()
      ## get the quantities needed for IRLS. 
      ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
      ## Dmu is deriv w.r.t. mu once, etc...
      r$Dmu <- 2 * wt * theta* (psi0.muth - psi0.onemuth - log.yoney)
      r$Dmu2 <- 2 * wt * theta^2*(psi1.muth+psi1.onemuth)
      r$EDmu2 <- r$Dmu2
      if (level>0) { ## quantities needed for first derivatives
        r$Dth <- 2 * wt *theta*(-mu*log.yoney - log1p(-y)+ mu*psi0.muth+onemu*psi0.onemuth -psi0.th) 
        r$Dmuth <- r$Dmu + 2 * wt * theta^2*(mu*psi1.muth -onemu*psi1.onemuth)
        r$Dmu3 <- 2 * wt *theta^3 * (psi2.muth - psi2.onemuth) 
        r$Dmu2th <- 2* r$Dmu2 + 2 * wt * theta^3* (mu*psi2.muth + onemu*psi2.onemuth)
      } 
      if (level>1) { ## whole lot
        r$Dmu4 <- 2 * wt *theta^4 * (psi3.muth+psi3.onemuth) 
        r$Dth2 <- r$Dth +2 * wt *theta^2* (mu^2*psi1.muth+ onemu^2*psi1.onemuth-psi1.th)
        r$Dmuth2 <- r$Dmuth + 2 * wt *theta^2* (mu^2*theta*psi2.muth+ 2*mu*psi1.muth -
                    theta*onemu^2*psi2.onemuth - 2*onemu*psi1.onemuth)
        r$Dmu2th2 <- 2*r$Dmu2th + 2* wt * theta^3* (mu^2*theta*psi3.muth +3*mu*psi2.muth+ 
                    onemu^2*theta*psi3.onemuth + 3*onemu*psi2.onemuth )
        r$Dmu3th <- 3*r$Dmu3 + 2 * wt *theta^4*(mu*psi3.muth-onemu*psi3.onemuth)
      }
      r
    }

    aic <- function(y, mu, theta=NULL, wt, dev) {
        if (is.null(theta)) theta <- get(".Theta")
        theta <- exp(theta)
        muth <- mu*theta
        term <- -lgamma(theta)+lgamma(muth)+lgamma(theta-muth)-(muth-1)*log(y)-
               (theta-muth-1)*log1p(-y) ## `-' log likelihood for each observation
        2 * sum(term * wt)
    }
    
    ls <- function(y,w,n,theta,scale) {
       ## the log saturated likelihood function.
       ## ls is defined as zero for REML/ML expression as deviance is defined as -2*log.lik 
       list(ls=0,## saturated log likelihood
            lsth1=0,  ## first deriv vector w.r.t theta - last element relates to scale
            lsth2=0) ##Hessian w.r.t. theta
     }

   
    ## preinitialization to reset G$y values of <=0 and >=1... 
    ## code to evaluate in estimate.gam...
    ## reset G$y values of <=0 and >= 1 to eps and 1-eps... 
    preinitialize <- NULL ## keep codetools happy
    eval(parse(text=paste("preinitialize <- expression({\n eps <- ",eps,
         "\n G$y[G$y >= 1-eps] <- 1 - eps\n  G$y[G$y<= eps] <- eps })")))
    
 #   preinitialize <- expression({
 #     eps <- 1e-7 
 #     G$y[G$y >= 1-eps] <- 1 - eps
 #     G$y[G$y<= eps] <- eps
 #   })

    saturated.ll <- function(y,wt,theta=NULL){
    ## function to find the saturated loglik by Newton method,
    ## searching for the mu (on logit scale) that max loglik given theta and data...

      gbh <- function(y,eta,phi,deriv=FALSE,a=1e-8,b=1-a) {
      ## local function evaluating log likelihood (l), gradient and second deriv 
      ## vectors for beta... a and b are min and max mu values allowed. 
      ## mu = (a + b*exp(eta))/(1+exp(eta))
        ind <- eta>0
        expeta <- mu <- eta;
        expeta[ind] <- exp(-eta[ind]);expeta[!ind] <- exp(eta[!ind])
        mu[ind] <- (a*expeta[ind] + b)/(1+expeta[ind])
        mu[!ind] <- (a + b*expeta[!ind])/(1+expeta[!ind])
        l <- dbeta(y,phi*mu,phi*(1-mu),log=TRUE)
        if (deriv) {
          g <- phi * log(y) - phi * log1p(-y) - phi * digamma(mu*phi) + phi * digamma((1-mu)*phi)
          h <- -phi^2*(trigamma(mu*phi)+trigamma((1-mu)*phi))
          dmueta2 <- dmueta1 <-  eta
          dmueta1 <- expeta*(b-a)/(1+expeta)^2 
          dmueta2 <- sign(eta)* ((a-b)*expeta+(b-a)*expeta^2)/(expeta+1)^3
          h <- h * dmueta1^2 + g * dmueta2
          g <- g * dmueta1
        } else g=h=NULL
        list(l=l,g=g,h=h,mu=mu)
      } ## gbh 
      ## now Newton loop...
      eps <- get(".betarEps")
      eta <- y
      a <- eps;b <- 1 - eps
      y[y<eps] <- eps;y[y>1-eps] <- 1-eps
      eta[y<=eps*1.2] <- eps *1.2
      eta[y>=1-eps*1.2] <- 1-eps*1.2
      eta <- log((eta-a)/(b-eta)) 
      mu <- LS <- ii <- 1:length(y)
      for (i in 1:200) {
        ls <- gbh(y,eta,theta,TRUE,a=eps/10)
        conv <- abs(ls$g)<mean(abs(ls$l)+.1)*1e-8
        if (sum(conv)>0) { ## some convergences occured
          LS[ii[conv]] <- ls$l[conv] ## store converged
          mu[ii[conv]] <- ls$mu[conv] ## store mu at converged
          ii <- ii[!conv] ## drop indices
          if (length(ii)>0) { ## drop the converged
            y <- y[!conv];eta <- eta[!conv]
            ls$l <- ls$l[!conv];ls$g <- ls$g[!conv];ls$h <- ls$h[!conv]
          } else break ## nothing left to do
        }
        h <- -ls$h
        hmin <- max(h)*1e-4 
        h[h<hmin] <- hmin ## make +ve def
        delta <- ls$g/h   ## step
        ind <- abs(delta)>2
        delta[ind] <- sign(delta[ind])*2 ## step length limit
        ls1 <- gbh(y,eta+delta,theta,FALSE,a=eps/10); ## did it work?
        ind <- ls1$l<ls$l ## failure index
        k <- 0
        while (sum(ind)>0&&k<20) { ## step halve only failed steps
          k <- k + 1
          delta[ind] <- delta[ind]/2
          ls1$l[ind] <- gbh(y[ind],eta[ind]+delta[ind],theta,FALSE,a=eps/10)$l
          ind <- ls1$l<ls$l
        }
        eta <- eta + delta
      } ## end newton loop
      if (length(ii)>0) { 
        LS[ii] <- ls$l
        warning("saturated likelihood may be inaccurate")
      }

      list(f=sum(wt*LS),term=LS,mu=mu) ## fields f (sat lik) and term (individual datum sat lik) expected
    } ## saturated.ll


    postproc <- expression({
    ## code to evaluate in estimate.gam, to find the saturated
    ## loglik by Newton method
    ## searching for the mu (on logit scale) that max loglik given theta...
      wts <- object$prior.weights
      theta <- object$family$getTheta(trans=TRUE) ## exp theta
      lf <- object$family$saturated.ll(G$y, wts,theta)
      ## storing the saturated loglik for each datum...
      object$family$data <- list(ls = lf$term,mu.ls = lf$mu)   
      l2 <- object$family$dev.resids(G$y,object$fitted.values,wts)
      object$deviance <- 2*lf$f + sum(l2)
      wtdmu <- if (G$intercept) sum(wts * G$y)/sum(wts) 
              else object$family$linkinv(G$offset)
      object$null.deviance <- 2*lf$f + sum(object$family$dev.resids(G$y, wtdmu, wts))
      object$family$family <- 
      paste("Beta regression(",round(theta,3),")",sep="")
    })

    initialize <- expression({
        n <- rep(1, nobs)
        mustart <- y 
    })

    residuals <- function(object,type=c("deviance","working","response","pearson")) {
      if (type == "working") { 
        res <- object$residuals 
      } else if (type == "response") {
        res <- object$y - object$fitted.values
      } else if (type == "deviance") { 
        y <- object$y
        mu <- object$fitted.values
        wts <- object$prior.weights
     #   sim <- attr(y,"simula")
     #   if (!is.null(sim)) {  ## if response values simulated, Newton search called to get saturated log.lik
           lf <- object$family$saturated.ll(y, wts,object$family$getTheta(TRUE))
           object$family$data$ls <- lf$term  
     #   }
        res <- 2*object$family$data$ls + object$family$dev.resids(y,mu,wts)
        res[res<0] <- 0
        s <- sign(y-mu)
        res <- sqrt(res) * s   
      } else if (type == "pearson") {
        mu <- object$fitted.values
        res <- (object$y - mu)/object$family$variance(mu)^.5
      }
      res
     } ## residuals

    rd <- function(mu,wt,scale) {
     ## simulate data given fitted latent variable in mu 
      Theta <- exp(get(".Theta"))
      r <- rbeta(mu,shape1=Theta*mu,shape2=Theta*(1-mu))
      eps <- get(".betarEps")
      r[r>=1-eps] <- 1 - eps
      r[r<eps] <- eps
      r
    }

    qf <- function(p,mu,wt,scale) {
      Theta <- exp(get(".Theta"))
      q <- qbeta(p,shape1=Theta*mu,shape2=Theta*(1-mu))
      eps <-  get(".betarEps")
      q[q>=1-eps] <- 1 - eps
      q[q<eps] <- eps
      q
    }

    environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
    environment(rd)<- environment(qf) <- environment(variance) <- environment(putTheta) <-
    environment(saturated.ll) <- env

    structure(list(family = "Beta regression", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd, variance=variance,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,ls=ls,
        preinitialize=preinitialize,postproc=postproc, residuals=residuals, saturated.ll=saturated.ll,
        validmu = validmu, valideta = stats$valideta, n.theta=n.theta,  
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,rd=rd,qf=qf), 
        class = c("extended.family","family"))
    
} ## betar


  
## scaled t (Natalya Pya) ...

scat <- function (theta = NULL, link = "identity") { 
## Extended family object for scaled t distribution
## length(theta)=2; log theta supplied. 
## Written by Natalya Pya.
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("identity", "log", "inverse")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for scaled t distribution; available links are \"identity\", \"log\",  and \"inverse\"")
    }
    ## Theta <-  NULL;
    n.theta <- 2
    if (!is.null(theta)&&sum(theta==0)==0) {
      if (abs(theta[1]<2)) stop("scaled t df must be >2")
      if (sum(theta<0)) { 
        iniTheta <- c(log(abs(theta[1])-2),log(abs(theta[2]))) ## initial theta supplied
      } else { ## fixed theta supplied
        iniTheta <- c(log(theta[1]-2),log(theta[2])) 
        n.theta <- 0 ## no thetas to estimate
      }
    } else iniTheta <- c(-2,-1) ## inital log theta value
               
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", iniTheta, envir = env)
    getTheta <- function(trans=FALSE) { 
    ## trans transforms to the original scale...
      th <- get(".Theta")
      if (trans) { th <- exp(th); th[1] <- th[1] + 2  }
      th
    }
    putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

    variance <- function(mu) { 
        th <- get(".Theta")
        nu <- exp(th[1])+2; sig <- exp(th[2])
        sig^2*nu/(nu-2)
    }

   validmu <- function(mu) all(is.finite(mu))

   dev.resids <- function(y, mu, wt,theta=NULL) {
      if (is.null(theta)) theta <- get(".Theta")
      nu <- exp(theta[1])+2; sig <- exp(theta[2])
      wt * (nu + 1)*log1p((1/nu)*((y-mu)/sig)^2)
    }
    
    Dd <- function(y, mu, theta, wt, level=0) {
    ## derivatives of the deviance...
      ## ltheta <- theta
      nu <- exp(theta[1])+2; sig <- exp(theta[2])
      nu1 <- nu + 1;  ym <- y - mu; nu2 <- nu - 2;
      a <- 1 + (ym/sig)^2/nu
      oo <- list()
      ## get the quantities needed for IRLS. 
      ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
      ## Dmu is deriv w.r.t. mu once, etc...
      nu1ym <- nu1*ym
      sig2a <- sig^2*a
      nusig2a <- nu*sig2a
      f <- nu1ym/nusig2a
      f1 <- ym/nusig2a
      oo$Dmu <- -2 * wt * f
      oo$Dmu2 <- 2 * wt * nu1*(1/nusig2a- 2*f1^2)  # - 2*ym^2/(nu^2*sig^4*a^2)
      term <- 2*nu1/sig^2/(nu+3)
      n <- length(y) 
      oo$EDmu2 <- rep(term,n)
      if (level>0) { ## quantities needed for first derivatives
        nu1nusig2a <- nu1/nusig2a
        nu2nu <- nu2/nu
        fym <- f*ym; ff1 <- f*f1; f1ym <- f1*ym; fymf1 <- fym*f1
        ymsig2a <- ym/sig2a
        oo$Dmu2th <- oo$Dmuth <- oo$Dth <- matrix(0,n,2)
        oo$Dth[,1] <- 1 * wt * nu2 * (log(a) - fym/nu) 
        oo$Dth[,2] <- -2 * wt * fym    
        oo$Dmuth[,1] <- 2 * wt *(f - ymsig2a - fymf1)*nu2nu
        oo$Dmuth[,2] <- 4* wt* f* (1- f1ym)
        oo$Dmu3 <- 4 * wt * f * (3/nusig2a - 4*f1^2) 
        oo$Dmu2th[,1] <- 2* wt * (-nu1nusig2a + 1/sig2a + 5*ff1- 2*f1ym/sig2a - 4*fymf1*f1)*nu2nu
        oo$Dmu2th[,2] <- 4*wt*(-nu1nusig2a + ff1*5 - 4*ff1*f1ym)
      } 
      if (level>1) { ## whole lot
        ## nu1nu2 <- nu1*nu2; 
        nu1nu <- nu1/nu
        fymf1ym <- fym*f1ym; f1ymf1 <- f1ym*f1
        oo$Dmu4 <- 12 * wt * (-nu1nusig2a/nusig2a + 8*ff1/nusig2a - 8*ff1 *f1^2) 
        n2d <- 3 # number of the 2nd order derivatives
        oo$Dmu3th <- matrix(0,n,2)
        oo$Dmu2th2 <- oo$Dmuth2 <- oo$Dth2 <- matrix(0,n,n2d)
        oo$Dmu3th[,1] <- 4*wt*(-6*f/nusig2a + 3*f1/sig2a + 18*ff1*f1 - 4*f1ymf1/sig2a - 12*nu1ym*f1^4)*nu2nu
        oo$Dmu3th[,2] <- 48*wt* f* (- 1/nusig2a + 3*f1^2 -  2*f1ymf1*f1)
        
        oo$Dth2[,1] <- 1*wt *(nu2*log(a) +nu2nu*ym^2*(-2*nu2-nu1+ 2*nu1*nu2nu - nu1*nu2nu*f1ym)/nusig2a) ## deriv of D w.r.t. theta1 theta1 
  
        oo$Dth2[,2] <- 2*wt*(fym - ym*ymsig2a - fymf1ym)*nu2nu  ## deriv of D wrt theta1 theta2
        oo$Dth2[,3] <- 4 * wt * fym *(1 - f1ym)  ## deriv of D wrt theta2 theta2

        term <- 2*nu2nu - 2*nu1nu*nu2nu -1 + nu1nu
        oo$Dmuth2[,1] <- 2*wt*f1*nu2*(term - 2*nu2nu*f1ym + 4*fym*nu2nu/nu - fym/nu - 2*fymf1ym*nu2nu/nu)
 
        oo$Dmuth2[,2] <- 4*wt* (-f + ymsig2a + 3*fymf1 - ymsig2a*f1ym - 2*fymf1*f1ym)*nu2nu

        oo$Dmuth2[,3] <- 8*wt* f * (-1 + 3*f1ym - 2*f1ym^2)

        oo$Dmu2th2[,1] <- 2*wt*nu2*(-term + 10*nu2nu*f1ym -
               16*fym*nu2nu/nu - 2*f1ym + 5*nu1nu*f1ym - 8*nu2nu*f1ym^2 + 26*fymf1ym*nu2nu/nu - 
               4*nu1nu*f1ym^2 - 12*nu1nu*nu2nu*f1ym^3)/nusig2a
       
        oo$Dmu2th2[,2] <- 4*wt*(nu1nusig2a - 1/sig2a - 11*nu1*f1^2 + 5*f1ym/sig2a + 22*nu1*f1ymf1*f1 - 
                    4*f1ym^2/sig2a - 12*nu1*f1ymf1^2)*nu2nu
        oo$Dmu2th2[,3] <- 8*wt * (nu1nusig2a - 11*nu1*f1^2 + 22*nu1*f1ymf1*f1 - 12*nu1*f1ymf1^2)

      }
      oo
    }  ## end of Dd

 
    aic <- function(y, mu, theta=NULL, wt, dev) {
        if (is.null(theta)) theta <- get(".Theta")
        nu <- exp(theta[1])+2; sig <- exp(theta[2])
        term <- -lgamma((nu+1)/2)+ lgamma(nu/2) + log(sig*(pi*nu)^.5) +
           (nu+1)*log1p(((y-mu)/sig)^2/nu)/2  ## `-'log likelihood for each observation
        2 * sum(term * wt)
    }
    
    ls <- function(y,w,n,theta,scale) {
       ## the log saturated likelihood function.
       nu <- exp(theta[1])+2; sig <- exp(theta[2]); nu2 <- nu-2;
       nu2nu <- nu2/nu; nu12 <- (nu+1)/2
       term <- lgamma(nu12) - lgamma(nu/2) - log(sig*(pi*nu)^.5)
       ls <- sum(term*w) 
       ## first derivative wrt theta...
       lsth <- rep(0,2) 
       lsth2 <- matrix(0,2,2)  ## rep(0, 3)
       term <- nu2 * digamma(nu12)/2- nu2 * digamma(nu/2)/2 - 0.5*nu2nu
       lsth[1] <- sum(w*term)
       lsth[2] <- sum(-1*w)
       
       ## second deriv...      
       term <-  nu2^2 * trigamma(nu12)/4 + nu2 * digamma(nu12)/2 -
           nu2^2 * trigamma(nu/2)/4 - nu2 * digamma(nu/2)/2 + 0.5*(nu2nu)^2 - 0.5*nu2nu
       lsth2[1,1] <- sum(term*w)
       lsth2[1,2] <- lsth2[2,1] <- lsth2[2,2] <- 0
       list(ls=ls,## saturated log likelihood
            lsth1=lsth, ## first derivative vector wrt theta
            lsth2=lsth2) ## Hessian wrt theta
    }

    preinitialize <- expression({
      ## initialize theta from raw observations..
       if (G$family$n.theta>0) {
         Theta <- c(-1, log(0.2*var(G$y)^.5))
         G$family$putTheta(Theta)
       } ## otherwise fixed theta supplied
    })

    initialize <- expression({
        if (any(is.na(y))) stop("NA values not allowed for the scaled t family")
        n <- rep(1, nobs)
        mustart <- y + (y == 0)*.1
    })
    postproc <- expression({
      object$family$family <- 
      paste("Scaled t(",paste(round(object$family$getTheta(TRUE),3),collapse=","),")",sep="")
    })
    rd <- function(mu,wt,scale) {
     ## simulate data given fitted latent variable in mu 
      theta <- get(".Theta")
      nu <- exp(theta[1])+2; sig <- exp(theta[2])
      n <- length(mu)
      stats::rt(n=n,df=nu)*sig + mu
    }

    environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
    environment(rd)<- environment(variance) <- environment(putTheta) <- env

    structure(list(family = "scaled t", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,variance=variance,postproc=postproc,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,ls=ls, preinitialize=preinitialize,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta,   
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta, rd=rd),
        class = c("extended.family","family"))
} ## scat



## zero inflated Poisson (Simon Wood)...

lind <- function(l,th,deriv=0,k=0) {
## evaluate th[1] + exp(th[2])*l and some derivs
  th[2] <- exp(th[2])
  r <- list(p = th[1] + (k+th[2])*l)
  r$p.l <- k + th[2]   ## p_l
  r$p.ll <- 0 ## p_ll
  if (deriv) {
    n <- length(l);  
    r$p.lllth <- r$p.llth <- r$p.lth <- r$p.th <- matrix(0,n,2)
    r$p.th[,1] <- 1   ## dp/dth1 
    r$p.th[,2] <- th[2]*l ## dp/dth2 
    r$p.lth[,2] <- th[2] ## p_lth2
    r$p.llll <- r$p.lll <- 0   ## p_lll,p_llll
    r$p.llth2 <- r$p.lth2 <- r$p.th2 <- matrix(0,n,3) ## ordered l_th1th1,l_th1th2,l_th2th2
    r$p.th2[,3] <- l*th[2] ## p_th2th2
    r$p.lth2[,3] <- th[2]  ## p_lth2th2
  }
  r
} ## lind

logid <- function(l,th,deriv=0,a=0,trans=TRUE) {
## evaluate exp(th[1]+th[2]*l)/(1+exp(th[1]+th[2]*l))
## and some of its derivatives
## if trans==TRUE then it is assumed that the 
## transformation th[2] = exp(th[2]) is applied on input
  b <- 1-2*a
  ## x is dth[2]/dth[2]' where th[2]' is input version, xx is second deriv over first
  if (trans) { xx <- 1; x <- th[2] <- exp(th[2])} else { x <- 1;xx <- 0}
  p <- f <- th[1] + th[2] * l
  ind <- f > 0; ef <- exp(f[!ind])
  p[!ind] <- ef/(1+ef); p[ind] <- 1/(1+exp(-f[ind]))  
  r <- list(p = a + b * p)
  
  a1 <- p*(1-p); a2 <- p*(p*(2*p-3)+1)
  r$p.l <- b * th[2]*a1;   ## p_l
  r$p.ll <- b * th[2]^2*a2 ## p_ll
  if (deriv>0) { 
    n <- length(l); r$p.lth <- r$p.th <- matrix(0,n,2) 
    r$p.th[,1] <- b * a1   ## dp/dth1
    r$p.th[,2] <- b * l*a1 * x ## dp/dth2
    r$p.lth[,1] <- b * th[2]*a2  ## p_lth1
    r$p.lth[,2] <- b * (l*th[2]*a2 + a1) * x ## p_lth2

    a3 <- p*(p*(p*(-6*p + 12) -7)+1)
    r$p.lll <- b * th[2]^3*a3   ## p_lll
    r$p.llth <- matrix(0,n,2) 
    r$p.llth[,1] <- b * th[2]^2 * a3 ## p_llth1
    r$p.llth[,2] <- b * (l*th[2]^2*a3 + 2*th[2]*a2) * x ## p_ppth2

    a4 <- p*(p*(p*(p*(p*24-60)+50)-15)+1)
    r$p.llll <- b * th[2]^4*a4 ## p_llll
    r$p.lllth <- matrix(0,n,2)
    r$p.lllth[,1] <- b * th[2]^3*a4 ## p_lllth1
    r$p.lllth[,2] <- b * (th[2]^3*l*a4 + 3*th[2]^2*a3) * x ## p_lllth2

    r$p.llth2 <- r$p.lth2 <- r$p.th2 <- matrix(0,n,3) ## ordered l_th1th1,l_th1th2,l_th2th2

    r$p.th2[,1] <- b * a2   ## p_th1th1
    r$p.th2[,2] <- b * l*a2 * x ## p_th1th2  
    r$p.th2[,3] <- b * l*l*a2 * x * x + xx* r$p.th[,2] ## p_th2th2

    r$p.lth2[,1] <- b * th[2]*a3  ## p_lth1th1
    r$p.lth2[,2] <- b * (th[2]*l*a3 + a2) * x ## p_lth1th2  
    r$p.lth2[,3] <- b * (l*l*a3*th[2] + 2*l*a2) *x * x + xx*r$p.lth[,2]  ## p_lth2th2

    r$p.llth2[,1] <- b * th[2]^2*a4  ## p_llth1th1
    r$p.llth2[,2] <- b * (th[2]^2*l*a4 + 2*th[2]*a3) *x ## p_llth1th2  
    r$p.llth2[,3] <- b * (l*l*th[2]^2*a4 + 4*l*th[2]*a3 + 2*a2) *x*x + xx*r$p.llth[,2] ## p_llth2th2
  }
  r
} ## logid



ziP <- function (theta = NULL, link = "identity",b=0) { 
## zero inflated Poisson parameterized in terms of the log Poisson parameter, gamma. 
## eta = theta[1] + exp(theta[2])*gamma), and 1-p = exp(-exp(eta)) where p is 
## probability of presence.

  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("identity")) { 
    stats <- make.link(linktemp)
  } else  stop(linktemp, " link not available for zero inflated; available link for `lambda' is only  \"loga\"")
  ## Theta <-  NULL;
  n.theta <- 2
  if (!is.null(theta)) {
      ## fixed theta supplied
      iniTheta <-  c(theta[1],theta[2])
      n.theta <- 0 ## no thetas to estimate
  } else iniTheta <- c(0,0) ## inital theta value - start at Poisson
  
  env <- new.env(parent = environment(ziP))# new.env(parent = .GlobalEnv) 
  
  if (b<0) b <- 0; assign(".b", b, envir = env)
  assign(".Theta", iniTheta, envir = env)
  getTheta <- function(trans=FALSE) { 
  ## trans transforms to the original scale...
    th <- get(".Theta")
    if (trans) {
      th[2] <- get(".b") + exp(th[2])
    }
    th
  }

  putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))
  
  validmu <- function(mu) all(is.finite(mu))
  
  dev.resids <- function(y, mu, wt,theta=NULL) {
    ## this version ignores saturated likelihood
    if (is.null(theta)) theta <- get(".Theta")
    b <- get(".b")
    p <- theta[1] + (b + exp(theta[2])) * mu ## l.p. for prob present
    -2*zipll(y,mu,p,deriv=0)$l
  }
  
  Dd <- function(y, mu, theta, wt=NULL, level=0) {
    ## here mu is lin pred for Poisson mean so E(y) = exp(mu)
    ## Deviance for log lik of zero inflated Poisson. 
    ## code here is far more general than is needed - could deal 
    ## with any 2 parameter mapping of lp of mean to lp of prob presence.
    if (is.null(theta)) theta <- get(".Theta")
    deriv <- 1; if (level==1) deriv <- 2 else if (level>1) deriv <- 4 
    b <- get(".b")
    g <- lind(mu,theta,level,b) ## the derviatives of the transform mapping mu to p
    z <- zipll(y,mu,g$p,deriv)
    oo <- list();n <- length(y)
    if (is.null(wt)) wt <- rep(1,n)
    oo$Dmu <- -2*wt*(z$l1[,1] + z$l1[,2]*g$p.l)
    oo$Dmu2 <- -2*wt*(z$l2[,1] + 2*z$l2[,2]*g$p.l + z$l2[,3]*g$p.l^2 + z$l1[,2]*g$p.ll)
    ## WARNING: following requires z$El1 term to be added if transform modified so 
    ##          that g$p.ll != 0....
    oo$EDmu2 <- -2*wt*(z$El2[,1] + 2*z$El2[,2]*g$p.l + z$El2[,3]*g$p.l^2)

    if (level>0) { ## l,p - ll,lp,pp -  lll,llp,lpp,ppp - llll,lllp,llpp,lppp,pppp
      oo$Dth <- -2*wt*z$l1[,2]*g$p.th ## l_p p_th
      oo$Dmuth <- -2*wt*(z$l2[,2]*g$p.th + z$l2[,3]*g$p.l*g$p.th + z$l1[,2]*g$p.lth) 
      oo$Dmu2th <- -2*wt*(z$l3[,2]*g$p.th + 2*z$l3[,3]*g$p.l*g$p.th + 2* z$l2[,2]*g$p.lth + 
       z$l3[,4]*g$p.l^2*g$p.th + z$l2[,3]*(2*g$p.l*g$p.lth + g$p.th*g$p.ll) + z$l1[,2]*g$p.llth)
      oo$Dmu3 <- -2*wt*(z$l3[,1] + 3*z$l3[,2]*g$p.l + 3*z$l3[,3]*g$p.l^2 + 3*z$l2[,2]*g$p.ll +
       z$l3[,4]*g$p.l^3 +3*z$l2[,3]*g$p.l*g$p.ll + z$l1[,2]*g$p.lll)
    } 
    if (level>1) {
      p.thth <- matrix(0,n,3);p.thth[,1] <- g$p.th[,1]^2
      p.thth[,2] <- g$p.th[,1]*g$p.th[,2];p.thth[,3] <- g$p.th[,2]^2
      oo$Dth2 <- -2*wt*(z$l2[,3]*p.thth + z$l1[,2]*g$p.th2)
      p.lthth <- matrix(0,n,3);p.lthth[,1] <- g$p.th[,1]*g$p.lth[,1]*2
      p.lthth[,2] <- g$p.th[,1]*g$p.lth[,2] + g$p.th[,2]*g$p.lth[,1];
      p.lthth[,3] <- g$p.th[,2]*g$p.lth[,2]*2
      oo$Dmuth2 <- -2*wt*( z$l3[,3]*p.thth + z$l2[,2]*g$p.th2 + z$l3[,4]*g$p.l*p.thth + 
        z$l2[,3]*(g$p.th2*g$p.l + p.lthth) + z$l1[,2]*g$p.lth2)
      p.lthlth <- matrix(0,n,3);p.lthlth[,1] <- g$p.lth[,1]*g$p.lth[,1]*2
      p.lthlth[,2] <- g$p.lth[,1]*g$p.lth[,2] + g$p.lth[,2]*g$p.lth[,1];
      p.lthlth[,3] <- g$p.lth[,2]*g$p.lth[,2]*2
      p.llthth <- matrix(0,n,3);p.llthth[,1] <- g$p.th[,1]*g$p.llth[,1]*2
      p.llthth[,2] <- g$p.th[,1]*g$p.llth[,2] + g$p.th[,2]*g$p.llth[,1];
      p.llthth[,3] <- g$p.th[,2]*g$p.llth[,2]*2

      oo$Dmu2th2 <- -2*wt*(z$l4[,3]*p.thth + z$l3[,2]*g$p.th2 + 2*z$l4[,4] * p.thth *g$p.l + 2*z$l3[,3]*(g$p.th2*g$p.l + p.lthth) +
        2*z$l2[,2]*g$p.lth2 + z$l4[,5]*p.thth*g$p.l^2 + z$l3[,4]*(g$p.th2*g$p.l^2 + 2*p.lthth*g$p.l + p.thth*g$p.ll) +
        z$l2[,3]*(p.lthlth + 2*g$p.l*g$p.lth2 + p.llthth + g$p.th2*g$p.ll) + z$l1[,2]*g$p.llth2)

      oo$Dmu3th <- -2*wt*(z$l4[,2]*g$p.th + 3*z$l4[,3]*g$p.th*g$p.l + 3*z$l3[,2]*g$p.lth + 2*z$l4[,4]*g$p.th*g$p.l^2 + 
        z$l3[,3]*(6*g$p.lth*g$p.l + 3*g$p.th*g$p.ll) + 3*z$l2[,2]*g$p.llth + z$l4[,4]*g$p.th*g$p.l^2 + 
        z$l4[,5]*g$p.th*g$p.l^3 + 3*z$l3[,4]*(g$p.l^2*g$p.lth + g$p.th*g$p.l*g$p.ll) +
        z$l2[,3]*(3*g$p.lth*g$p.ll + 3*g$p.l*g$p.llth + g$p.th*g$p.lll) + z$l1[,2]*g$p.lllth)

      oo$Dmu4 <- -2*wt*(z$l4[,1] + 4*z$l4[,2]*g$p.l + 6*z$l4[,3]*g$p.l^2 + 6*z$l3[,2]*g$p.ll + 
        4*z$l4[,4]*g$p.l^3 + 12*z$l3[,3]*g$p.l*g$p.ll + 4*z$l2[,2]*g$p.lll + z$l4[,5] * g$p.l^4 +
        6*z$l3[,4]*g$p.l^2*g$p.ll + z$l2[,3] *(4*g$p.l*g$p.lll + 3*g$p.ll^2) + z$l1[,2]*g$p.llll)

    }
    oo
  } ## end Dd for ziP
  
  aic <- function(y, mu, theta=NULL, wt, dev) {
    if (is.null(theta)) theta <- get(".Theta")
    b <- get(".b")
    p <- theta[1] + (b+ exp(theta[2])) * mu ## l.p. for prob present
    sum(-2*wt*zipll(y,mu,p,0)$l)
  }

  ls <- function(y,w,n,theta,scale) {
       ## the log saturated likelihood function.
       ## ls is defined as zero for REML/ML expression as deviance is defined as -2*log.lik 
       list(ls=0,## saturated log likelihood
            lsth1=c(0,0),  ## first deriv vector w.r.t theta - last element relates to scale
            lsth2=matrix(0,2,2)) ##Hessian w.r.t. theta
  }


    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the zero inflated Poisson family")
        if (all.equal(y,round(y))!=TRUE) {
          stop("Non-integer response variables are not allowed with ziP ")
        }
        if ((min(y)==0&&max(y)==1)) stop("Using ziP for binary data makes no sense")
        n <- rep(1, nobs)
        mustart <- log(y + (y==0)/5) 
    })

    postproc <- expression({
      object$family$family <- 
      paste("Zero inflated Poisson(",paste(round(object$family$getTheta(TRUE),3),collapse=","),")",sep="")
      ## need to fix deviance here!!
      ## wts <- object$prior.weights
      lf <- object$family$saturated.ll(G$y,family, object$prior.weights)
      ## storing the saturated loglik for each datum...
      object$family$data <- list(ls = lf)   
      l2 <- object$family$dev.resids(G$y,object$linear.predictors,object$prior.weights)
      object$deviance <- sum(l2-lf)
      fnull <- function(gamma,object) {
        ## evaluate deviance for single parameter model
        sum(object$family$dev.resids(object$y, rep(gamma,length(object$y)), object$prior.weights))
      }
      meany <- mean(object$y)
      object$null.deviance <- optimize(fnull,interval=c(meany/5,meany*3),object=object)$objective - sum(lf)
 
      ## object$weights <- pmax(0,object$working.weights) ## Fisher can be too extreme
      ## E(y) = p * E(y) - but really can't mess with fitted.values if e.g. rd is to work.

    })

#   fv <- function(lp,theta=NULL) {
#    ## optional function to give fitted values... 
#      if (is.null(theta)) theta <- get(".Theta")
#      th1 <- theta[1]; th2 <- exp(theta[2]); 
#      eta <- th1 + th2*lp
#      p <- 1 - exp(-exp(eta))
#      fv <- lambda <- exp(lp)
#      ind <- lp < log(.Machine$double.eps)/2
#      fv[!ind] <- p[!ind] * lambda[!ind]/(1-exp(-lambda[!ind]))
#      fv[ind] <- p[ind]
#      fv      
#    } ## fv

    rd <- function(mu,wt,scale) {
    ## simulate data given fitted latent variable in mu 
      rzip <- function(gamma,theta) { ## generate ziP deviates according to model and lp gamma
        y <- gamma; n <- length(y)
        lambda <- exp(gamma)
        mlam <- max(c(lambda[is.finite(lambda)],.Machine$double.eps^.2))
        lambda[!is.finite(lambda)] <- mlam
        b <- get(".b")
        eta <- theta[1] + (b+exp(theta[2]))*gamma
        p <- 1- exp(-exp(eta))
        ind <- p > runif(n)
        y[!ind] <- 0
        #np <- sum(ind)
        ## generate from zero truncated Poisson, given presence...
        lami <- lambda[ind]
        yi <- p0 <- dpois(0,lami)
        nearly1 <- 1 - .Machine$double.eps*10
        ii <- p0 > nearly1 
        yi[ii] <- 1 ## lambda so low that almost certainly y=1
        yi[!ii] <- qpois(runif(sum(!ii),p0[!ii],nearly1),lami[!ii])
        y[ind] <- yi 
        y
      } 
      rzip(mu,get(".Theta"))
    }
   
   saturated.ll <- function(y,family,wt=rep(1,length(y))) {
      ## function to get saturated ll for ziP - 
      ## actually computes -2 sat ll.
      pind <- y>0 ## only these are interesting
      wt <- wt[pind]
      y <- y[pind];
      mu <- log(y)  
      keep.on <- TRUE
      theta <- family$getTheta() 
      r <- family$Dd(y,mu,theta,wt)
      l <- family$dev.resids(y,mu,wt,theta)
      lmax <- max(abs(l))
      ucov <- abs(r$Dmu) > lmax*1e-7
      k <- 0
      while (keep.on) { 
        step <- -r$Dmu/r$Dmu2
        step[!ucov] <- 0
        mu1 <- mu + step
        l1 <- family$dev.resids(y,mu1,wt,theta)
        ind <- l1>l & ucov
        kk <- 0
        while (sum(ind)>0&&kk<50) {
            step[ind] <- step[ind]/2
            mu1 <- mu + step
            l1 <- family$dev.resids(y,mu1,wt,theta) 
            ind <- l1>l & ucov
            kk <- kk + 1
        }       
        mu <- mu1;l <- l1
        r <- family$Dd(y,mu,theta,wt)
        ucov <- abs(r$Dmu) > lmax*1e-7
        k <- k + 1
        if (all(!ucov)||k==100) keep.on <- FALSE
      }
      l1 <- rep(0,length(pind));l1[pind] <- l
      l1
    } ## saturated.ll

  residuals <- function(object,type=c("deviance","working","response")) {
    if (type == "working") { 
      res <- object$residuals 
    } else if (type == "response") {
      res <- object$y - predict.gam(object,type="response")
    } else if (type == "deviance") { 
      y <- object$y
      mu <- object$linear.predictors
      wts <- object$prior.weights
      res <- object$family$dev.resids(y,mu,wts)
      res <- res - object$family$saturated.ll(y,object$family,wts)
      fv <- predict.gam(object,type="response")
      s <- attr(res,"sign")
      if (is.null(s)) s <- sign(y-fv)
      res <- as.numeric(sqrt(pmax(res,0)) * s) 
    }
    res
  } ## residuals

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                beta=NULL,off=NULL,Vb=NULL) {
  ## optional function to give predicted values - idea is that 
  ## predict.gam(...,type="response") will use this, and that
  ## either eta will be provided, or {X, beta, off, Vb}. family$data
  ## contains any family specific extra information. 
 
    theta <- family$getTheta()

    if (is.null(eta)) { ## return probabilities
      gamma <- drop(X%*%beta + off) ## linear predictor for poisson parameter 
      se <- if (se) drop(sqrt(pmax(0,rowSums((X%*%Vb)*X)))) else NULL ## se of lin pred
    } else { se <- NULL; gamma <- eta}
    ## now compute linear predictor for probability of presence...
    b <- get(".b")
    eta <- theta[1] + (b+exp(theta[2]))*gamma
    et <- exp(eta)
    mu <- p <- 1 - exp(-et)
    fv <- lambda <- exp(gamma)  
    ind <- gamma < log(.Machine$double.eps)/2
    mu[!ind] <- lambda[!ind]/(1-exp(-lambda[!ind]))
    mu[ind] <- 1
    fv <- list(p*mu)    ## E(y)    
    if (is.null(se)) return(fv) else {
      dp.dg <- p  
      ind <- eta < log(.Machine$double.xmax)/2
      dp.dg[!ind] <- 0
      dp.dg <- exp(-et)*et*exp(theta[2])
      dmu.dg <- (lambda + 1)*mu - mu^2
      fv[[2]] <- abs(dp.dg*mu+dmu.dg*p)*se   
      names(fv) <- c("fit","se.fit")
      return(fv)
    }
  } ## predict


   
  environment(saturated.ll) <- environment(dev.resids) <- environment(Dd) <-
  environment(aic) <- environment(getTheta) <- environment(rd) <- environment(predict) <-
  environment(putTheta) <- env

  structure(list(family = "zero inflated Poisson", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd, rd=rd,residuals=residuals,
        aic = aic, mu.eta = stats$mu.eta, g2g = stats$g2g,g3g=stats$g3g, g4g=stats$g4g, 
        #preinitialize=preinitialize,
        initialize = initialize,postproc=postproc,ls=ls,no.r.sq=TRUE,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta,predict=predict,
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,saturated.ll = saturated.ll),
        class = c("extended.family","family"))

} ## ziP


