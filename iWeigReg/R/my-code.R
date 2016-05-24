
#source R codes:

#histw()

#myinv()
#loglik()
#loglik.g()

#---
#mn.lik()
#mn.clik()

#mn.HT()
#mn.reg()
#mn.creg()

#---
#ate.lik()
#ate.clik()

#ate.HT()
#ate.reg()
#ate.creg()

############################################ basic

#histw(): This function plots a weighted histogram.

histw <- function(x, w, xaxis, xmin, xmax, ymax, bar=TRUE, add=FALSE, col="black", dens=TRUE) {
     #x: data vector
     #w: inverse weight
     #xaxis: vector of cut points
     #xmin,xmax: the range of x coordinate
     #ymax: the maximum of y coordinate

     #bar: bar plot (if TRUE) or line plot
     #add: if TRUE, the plot is added to an existing plot
     #col: color of lines
     #dens: if TRUE, the histogram has a total area of one

     nbin <- length(xaxis)
     xbin <- cut(x, breaks=xaxis, include.lowest=T, labels=1:(nbin-1))

     y <- tapply(w, xbin, sum)
     y[is.na(y)] <- 0
     y <- y/sum(w)
     if (dens) y <- y/ (xaxis[-1]-xaxis[-nbin])

    if (!add) {
     plot.new()
     plot.window(xlim=c(xmin,xmax), ylim=c(0,ymax))
     axis(1, pos=0)
     axis(2, pos=xmin)
    }

     if (bar==1) {
        rect(xaxis[-nbin], 0, xaxis[-1], y)
     }
     else {
        xval <- as.vector(rbind(xaxis[-nbin],xaxis[-1])) 
        yval <- as.vector(rbind(y,y))
        lines(c(min(xmin,xaxis[1]), xval, max(xmax,xaxis[length(xaxis)])), 
              c(0,yval,0), lty="11", lwd=2, col=col)
     }
     invisible()
}


#Used for regression
myinv <- function(A, type="solve") {
   if (type=="solve")
     solve(A)
   else if (type=="ginv")
     ginv(A, tol=sqrt(.Machine$double.eps))
}


#loglik(): This function uses a constraint matrix and an indicator vector 
#to compute log-likelihood, its gradient and its hessian matrix. 

loglik <- function(lam, tr, h) 
{
   #lam: argument
   #tr: non-missingness/treatment indicator
   #h: constraint matrix

   n <- length(tr)
   k <- dim(h)[2]

   w <- as.vector(h%*%c(1,lam))
   w[tr==0] <- 1-w[tr==0]

   if (!any(is.na(lam)) & sum(w>0)==n) {
     val <- -sum(log(w))/n

     grad <- h[,-1, drop=FALSE]/w
     grad[tr==0,] <- -grad[tr==0,]

     gradient <- -apply(grad,2,sum)/n
     hessian <- t(grad)%*%grad/n
   }else {
     val <- Inf

     gradient <- rep(NA, k-1)
     hessian <- matrix(NA, k-1, k-1)
   }

   list(value=val, gradient=gradient, hessian=hessian)
}


#loglik.g(): This function uses a constraint matrix and an indicator vector 
#to compute log-likelihood, its gradient and its hessian matrix for calibrated 
#(or double robust) likelihood estimator in Tan (2010), Biometrika. 

loglik.g <- function(lam, tr, h, pr, g) 
{
   #lam: argument
   #tr: non-missingness/treatment indicator
   #h: constraint matrix
   #pr: fitted propensity score
   #g: control variate matrix

   n <- length(tr)

   g <- cbind(g)   #g may consist of one column
   k <- dim(g)[2]

   w <- rep(1/2,n)
   w[tr==1] <- h[tr==1,]%*%c(1,lam)

   if (!any(is.na(lam)) & sum(w[tr==1]>0)+sum(tr==0)==n) {
     val <- -sum(log(w[tr==1])/(1-pr[tr==1]))/n + sum(g%*%lam)/n

     mgrad <- matrix(0,n,k)
     mgrad[tr==1,] <- g[tr==1,]/w[tr==1]

     grad <- mgrad-g

     gradient <- -apply(grad,2,sum)/n
     hessian <- t((1-pr)*mgrad)%*%(mgrad)/n
   }else {
     val <- Inf

     gradient <- rep(NA, k)
     hessian <- matrix(NA, k,k)
   }

   list(value=val, gradient=gradient, hessian=hessian)
}


############################################ missing-data setup

#mn.lik(): This function implements the non-calibrated (or non-doubly robust) likelihood estimator 
#in Tan (2006), JASA, and returns the population mean for missing-data setup.

mn.lik <- function(y, tr, p, g, X=NULL, evar=TRUE, inv="solve") 
{
   #y: outcome
   #tr: non-missingness indicator
   #p: fitted propensity score
   #g: control variate matrix
   #X: propensity score model matrix from logistic regression (=NULL if p is known or treated to be so)
   #evar: if yes, variance estimation
   #inv: type of matrix inversion

   n <- length(tr)

   h <- cbind(p, (1-p)*g)  
   if (!is.null(X))
      h <- cbind(h, p*(1-p)*X) 

   solv <- trust(loglik, rep(0,dim(h)[2]-1), rinit=1, rmax=100, iterlim=1000, 
                 tr=tr, h=h)   
   lam <- solv$argument
   norm <- max(abs(solv$gradient))
   conv <- solv$converged

   w <- as.vector(h%*%c(1,lam))

   mu <- sum(y[tr==1]/w[tr==1])/n

   #variance estimation
  if (!evar) {
   list(mu=mu, w=w,
        lam=lam, norm=norm, conv=conv)

  } else {
   y1 <- y*tr/w
   yc1 <- y1-mu

   z <- h[,-1 ,drop=FALSE]* (tr/w-(1-tr)/(1-w))
   #zero <- apply(z, 2, mean)
   #equal to -solv$gradient

   B <- myinv(t(z)%*%z/n, type=inv)
   b1 <- B%*%(t(z)%*%y1/n)

   psi <- yc1 - as.vector(z%*%b1)
   v <- mean( psi^2 )/n

   list(mu=mu, v=v, w=w,
        lam=lam, norm=norm, conv=conv)
  }
}


#mn.clik(): This function implements the calibrated (or doubly robust) likelihood estimator 
#in Tan (2010), Biometrika, and returns the population mean for missing-data setup.

mn.clik <- function(y, tr, p, g, X=NULL, evar=TRUE, inv="solve") 
{
   #y: outcome
   #tr: non-missingness indicator
   #p: fitted propensity score
   #g: control variate matrix
   #X: propensity score model matrix from logistic regression (=NULL if p is known or treated to be so)
   #evar: if yes, variance estimation
   #inv: type of matrix inversion

   n <- length(tr)

   g <- cbind(g)   #g may consist of one column

   if (is.null(X)) {
     h1 <- cbind(p, (1-p)*g)

   } else {
     #first step of the two-step procedure in Tan (2010)
     out <- mn.lik(y, tr, p, g, X)

     h1 <- cbind(out$w, (1-p)*g)
   }

   solv <- trust(loglik.g, rep(0,dim(g)[2]), rinit=1, rmax=100, iterlim=1000,
                 tr=tr, h=h1, pr=p, g=g)
   lam <- solv$argument
   norm <- max(abs(solv$gradient))
   conv <- solv$converged

   w <- as.vector(h1%*%c(1,lam))

   mu <- sum(y[tr==1]/w[tr==1])/n

   #variance estimation
  if (!evar) {
   list(mu=mu, w=w,
        lam=lam, norm=norm, conv=conv)

  } else {
   y1 <- y*tr/w
   yc1 <- y1-mu

   z1 <- g*(tr/w-1)
   #zero1 <- apply(z1, 2, mean)
   #equal to -solv$gradient

   ze1 <- h1[,-1, drop=FALSE]/w
   zm1 <- g*tr/w

   B1 <- myinv(t(ze1)%*%zm1/n, type=inv)
   b1 <- B1%*%(t(ze1)%*%y1/n)

   psi <- yc1 - as.vector(z1%*%b1)

   if (!is.null(X)) {
      h <- cbind(p, (1-p)*g, p*(1-p)*X)

      z <- h[,-1 ,drop=FALSE]* (tr/out$w-(1-tr)/(1-out$w))
      B <- myinv(t(z)%*%z/n, type=inv)

      ze <- h[,-1, drop=FALSE]/w

      psi <- psi -
             as.vector(z%*%B%*%(t(ze)%*%(y1-as.vector(zm1%*%b1))/n)) 
   }
   v <- mean( psi^2 )/n      

   list(mu=mu, v=v, w=w,
        lam=lam, norm=norm, conv=conv)
  }
}


#mn.HT(): This function implements the Horvitz-Thompson estimator,
#and returns the population mean for missing-data setup. 

mn.HT <- function(y, tr, p, X=NULL, bal=FALSE) {
   #y: outcome
   #tr: non-missingness indicator
   #p: fitted propensity score
   #X: propensity score model matrix from logistic regression (=NULL if p is known or treated to be so)
   #bal: if TRUE, the function is used for checking balance

   n <- length(tr)
  
   y <- cbind(y)
   y1 <- y*tr/p

   if (bal)
      y1 <- y1- y 

   mu <- apply(y1, 2, mean)
   yc1 <- t(t(y1)-mu)

   #variance estimation
   if (!is.null(X)) {
      s <- X*(tr-p)
      V <- #solve(t(s)%*%s/n)
           solve(t(X*p)%*%((1-p)*X)/n)   #misspecification

      yc1 <- yc1 - s %*% V %*% (t(s)%*%yc1/n)
   }

   v <- apply( yc1^2, 2, mean )/n 

   list(mu=mu, v=v)
}


#mn.reg(): This function implements the non-calibrated (or non-doubly robust) regression estimator,
#and returns the population mean for missing-data setup. 

mn.reg <- function(y, tr, p, g, X=NULL, evar=TRUE, inv="solve") {
   #y: outcome
   #tr: non-missingness indicator
   #p: fitted propensity score
   #g: control variate matrix
   #X: propensity score model matrix from logistic regression (=NULL if p is known or treated to be so)
   #evar: if yes, variance estimation
   #inv: type of matrix inversion
 
   n <- length(tr)
  
   g <- cbind(g)   #g may consist of one column

   #raw
   y1 <- y*tr/p
   mu.raw <- mean(y1)  
   yc1 <- y1-mu.raw   #centered

   #regression
   z <- g*(tr/p-1)
   zero <- apply(z, 2, mean)
   zc <- t(t(z)-zero)   #centered

   if (!is.null(X)) {
      s <- X*(tr-p)

      z <- cbind(z, s)
      zero <- c(zero, rep(0,dim(X)[2]))
      zc <- cbind(zc, s)  
   }

   B <- myinv(t(zc)%*%zc/n, type=inv)
   b1 <- B%*%(t(zc)%*%yc1/n)

   mu <- mu.raw - as.vector(t(zero)%*%b1)

   #variance estimation
  if (!evar) {
   list(mu=mu, 
        b=b1)

  } else {
   psi <- (yc1 - as.vector(zc%*%b1) -
           as.vector((zc*(yc1-as.vector(zc%*%b1)))%*%(B%*%zero)))

   v <- mean( psi^2 )/n      

   list(mu=mu, v=v, 
        b=b1)
  }
}


#mn.creg(): This function implements the calibrated (or doubly robust) regression estimator 
#in Tan (2006), JASA, and returns the population mean for missing-data setup. 

mn.creg <- function(y, tr, p, g, X=NULL, evar=TRUE, inv="solve") {
   #y: outcome
   #tr: non-missingness indicator
   #p: fitted propensity score
   #g: control variate matrix
   #X: propensity score model matrix from logistic regression (=NULL if p is known or treated to be so)
   #evar: if yes, variance estimation
   #inv: type of matrix inversion

   n <- length(tr)
    
   g <- cbind(g)   #g may consist of one column

   #raw
   y1 <- y*tr/p
   mu.raw <- mean(y1)  
   yc1 <- y1-mu.raw   #centered

   #regression
   z <- g*(tr/p-1)
   zero <- apply(z, 2, mean)
   zc <- t(t(z)-zero)   #centered

   zm1 <- g*tr/p

   if (!is.null(X)) {
      s <- X*(tr-p)
      sm1 <- X*tr

      z <- cbind(z, s)
      zero <- c(zero, rep(0,dim(X)[2]))
      zc <- cbind(zc, s)  

      zm1 <- cbind(zm1, sm1)
   }

   B1 <- myinv(t(z)%*%zm1/n, type=inv)
   b1 <- B1%*%(t(z)%*%y1/n)

   mu <- mu.raw - as.vector(t(zero)%*%b1)

   #variance estimation
  if (!evar) {
   list(mu=mu, 
        b=b1)

  } else {
   psi <- yc1 - as.vector(zc%*%b1) - 
          as.vector((z*(y1-as.vector(zm1%*%b1)))%*%(B1%*%zero)) 

   v <- mean( psi^2 )/n

   list(mu=mu, v=v,
        b=b1)
  }
}


############################################ causal-inference setup

#ate.lik(): This function implements the non-calibrated (or non-doubly robust) likelihood estimator 
#in Tan (2006), JASA, and returns the population means for causal-inference setup.

ate.lik <- function(y, tr, p, g0,g1, X=NULL, evar=TRUE, inv="solve") {
   #y: outcome
   #tr: treatment indicator
   #p: fitted propensity score
   #g0,g1: control variate matrix
   #X: propensity score model matrix from logistic regression (=NULL if p is known or treated to be so)
   #evar: if yes, variance estimation
   #inv: type of matrix inversion

   n <- length(tr)

   h <- cbind(p, (1-p)*g1, p*g0)  
   if (!is.null(X))
      h <- cbind(h, p*(1-p)*X) 

   solv <- trust(loglik, rep(0,dim(h)[2]-1), rinit=1, rmax=100, iterlim=1000, 
                 tr=tr, h=h)   
   lam <- solv$argument
   norm <- max(abs(solv$gradient))
   conv <- solv$converged

   w <- as.vector(h%*%c(1,lam))

   mu1 <- sum(y[tr==1]/w[tr==1])/n
   mu0 <- sum(y[tr==0]/(1-w[tr==0]))/n
   diff <- mu1 - mu0

   mu <- c(mu1, mu0)
   names(mu) <- c("treat 1", "treat 0")

   #variance estimation
  if (!evar) {
   list(mu=mu, diff=diff, 
        w=w,
        lam=lam, norm=norm, conv=conv)

  } else {
   y1 <- y*tr/w
   y0 <- y*(1-tr)/(1-w)

   yc1 <- y1-mu1
   yc0 <- y0-mu0

   z <- h[,-1 ,drop=FALSE]* (tr/w-(1-tr)/(1-w))

   B <- myinv(t(z)%*%z/n, type=inv)
   b1 <- B%*%(t(z)%*%y1/n)
   b0 <- B%*%(t(z)%*%y0/n)

   psi1 <- yc1 - as.vector(z%*%b1)
   psi0 <- yc0 - as.vector(z%*%b0)

   v1 <- mean( psi1^2 )/n      
   v0 <- mean( psi0^2 )/n      
   v.diff <- mean( (psi1-psi0)^2 )/n      

   v <- c(v1, v0)
   names(v) <- c("treat 1", "treat 0")

   list(mu=mu, diff=diff, 
        v=v, v.diff=v.diff,
        w=w,
        lam=lam, norm=norm, conv=conv)
  }
}


#ate.clik(): This function implements the calibrated (or doubly robust) likelihood estimator 
#in Tan (2010), Biometrika, and returns the population means for causal-inference setup.

ate.clik <- function(y, tr, p, g0,g1, X=NULL, evar=TRUE, inv="solve") 
{
   #y: outcome
   #tr: treatment indicator
   #p: fitted propensity score
   #g0,g1: control variate matrix
   #X: propensity score model matrix from logistic regression (=NULL if p is known or treated to be so)
   #evar: if yes, variance estimation
   #inv: type of matrix inversion

   n <- length(tr)

   g1 <- cbind(g1)   #g1 may consist of one column
   g0 <- cbind(g0)

   #point estimation 1
   if (is.null(X)) {
     h1 <- cbind(p, (1-p)*g1)

   } else { 
     #first step of the two-step procedure in Tan (2010)
     out <- ate.lik(y, tr, p, g0,g1, X)

     h1 <- cbind(out$w, (1-p)*g1)
   }

   solv1 <- trust(loglik.g, rep(0,dim(g1)[2]), rinit=1, rmax=100, iterlim=1000,
                 tr=tr, h=h1, pr=p, g=g1)
   lam1 <- solv1$argument
   norm1 <- max(abs(solv1$gradient))
   conv1 <- solv1$converged

   w1 <- as.vector(h1%*%c(1,lam1))
   mu1 <- sum(y[tr==1]/w1[tr==1])/n

   #point estimation 0
   if (is.null(X)) {
     h0 <- cbind(1-p, p*g0)

   } else { 
     h0 <- cbind(1-out$w, p*g0)
   }

   solv0 <- trust(loglik.g, rep(0,dim(g0)[2]), rinit=1, rmax=100, iterlim=1000,
                 tr=1-tr, h=h0, pr=1-p, g=g0)
   lam0 <- solv0$argument
   norm0 <- max(abs(solv0$gradient))
   conv0 <- solv0$converged

   w0 <- as.vector(h0%*%c(1,lam0))
   mu0 <- sum(y[tr==0]/w0[tr==0])/n
   diff <- mu1 - mu0

   mu <- c(mu1, mu0)
   names(mu) <- c("treat 1", "treat 0")
  
  if (!evar) {
   list(mu=mu, diff=diff, 
        w=cbind(w1, w0),
        lam=c(lam1, lam0), norm=c(norm1, norm0), conv=c(conv1, conv0))

  } else {
   #variance estimation 1
   y1 <- y*tr/w1
   yc1 <- y*tr/w1-mu1

   z1 <- g1*(tr/w1-1)

   ze1 <- h1[,-1, drop=FALSE]/w1
   zm1 <- g1*tr/w1

   B1 <- myinv(t(ze1)%*%zm1/n, type=inv)
   b1 <- B1%*%(t(ze1)%*%y1/n)

   psi1 <- yc1 - as.vector(z1%*%b1)

   if (!is.null(X)) {
      h <- cbind(p, (1-p)*g1, p*g0, p*(1-p)*X)

      z <- h[,-1 ,drop=FALSE]* (tr/out$w-(1-tr)/(1-out$w))
      B <- myinv(t(z)%*%z/n, type=inv)

      ze <- h[,-1, drop=FALSE]/w1

      psi1 <- psi1 -
             as.vector(z%*%B%*%(t(ze)%*%(y1-as.vector(zm1%*%b1))/n)) 
   }

   v1 <- mean( psi1^2 )/n      

   #variance estimation 0
   y0 <- y*(1-tr)/w0
   yc0 <- y*(1-tr)/w0-mu0

   z0 <- g0*((1-tr)/w0-1)

   ze0 <- h0[,-1, drop=FALSE]/w0
   zm0 <- g0*(1-tr)/w0

   B0 <- myinv(t(ze0)%*%zm0/n, type=inv)
   b0 <- B0%*%(t(ze0)%*%y0/n)

   psi0 <- yc0 - as.vector(z0%*%b0)

   if (!is.null(X)) {
      z <- -z

      ze <- h[,-1, drop=FALSE]/w0

      psi0 <- psi0 -
             as.vector(z%*%B%*%(t(ze)%*%(y0-as.vector(zm0%*%b0))/n)) 
   }

   v0 <- mean( psi0^2 )/n      
   v.diff <- mean( (psi1-psi0)^2 )/n      

   v <- c(v1, v0)
   names(v) <- c("treat 1", "treat 0")

   list(mu=mu, diff=diff, 
        v=v, v.diff=v.diff,
        w=cbind(w1, w0),
        lam=c(lam1, lam0), norm=c(norm1, norm0), conv=c(conv1, conv0))
  }
}


#ate.HT(): This function implements the Horvitz-Thompson estimator,
#and returns the population means for causal-inference setup. 

ate.HT <- function(y, tr, p, X=NULL, bal=FALSE) {
   #y: outcome
   #tr: non-missingness indicator
   #p: fitted propensity score
   #X: propensity score model matrix from logistic regression (=NULL if p is known or treated to be so)
   #bal: if TRUE, the function is used for checking balance

   n <- length(tr)

   y <- cbind(y)
   y1 <- y*tr/p
   y0 <- y*(1-tr)/(1-p)

   if (bal) {
      y1 <- y1- y 
      y0 <- y0- y 
   }

   mu1 <- apply(y1, 2, mean)
   mu0 <- apply(y0, 2, mean)
   diff <- mu1 - mu0

   yc1 <- t(t(y1)-mu1)
   yc0 <- t(t(y0)-mu0)

   #variance estimation
   if (!is.null(X)) {
      s <- X*(tr-p)
      V <- #solve(t(s)%*%s/n)
           solve(t(X*p)%*%((1-p)*X)/n)   #misspecification

      yc1 <- yc1 - s %*% V %*% (t(s)%*%yc1/n)
      yc0 <- yc0 - s %*% V %*% (t(s)%*%yc0/n)
   }

   v1 <- apply( yc1^2, 2, mean )/n 
   v0 <- apply( yc0^2, 2, mean )/n 
   v.diff <- apply( (yc1-yc0)^2, 2, mean )/n 

   mu <- cbind(mu1, mu0)
   v <- cbind(v1, v0)
 
   list(mu=mu, diff=diff,
        v=v, v.diff=v.diff)
}


#ate.reg(): This function implements the non-calibrated (or non-doubly robust) regression estimator,
#and returns the population means for causal-inference setup. 

ate.reg <- function(y, tr, p, g0,g1, X=NULL, evar=TRUE, inv="solve") {
   #y: outcome
   #tr: treatment indicator
   #p: fitted propensity score
   #g0,g1: control variate matrix
   #X: propensity score model matrix from logistic regression (=NULL if p is known or treated to be so)
   #evar: if yes, variance estimation
   #inv: type of matrix inversion

   n <- length(tr)
  
   #raw
   y1 <- y*tr/p
   y0 <- y*(1-tr)/(1-p)

   mu1.raw <- mean(y1)
   mu0.raw <- mean(y0)

   yc1 <- y1-mu1.raw   #centered
   yc0 <- y0-mu0.raw

   #regression
   z <- cbind(g1*(tr/p-1), g0*((1-tr)/(1-p)-1))
   zero <- apply(z, 2, mean)
   zc <- t(t(z)-zero)   #centered

   if (!is.null(X)) {
      s <- X*(tr-p)

      z <- cbind(z, s)
      zero <- c(zero, rep(0,dim(X)[2]))
      zc <- cbind(zc, s)  
   }

   B <- myinv(t(zc)%*%zc/n, type=inv)
   b1 <- B%*%(t(zc)%*%yc1/n)
   b0 <- B%*%(t(zc)%*%yc0/n)

   mu1 <- mu1.raw - as.vector(t(zero)%*%b1)
   mu0 <- mu0.raw - as.vector(t(zero)%*%b0)
   diff <- mu1 - mu0

   mu <- c(mu1, mu0)
   names(mu) <- c("treat 1", "treat 0")

   #variance estimation
  if (!evar) {
   list(mu=mu, diff=diff,
        b=cbind(b1, b0))

  } else {
   psi1 <- yc1 - as.vector(zc%*%b1) - 
           as.vector((zc*(yc1-as.vector(zc%*%b1)))%*%(B%*%zero))
   psi0 <- yc0 - as.vector(zc%*%b0) - 
           as.vector((zc*(yc0-as.vector(zc%*%b0)))%*%(B%*%zero))

   v1 <- mean( psi1^2 )/n      
   v0 <- mean( psi0^2 )/n      
   v.diff <- mean( (psi1-psi0)^2 )/n      

   v <- c(v1, v0)
   names(v) <- c("treat 1", "treat 0")

   list(mu=mu, diff=diff,
        v=v, v.diff=v.diff, 
        b=cbind(b1, b0))
  }
}


#ate.creg(): This function implements the calibrated (or doubly robust) regression estimator 
#in Tan (2006), JASA, and returns the population means for causal-inference setup. 

ate.creg <- function(y, tr, p, g0,g1, X=NULL, evar=TRUE, inv="solve") {
   #y: outcome
   #tr: treatment indicator
   #p: fitted propensity score
   #g0,g1: control variate matrix
   #X: propensity score model matrix from logistic regression (=NULL if p is known or treated to be so)
   #evar: if yes, variance estimation
   #inv: type of matrix inversion

   n <- length(tr)
  
   #raw
   y1 <- y*tr/p
   y0 <- y*(1-tr)/(1-p)

   mu1.raw <- mean(y1)
   mu0.raw <- mean(y0)

   yc1 <- y1-mu1.raw   #centered
   yc0 <- y0-mu0.raw

   #regression
   z <- cbind(g1*(tr/p-1), g0*((1-tr)/(1-p)-1))
   zero <- apply(z, 2, mean)
   zc <- t(t(z)-zero)   #centered

   zm1 <- cbind(g1/p, -g0/(1-p))*tr
   zm0 <- cbind(-g1/p, g0/(1-p))*(1-tr)

   if (!is.null(X)) {
      s <- X*(tr-p)
      sm1 <- X*tr
      sm0 <- -X*(1-tr)

      z <- cbind(z, s)
      zero <- c(zero, rep(0,dim(X)[2]))
      zc <- cbind(zc, s)  

      zm1 <- cbind(zm1, sm1)
      zm0 <- cbind(zm0, sm0)
   }

   B1 <- myinv(t(z)%*%zm1/n, type=inv)
   b1 <- B1%*%(t(z)%*%y1/n)

   B0 <- myinv(t(z)%*%zm0/n, type=inv)
   b0 <- B0%*%(t(z)%*%y0/n)

   mu1 <- mu1.raw - as.vector(t(zero)%*%b1)
   mu0 <- mu0.raw - as.vector(t(zero)%*%b0)
   diff <- mu1 - mu0

   mu <- c(mu1, mu0)
   names(mu) <- c("treat 1", "treat 0")
 
   #variance estimation
  if (!evar) {
   list(mu=mu, diff=diff, 
        b=cbind(b1, b0))

  } else {
   psi1 <- yc1 - as.vector(zc%*%b1) - 
           as.vector((z*(y1-as.vector(zm1%*%b1)))%*%(B1%*%zero))
   psi0 <- yc0 - as.vector(zc%*%b0) - 
           as.vector((z*(y0-as.vector(zm0%*%b0)))%*%(B0%*%zero))

   v1 <- mean( psi1^2 )/n      
   v0 <- mean( psi0^2 )/n      
   v.diff <- mean( (psi1-psi0)^2 )/n      

   v <- c(v1, v0)
   names(v) <- c("treat 1", "treat 0")

   list(mu=mu, diff=diff, 
        v=v, v.diff=v.diff,
        b=cbind(b1, b0))
  }
}

#save(list=ls(), file="iWeigReg.R")