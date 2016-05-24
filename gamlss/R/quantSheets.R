# Paul Eilers, Mikis Stasinopoulos, Bob Rigby, Vlasios Voudouris 
# there are one main functions here 
# quantileSheets() 
#-------------------------------------------------------------------------------
# TO DO
# i) define residuals for quantileSheets using the a) van Buuren code in code and data file (at the desk top) or create the flexDist() and use it for calculation the cdf of y 
# ii) the fitted values should be for all the x_is not only for the one used in the fitting
# iii) (print OK), (fitted is OK),  (predict is OK)  
#      (residuals OK)  (z.scores OK) 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
quantSheets <- function(y, x,
                    x.lambda = 1,
                    p.lambda = 1,
                        data = NULL,                
                        cent = 100 * pnorm((-4:4) * 2/3),
                     control = quantSheets.control(...),
                       print = TRUE,
                     ...      
                    )
{
#-------------------------------------------------------------------------------
# local functions
#-------------------------------------------------------------------------------
  tpower <- function(x, t, p = 1)
    # Truncated p-th power function
    (x - t) ^ p * (x > t)
#-------------------------------------------------------------------------------
  bbase <- function(x, xl = min(x), xr = max(x), ndx = 10, deg = 3){
    # Construct B-spline basis
    dx <- (xr - xl) / ndx
    kts <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    P <- outer(x, kts, FUN = tpower, deg)
    n <- dim(P)[2]
    D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
    B <- (-1) ^ (deg + 1) * P %*% t(D)
    B 
  } 
#-------------------------------------------------------------------------------
  rowtens = function(X){
    # Row-wise tensor products
    one <- matrix(1, nrow = 1, ncol = ncol(X))
    kronecker(X, one) * kronecker(one, X)
  }
#-------------------------------------------------------------------------------
   ptrans <- function(x, p) if (abs(p)<=0.0001) log(x) else I(x^p)
#   ptrans <- function(x, p) if (p==0) log(x) else x^p
invptrans <- function(x, p) if (abs(p)<=0.0001) exp(x) else x^(1/p)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
   scall <- deparse(sys.call())
    ylab <- deparse(substitute(y))
    xlab <- deparse(substitute(x))
       y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
       x <- if (!is.null(data)) get(deparse(substitute(x)), envir=as.environment(data)) else x
#  if (!is.null(data)) {attach(data); on.exit(detach(data))}
  if (!is.null(control$power))
    {
      ox <- x
       x <-  ptrans(x,control$power)
    } 
       m <- length(x)
      xl <- min(x)
      xr <- max(x)
   nsegx <- control$x.inter # this has to be an argument  
   nsegp <- control$p.inter  # this has to be an argument
    bdeg <- control$degree
       p <- cent/100  #seq(0.05, 0.95, by = 0.1) 
       n <- length(p)
      Bx <- bbase(x, xl, xr, nsegx, bdeg)  ## basis for X
   if (control$logit)                      ## basis for p
   {
     logitp <- log(p/(1-p))
         Bp <- bbase(logitp, -20, 20, nsegp, bdeg)
   }
  else
   {
         Bp <- bbase(p, 0, 1, nsegp, bdeg)  
  }
     nbx <- ncol(Bx)
     nbp <- ncol(Bp)
      Tx <- rowtens(Bx)
      Tp <- rowtens(Bp)
# -----     the   penalties   -------------------------------------------------      
      Dx <- diff(diag(nbx), diff = control$order)
      Dp <- diff(diag(nbp), diff = control$order)
      Px <- x.lambda * t(Dx) %*% Dx
      Pp <- p.lambda * t(Dp) %*% Dp
       P <- kronecker(Pp, diag(nbx)) + kronecker(diag(nbp), Px)
 #  kappa <- 0 # what is kappa is it ?? 
       P <- P + control$kappa * diag(nrow(P))
# Initialize
       Y <- outer(y, rep(1, n))
       Z <- 0 * Y + mean(Y)
      OP <- outer(rep(1, m), p)
# Iterate
       b <- 0.001
  for (it in 1:control$n.cyc) 
  {
       R <- Y - Z                 
       W <- ifelse(R > 0, OP, 1- OP) / sqrt(b + R ^ 2)
       Q <- t(Tx) %*% W %*% Tp
  dim(Q) <- c(nbx, nbx, nbp, nbp)
       Q <- aperm(Q, c(1, 3, 2, 4))
  dim(Q) <- c(nbx * nbp, nbx * nbp)
       r <- t(Bx) %*% (Y * W) %*% Bp
  dim(r) <- c(nbx * nbp, 1)
       A <- solve(Q + P, r)
  dim(A) <- c(nbx, nbp)
    Znew <- Bx %*% A %*% t(Bp)
      dz <- sum(abs(Z - Znew))
    # cat(it, dz, '\n')
  if (dz < control$c.crit) break
       Z <- Znew
  }
      xg <- seq(xl, xr, length = 100)
      Bg <- bbase(xg, xl, xr, nsegx, bdeg)
      Zg <- Bg %*% A %*% t(Bp)
if (!is.null(control$power))
    { 
        x <- ox
       xg <- invptrans(xg, control$power)
    }
 if (control$plot)
  {
       plot(x, y, pch = 15, cex = 0.5, col = gray(0.7), ylab=ylab, xlab=xlab)
   matlines(x[order(x)], Z[order(x),], type = 'l', lty = 1, lwd = 1)
#matplot(xg, Zg, type = 'l', lty = 1, lwd = 2)
  } 
#  per[ii] <- (1-sum(oyvar>ll)/length(oyvar))*100
#  if (!save) cat("% of cases below ", var,"centile is ", per[ii], "\n" )
  colnames(Zg) <- as.character(round(cent,2))
# calculating percentages of the sample
  per <- rep(0, length(cent))
  quantFun <- list()
  for (i in 1:length(cent))
  {
    quantFun[[i]] <- splinefun(xg,Zg[,i], method="natural")  
              ll <- quantFun[[i]](x)
          per[i] <- (1-sum(y>ll)/length(y))*100
 if (print)   cat("% of cases below ", cent[i],"centile is ", per[i], "\n" )
  }
  names(quantFun) <- namesFun <- as.character(round(cent,2)) # put less digits
  out <- list(y=y, x=x, knots=xg, fitted.values=Zg, cent=cent, sample.perc=per, quantFun=quantFun, call=scall, ylab=ylab, xlab=xlab, namesFun=namesFun, noObs=length(y))  
  class(out) <- "quantSheets"
return(invisible(out))
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
quantSheets.control <- function( x.inter = 10,  p.inter = 10,  degree= 3, 
    logit = FALSE, order = 2, kappa = 0,  n.cyc = 100, c.crit = 1e-5,  
    plot = TRUE , power = NULL, ...)
{ 
# x.inter: number of intervals in x-variable
# p.inter: number of intrrvals in the probability direction
# degree : degree of the polynomial in B-spline
# logit : whether to use logit(p) instead of p (probabilities)
# order : the order of the penalty   
# kappa : is a ridge parameter set to zero (for no ridge efferct) 
# n.cyc :  number of cycles 
# c.crit :  convergence criterion
# plot ; whether to plot the resultsE
  if(x.inter <= 0) {
    warning("the value of x.inter supplied is less than 0, the value of 10 was used instead")
    x.inter <- 10 }  
  if(p.inter <= 0) {
    warning("the value of p.inter supplied is less than 0, the value of 10 was used instead")
    p.inter <- 10 }
  if(degree <= 0) {
    warning("the value of degree supplied is less than zero or negative the default value of 3 was used instead")  
    degree <- 3}
  if(order < 0) {
    warning("the value of order supplied is zero or negative the default value of 2 was used instead")
    order <- 2}
  if(kappa < 0) {
    warning("the value of kapa supplied is less than 0, the value of zero was used instead")
    kappa <- 0 }
  if(n.cyc < 0) {
    warning("the value of n.cyc is less than zero the default value of 100 was used instead")
    n.cyc  <- 100}              
  
  if(c.crit <= 0) {
    warning("the value of c.crit is less or equal than zero the default value of 1e-5 was used instead")
    c.crit <- 1e-5}                        
 out <- list(x.inter = x.inter, p.inter = p.inter, degree = degree,  logit = as.logical(logit)[1], order = order, kappa=kappa, n.cyc = n.cyc, c.crit = c.crit, plot = as.logical(plot)[1], power=power)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# methods for quantSheets
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
print.quantSheets <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{  
  cat("\nQuantile Sheets fit of", x$xlab, "against",x$ylab,   "\n") 
 
  # cat("\nCall: ", deparse(x$call), "\n\n")
  cat("\nCall: ", deparse(x$call),  "\n", fill=TRUE)
  cat("Estimation is done on the following centile points", "\n", fill=TRUE)
  cat(x$namesFun, "\n")
  invisible(x)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# fitted()
# the x argument produce fitted values at x and knots arguments produce fitted
# values at the evaluating function min(x): max(x), 100 length
fitted.quantSheets <- function(object, type=c( "x", "knots"), cent=NULL,...)
{
  type <- match.arg(type) 
  switch(type, 
"x"={ 
     if (!is.null(cent))
     {
       if (as.character(cent) %in% object$namesFun ) 
       {
         mat<- object$quantFun[[as.character(cent)]](object$x) ; return(mat)
       }
       else stop("the values in cent is not in the estimated quantile list")     
     }
     else
     {
      mat <- matrix(0,  nrow=length(object$x), ncol=length(object$cent))
      for (i in 1:length(object$cent))
      {
        mat[,i] <- object$quantFun[[i]](object$x)    
      }
      colnames(mat) <- object$namesFun
      return(mat)
     }
  }, 
"knots"={
        if (!is.null(cent)) 
        {
          if (as.character(cent) %in% object$namesFun) 
          {
            mat<-  object$fitted.values[, as.character(cent)] ; return(mat)
          }
          else stop("the values in cent is not in the estimated quantile list")
        }
       else 
       {
         mat <- object$fitted.values
         return(mat)
       }
 })  
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# predict()
predict.quantSheets <- function(object, 
                        newdata = NULL,
                        ...
)
{
  # what <- match.arg(what)
  # if newdata is not null check the class
  # if data.frame check for  object$xlab
  # other if vector use it 
  if (!is.null(newdata))
  {
    if (class(newdata)=="data.frame")
    {
      if (!object$xlab%in%names(newdata)) stop("the name in the data.frame do not much the x-variable in the model")
      x <- newdata[[object$xlab]]
    } else 
    {
      x <- newdata  
    }
    out <- matrix(0, nrow = length(x),  ncol= length(object$cent)) 
    colnames(out) <- object$namesFun
    for (i in 1:length(object$cent))
    {
      out[,i] <-  object$quantFun[[i]](x)
    }   
  } 
if (is.null(newdata)) out <- fitted(object) 
  out    
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
residuals.quantSheets <- function(object, inter=100, all=FALSE, ...)
{
  if (all==FALSE)
  {
    xinter <- seq(min(object$x), max(object$x), length=inter)
    fx <- findInterval(object$x, quantile(object$x, seq(0,1, length=101)), 
                       rightmost.closed = T)
    ymin <- tapply(object$y, fx, "min")
    ymax <- tapply(object$y, fx, "max")
    tol <- (ymax-ymin)*.10
    midInter <- tapply(object$x, fx, "mean")  
    predMat <- predict(object, newdata=midInter)
    minboth <- apply(cbind(ymin,predMat[,1]),1,"min")
    maxboth <- apply(cbind(ymax,predMat[,dim(predMat)[2]]),1,"max")
    Ires <- rep(0, length(object$y))
    for (i in 1:inter)
    {
      #  if (i==4) browser()
      FDIST <- flexDist(quantiles=list(values=predMat[i,], prob=(object$cent/100)), 
                        plot=FALSE, lower=minboth[i]-tol[i], upper=maxboth[i]+tol[i])
      Ires[fx==i] <- FDIST$pFun(object$y[fx==i])  
    }
    Ires <- ifelse(Ires> 0.999999999, 0.999999999, Ires)
    Ires <- ifelse(Ires<=0.000000001,0.000000001, Ires)
    res <- qnorm(Ires) 
  } else
  {
    res <-   z.scoresQS(object, y=object$y, x=object$x)
  }   
  res  
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#===============================================================================
#===============================================================================
findPower <- function(y, x, data = NULL,  lim.trans = c(0, 1.5), prof=FALSE, k=2,  c.crit = 0.01, step=0.1)  
{
  ylab <- deparse(substitute(y))
  xlab <- deparse(substitute(x))
  y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
  x <- if (!is.null(data)) get(deparse(substitute(x)), envir=as.environment(data)) else x
  ##  checking for transformation in x        
  cat("*** Checking for transformation for x ***", "\n") 
  ptrans <- function(x, p) if (abs(p)<=0.0001) log(x) else I(x^p) 
  fn <- function(p) GAIC(gamlss(y~pb(ptrans(x,p)), c.crit = c.crit, trace=FALSE), k=k)
  if (prof) # profile dev
  {
    pp <- seq(lim.trans[1],lim.trans[2], step) 
    pdev <- rep(0, length(pp)) 
    for (i in 1:length(pp)) 
    {
      pdev[i] <- fn(pp[i])  
      #   cat(pp[i], pdev[i], "\n")
    }
    plot(pdev~pp, type="l")
    points(pdev~pp,col="blue")
    par <- pp[which.min(pdev)]
    cat('*** power parameters ', par,"***"," \n") 
  } else
  {
    #     fn <- function(p) GAIC(gamlss(y~pb(ptrans(x,p)),sigma.fo=~pb(ptrans(x,p)),nu.fo=~pb(ptrans(x,p)), data=data,tau.fo=~pb(ptrans(x,p)), c.crit = c.crit, trace=FALSE, family=BCT), k=k)
    #    fn <- function(p) GAIC(gamlss(y~pb(ptrans(x,p)),sigma.fo=~pb(ptrans(x,p)), data=data, c.crit = c.crit, trace=FALSE), k=k)
    fn <- function(p) GAIC(gamlss(y~pb(ptrans(x,p)),  c.crit = c.crit, trace=FALSE), k=k)
    par <- optimise(fn, lower=lim.trans[1], upper=lim.trans[2])$minimum
    # browser()
    #  par <- optim(.5, fn, lower=lim.trans[1], upper=lim.trans[2], method="L-BFGS-B")$par
    cat('*** power parameters ', par,"***"," \n") 
  }  
  par
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
z.scoresQS <- function(object, y, x, plot=FALSE, tol = NULL )
{
  if (!is(object,"quantSheets"))  stop(paste("This is not an quantSheets  object", "\n", "")) 
  if (is.null(y)) stop("the y values should be set for z-scores")
  if (is.null(x)) stop("the x values should be set for z-scores")
  if (length(y)!= length(x)) stop("length of x and y is not the same")
  tol <- if (is.null(tol))
  {
    if (length(y)==1) 1 else (max(y)-min(y))*0.1
  } else tol
  pred <- predict(object, newdata=x)
  rqres <- rep(0, length(x))
  for (i in 1:length(x))
  {
    #  cat(i, "\n")
    # if (i==267) browser()
    FDIST <- flexDist(quantiles=list(values=pred[i,], prob=(object$cent/100)), 
                      plot=plot) 
    if (plot) abline(v=y[i], col="blue")
    rqres[i] <- FDIST$pFun(y[i]) 
    #     if(rqres[i]==0)
    #     {
    #       FDIST <-  flexDist(quantiles=list(values=pred[i,], prob=(object$cent/100)), 
    #                plot=plot, lower=y[i]-tol)
    #    rqres[i] <- FDIST$pFun(y[i]) 
    #     }
    #    if(rqres[i]==1)
    #    {
    #      FDIST <-  flexDist(quantiles=list(values=pred[i,], prob=(object$cent/100)), 
    #                         plot=plot, upper=y[i]+tol)
    #      rqres[i] <- FDIST$pFun(y[i]) 
    #    }           
  }  
  rqres <- ifelse( rqres> 0.999999999, 0.999999999,  rqres)
  rqres <- ifelse( rqres<=0.000000001,0.000000001,  rqres)
  rqres <- qnorm(rqres)
  rqres  
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------