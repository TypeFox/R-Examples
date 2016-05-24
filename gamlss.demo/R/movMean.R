Locmean <- function(y,x=seq(1,length(y)),w=rep(1,length(y)), span=.5)
 { 
          n <- length(y)  
        o.x <- order(x)
        x.o <- x[o.x]
        y.o <- y[o.x]
        w.o <- w[o.x]
        if (any(w.o<0)) stop("weights should be positive")
         fv <- rep(0,n)
          k <- trunc((n*span-1)/2)   
          S <- matrix(0, ncol=n, nrow=n)         
      for (i in 1:n)
        {
                  s <- rep(0,n)
                 lo <- max(i-k,1)
                 up <- min(i+k,n)
        s[c(lo:up)] <-1/(up-lo+1)
              S[i,] <- s
              fv[i] <-weighted.mean(window(y.o, start=lo, end=up), window(w.o, start=lo, end=up), na.rm = FALSE)
        }
         fv <- fv[o.x]
        rss <- sum((w*(y-fv))^2)
        out <- list(fitted.values=fv, residuals=y-fv[o.x], edf=sum(diag(S)), rss=rss, lambda=span, span=span, 
                       k=k, y=y, x=x, w=w) 
    class(out) <- "locW"
    out
 }
#----------------------------------------------------------------------------------------
Locpoly <- function(y,x=seq(1,length(y)),w=rep(1,length(y)), span=.5, order=1)
 { 
#-------------------------------------------------------------------------------  
  .hat.WX<-function (w,x) 
  {
    p <- length(x)
    X  <- if (!is.matrix(x)) matrix(cbind(1,x), ncol=2)
    else x 
    k <- length(w)
    p <- dim(X)[1]
    if (p != k) 
      stop("`w' and 'x' are not having the same length")
    Is <- sqrt(w)
    if (any(!is.finite(Is))) 
      warning("diagonal weights has non-finite entries")
    WX <- X 
    WX[] <- Is * X
    h<-hat(qr(WX))
    h
  }
#-------------------------------------------------------------------------------  
  
          n <- length(y)  
        o.x <- order(x)
        x.o <- x[o.x]
        y.o <- y[o.x]
        w.o <- w[o.x]
        if (any(w.o<0)) stop("weights should be positive")
         fv <- hatv<- rep(0,n)
          k <- trunc((n*span-1)/2)   
          S <- matrix(0, ncol=n, nrow=n)         
      for (i in 1:n)
        {               
                  s <- rep(0,n)
                 lo <- max(i-k,1)
                 up <- min(i+k,n)
        s[c(lo:up)] <-1
                fit <- lm(y.o~poly(x.o, order), weights=w*s)
            hatv[i] <- .hat.WX(w*s,x.o)[i]
              fv[i] <- fitted(fit)[i]
        }
         fv <- fv[o.x]
        rss <- sum((w*(y-fv))^2)
        out <- list(fitted.values=fv, residuals=y-fv, edf=sum(sum(hatv)), rss=rss, lambda=span, span=span, 
                     k=k, y=y, x=x, w=w)   
    class(out) <- "locW"
    out
 }
#----------------------------------------------------------------------------------------
WLocmean <- function(y,x=seq(1,length(y)),w=rep(1,length(y)), lambda=.5)
 { 
          n <- length(y)  
        o.x <- order(x)
        x.o <- x[o.x]
        y.o <- y[o.x]
        w.o <- w[o.x]
        if (any(w.o<0)) stop("weights should be positive")
         fv <- rep(0,n)
          S <- matrix(0, ncol=n, nrow=n)        
      for (i in 1:n)
        {
                  s <- dNO((x.o-x.o[i]), mu=0, sigma=lambda)
              S[i,] <- s/sum(s)
              fv[i] <-weighted.mean(y.o, w.o*s, na.rm = FALSE)
        }
      ## this function is for plotting specific rows of the S matrix  used in the book
      # PlotSrows(S,x.o)
         fv <- fv[o.x]
        rss <- sum((w*(y-fv))^2)
        out <- list(fitted.values=fv, residuals=y-fv, edf=sum(diag(S)), rss=rss, lambda=lambda, 
                     y=y, x=x, w=w)   
    class(out) <- "locW"
    out
 }
#----------------------------------------------------------------------------------------
PlotSrows <- function(S, x)
{
#pdf("C:/gamlss/book-2009/Chapter9/pictures/RowsOftheSmatrix.pdf")
#postscript("C:/gamlss/book-2009/Chapter9/pictures/RowsOftheSmatrix.ps")
 op <- par( mfrow=c(4,1), mar=par("mar")+c(0,1,0,0), pch="+", cex=0.45, cex.lab=1.8, cex.axis=1.6)
      page <- c("S[1,]", "S[40,]", "S[90,]", "S[125,]")
      j =1
      for (i in c(1, 40, 90, 125))
      {
       plot(S[i,]~x, type="h", xlab="times (ordered)", ylab=page[j], ylim=c(0, 0.30) )
      j=j+1
      }
      par(op)
#dev.off()
}
#op <- par(mfrow = c(3, 4),   mar=par("mar")+c(0,1,0,0),
#         pch="+", cex=0.45, cex.lab=1.8, cex.axis=1.6)
#page <- c("age^-0.5", "log(age)", "age^.5", "age")
#pcd4 <- c("cd4^-0.5", "log(cd4+0.1)", "cd4^.5")
#for (i in 1:3)
# {
# yy <- with(CD4, eval(parse(text=pcd4[i])))
# for (j in 1:4)
#   {
#   xx<- with(CD4, eval(parse(text=page[j])))
#     plot(yy~xx, xlab=page[j], ylab=pcd4[i] )
#    }
#  }
#par(op)
#----------------------------------------------------------------------------------------
WLocpoly <- function(y,x=seq(1,length(y)),w=rep(1,length(y)), lambda=.5, order=1)
 { 
  #-------------------------------------------------------------------------------  
  .hat.WX<-function (w,x) 
  {
    p <- length(x)
    X  <- if (!is.matrix(x)) matrix(cbind(1,x), ncol=2)
    else x 
    k <- length(w)
    p <- dim(X)[1]
    if (p != k) 
      stop("`w' and 'x' are not having the same length")
    Is <- sqrt(w)
    if (any(!is.finite(Is))) 
      warning("diagonal weights has non-finite entries")
    WX <- X 
    WX[] <- Is * X
    h<-hat(qr(WX))
    h
  }
  #------------------------------------------------------------------------------- 
          n <- length(y)  
        o.x <- order(x)
        x.o <- x[o.x]
        y.o <- y[o.x]
        w.o <- w[o.x]
        if (any(w.o<0)) stop("weights should be positive")
         fv <- hatv<- rep(0,n)
         # k <- trunc((n*lambda-1)/2)   
          S <- matrix(0, ncol=n, nrow=n)         
      for (i in 1:n)
        {               
                  s <- dNO((x.o-x.o[i]), mu=0, sigma=lambda)
                fit <- lm(y.o~poly(x.o, order), weights=w*s)
            hatv[i] <- .hat.WX(w*s,x.o)[i]
              fv[i] <- fitted(fit)[i]
        }
         fv <- fv[o.x]
        rss <- sum((w*(y-fv))^2)
        out <- list(fitted.values=fv, residuals=y-fv, edf=sum(hatv), rss=rss, lambda=lambda, 
                     y=y, x=x, w=w)   
    class(out) <- "locW"
    out
 }
#----------------------------------------------------------------------------------------
fitted.locW<-function(object,...) 
{
object$fitted.values
}
#------------------------------------------
residuals.locW<-function(object,...) 
{
object$residuals
}
#------------------------------------------
#---------------------------------------------------------------------------------------
