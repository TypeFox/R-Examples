

pk<- function(z, nx, ny = Inf){ # exaxt from R "stats"
   if(nx==Inf || ny==Inf){
      pp <- .C("pkolmogorov2x",
               p = as.double(z),
               as.integer(min(nx,ny)),
               PACKAGE="QTLRel")$p
   }else{
      pp <- .C("psmirnov2x",
               p = as.double(z),
               as.integer(nx), as.integer(ny),
               PACKAGE="QTLRel")$p
   }

   max(min(pp,1),0)
}

qk<- function(p, nx, ny = Inf){
   if(p<0 || p>1){
      stop("probability should be between 0 and 1.")
   }
   func<- function(x) pk(x, nx, ny)-p

   xx <- seq(-0.2,2.2,by=0.1)
   for(i in 1:length(xx)) if(func(xx[i]) >= 0) break
      xx <- seq(xx[i-1],xx[i+1],by=0.01)
      for(i in 1:length(xx)) if(func(xx[i]) >= 0) break

   qq <- uniroot(func, c(xx[i-1],xx[i+1]),
      tol = .Machine$double.eps, maxiter = 10000)$root
   max(min(qq,2),0)
}

pkolm <- function(z, nx, ny = Inf){
   z <- as.double(z)/sqrt(1/nx+1/ny)
   pp<- .C("kolm",
           z = as.double(z),
           as.integer(length(z)),
           PACKAGE="QTLRel")$z

   max(min(pp,1),0)
}

qkolm <- function(p, nx, ny = Inf){
   if(p<0 || p>1){
      stop("probability should be between 0 and 1.")
   }
   func <-  function(y)
     .C("kolm", y = as.double(y), as.integer(length(y)),
        PACKAGE="QTLRel")$y - p

   xx <- (-1:100)^3
   for(i in 1:length(xx)) if(func(xx[i]) >= 0) break
      xx <- (xx[i-1]:xx[i+1])
      for(i in 1:length(xx)) if(func(xx[i]) >= 0) break

   qq<- uniroot(func, c(xx[i-1],xx[i+1]),tol = .Machine$double.eps^0.35)$root
   max(min(qq * sqrt(1/nx+1/ny), 2), 0)
}

Fn.0 <- function(x){
   x <- sort(x)
   nx <- length(x)
   p <- c(1:nx)/nx
   data.frame(x=c(-Inf,x),p=c(0,p))
}

Fn. <- function(t,x){
   nt<- length(t)
   xm<- matrix(x,nrow=length(x),ncol=nt)
      xm<- sweep(xm,2,t,"<=")
   colMeans(xm)
}

qFn. <- function(t,x){
   if(any(t<0 | t>1))
      stop("Probability should between 0 and 1.")

   tmp<- Fn.0(x)
   nt<- length(t)
   pm <- matrix(tmp$p,nrow=length(tmp$p),ncol=nt)
      pm<- sweep(pm,2,t,">=")
   xm<- matrix(tmp$x,nrow=length(tmp$x),ncol=nt)
      xm[!pm]<- Inf
   qq <- apply(xm,2,min)
      qq[t==1] <- 1e+300
      qq[qq < -1e+300] <- -1e+300
   qq
}

Fn <- function(t,x){
   t<- as.double(t)
   x<- as.double(x)
   .C("Fn",
      t = as.double(t),
      as.integer(length(t)),
      as.double(x),
      as.integer(length(x)),
      PACKAGE="QTLRel")$t
}

qFn <- function(t,x){
   if(any(t<0 | t>1))
      stop("Probability should between 0 and 1.")

   t<- as.double(t)
   x<- as.double(x)
      x<- sort(x, decreasing = FALSE) # must sort to call "qFn"
   .C("qFn",
      t = as.double(t),
      as.integer(length(t)),
      as.double(x),
      as.integer(length(x)),
      PACKAGE="QTLRel")$t
}

qqPlot.default <- function(y, x = "norm", ...,
   xlab = if(is.numeric(x)) deparse(substitute(x)) else x,
   ylab = deparse(substitute(y)),
   main = "Q-Q Plot",
   type = "p", col = 1, lty = 2, lwd=1, pch=1, cex = 0.7,
   plot.it = TRUE, confidence = .95, qqline = c("observed", "expected", "none"))
{
   xlab<- xlab
   qqline<- match.arg(qqline)
   sy <- sort(y); nsy<- length(sy)
   if(!is.numeric(x)){
      q.function <- eval(parse(text=paste("q", x, sep="")))
      p.function <- eval(parse(text=paste("p", x, sep="")))
      px <- ppoints(nsy)
      x <- q.function(px, ...)
      sx <- sort(x); nsx <- Inf
   }else{
      sx <- sort(x); nsx <- length(x)
      if(nsy > nsx){
         sx <- approx(1L:nsx, sx, n = nsy)$y
      }else if(nsy < nsx){
         sy <- approx(1L:nsy, sy, n = nsx)$y
      }
      px<- ppoints(max(nsx,nsy))

      sx <- qFn(px, sx)
   }
   sy <- qFn(px, sy)

   if(confidence){
      if(confidence<0 || confidence>1)
         stop("confidence should be between 0 and 1")
      if(nsx*nsy > 5e+6){
         ka <- qkolm(confidence, nsx, nsy)
      }else ka<- qk(confidence, nsx, nsy) # exact
      pxL <- px - ka; pxL[pxL<0] <- 0
      pxU <- px + ka; pxU[pxU>1] <- 1
      if(nsx==Inf){
         qL <- q.function(pxL, ...); qL[qL == -Inf] <- -1e+38
         qU <- q.function(pxU, ...); qU[qU == Inf] <- 1e+38
      }else{
         qL <- qFn(pxL,sx)
         qU <- qFn(pxU,sx)
      }
   }

   if (qqline != "none"){
      Q.x <- quantile(sx, c(.25,.75))
      Q.y <- quantile(sy, c(.25,.75))
      b <- (Q.y[2] - Q.y[1])/(Q.x[2] - Q.x[1])
      a <- Q.y[1] - b*Q.x[1]
   }
   if (plot.it){
      if(confidence){
         plot(sx, sy, xlab = xlab, ylab = ylab, main = main,
            type = type, pch = pch, col = col, cex = cex, lty = lty, lwd = lwd)
         idx <- (sy < qL) | (sy> qU)
#         lines(sx[idx], sy[idx], type = type, pch = 4, cex = cex, lty = lty, lwd = lwd, col = col+1)
         points(sx[idx], sy[idx], pch = pch, cex = cex, col = col+1)
         if(qqline == "observed"){
            abline(a, b, col=col, lty = lty, lwd=lwd)
         }else if(qqline == "expected"){
            lines(sx, sx, col = col, lty = lty, lwd = lwd)
         }
         lines(sx, qL, col = col, lty = lty, lwd = lwd)
         lines(sx, qU, col = col, lty = lty, lwd = lwd)
      }else{
         plot(sx, sy, xlab = xlab, ylab = ylab, main = main,
            type = type, pch = pch, col = col, cex = cex, lty = lty, lwd = lwd)
         if(qqline == "expected" && !confidence){
            qqline<- "observed"
            warning("Observed qqline was drawn instead.")
         }
         if(qqline != "none") abline(a, b, col = col, lty = lty, lwd = lwd)
      }
   }
   if(confidence){
      invisible(list(x = sx, y = sy, lower = qL, upper = qU))
   }else
      invisible(list(x = sx, y = sy))
}

qqPlot <- function(y, x = "norm", ...,
   xlab = if(is.numeric(x)) deparse(substitute(x)) else x,
   ylab = deparse(substitute(y)),
   main = "Q-Q Plot",
   type = "p", col = 1, lty = 2, lwd=1, pch=1, cex = 0.7,
   plot.it = TRUE, confidence = .95, qqline = c("observed", "expected", "none"))
{
   UseMethod("qqPlot")
}
