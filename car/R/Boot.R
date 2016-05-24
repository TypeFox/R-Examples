# Boot:  A reimplementation of bootCase using the 'boot' package to do the
#  work.  The main function 'Boot' creates the 'statistic' argument to
#  'boot', and passes this function to 'boot'
# For the call b1 <- Boot(m1) and b2 <- bootCase(m1),
#   b2 was the returned bootstaps; this is in b1$t
# b1 is of class c("Boot", "boot", so ALL the 'boot' generic methods work
# 'Boot' has new generic methods 'summary', 'confint' and 'hist'

# notes:  See Davison and Hinkley Chapters 6 and 7.
# Boot.lm, method="case" is the simple case resampling
#             method="residual" uses algorithm 6.3, p. 271
#               The use of weights comes from using 'pearson' residuals
#               This is equivalent to alg. 6.1, p262, unweighted
# Boot.glm method="case" as for lm
#             method="residual" not implemented.  Too problematic.
# May 23, 2012 Sanford Weisberg sandy@umn.edu
# June 1, 2012:  changed from class c("Boot", "boot") to just class "boot"
# 2012-12-10 replaced .GlobalEnv with .carEnv to avoid warnings
# 2013-07-08 changed .carEnv to car:::.carEnv so 'boot' could find the environment
# 4014-08-17: added calls to requireNamespace() and :: where necessary. J. Fox
# 2015-01-27 .carEnv now in global environment. John
# 2015-02-20: fixed coding error in Boot.nls(). John

Boot <- function(object, f, labels, R=999, method){UseMethod("Boot")}

Boot.default <- function(object, f=coef, labels=names(coef(object)),
                     R=999, method=c("case", "residual")) {
  if(!(requireNamespace("boot"))) stop("The 'boot' package is missing")
  f0 <- f(object)
  if(length(labels) != length(f0)) labels <- paste("V", seq(length(f0)), sep="")
  method <- match.arg(method)
  if(method=="case") {
     boot.f <- function(data, indices, .fn) {
      assign(".boot.indices", indices, envir=.carEnv)
      mod <- update(object, subset=get(".boot.indices", envir=.carEnv))
      if(mod$qr$rank != object$qr$rank){
            out <- .fn(object)
            out <- rep(NA, length(out)) } else  {out <- .fn(mod)}
     out
     }
    } else {
    boot.f <- function(data, indices, .fn) {
      first <- all(indices == seq(length(indices)))
      res <- if(first) object$residuals else
                  residuals(object, type="pearson")/sqrt(1 - hatvalues(object))
      res <- if(!first) (res - mean(res)) else res
      val <- fitted(object) + res[indices]
      if (!is.null(object$na.action)){
            pad <- object$na.action
            attr(pad, "class") <- "exclude" 
            val <- naresid(pad, val)
            }
      assign(".y.boot", val, envir=.carEnv)
      mod <- update(object, get(".y.boot", envir=.carEnv) ~ .)
      if(mod$qr$rank != object$qr$rank){
            out <- .fn(object)
            out <- rep(NA, length(out)) } else  {out <- .fn(mod)}
      out
      }
  } 
  b <- boot::boot(data.frame(update(object, model=TRUE)$model), boot.f, R, .fn=f)
  colnames(b$t) <- labels
  if(exists(".y.boot", envir=.carEnv))
     remove(".y.boot", envir=.carEnv)
  if(exists(".boot.indices", envir=.carEnv))
     remove(".boot.indices", envir=.carEnv)
  b
  }
  
Boot.lm <- function(object, f=coef, labels=names(coef(object)),
                     R=999, method=c("case", "residual"))
   Boot.default(object, f, labels, R, method)
  
Boot.glm <- function(object, f=coef, labels=names(coef(object)),
                     R=999, method=c("case", "residual")) {
  method <- match.arg(method)
  if(method=="case") { Boot.default(object, f, labels, R, method)
    } else {
    stop("Residual bootstrap not implemented in the 'car' function Boot.
  Use the 'boot' function in the 'boot' package to write
  your own version of residual bootstrap for a glm.")
   }
  }
  
Boot.nls <- function(object, f=coef, labels=names(coef(object)),
                     R=999, method=c("case", "residual")) {
  if(!(requireNamespace("boot"))) stop("The 'boot' package is missing")
  f0 <- f(object)
  if(length(labels) != length(f0)) labels <- paste("V", seq(length(f0)), sep="")
  method <- match.arg(method)
  opt<-options(show.error.messages = FALSE)
  if(method=="case") {
     boot.f <- function(data, indices, .fn) {
         assign(".boot.indices", indices, envir=.carEnv)
         mod <- try(update(object, subset=get(".boot.indices", envir=.carEnv),
                 start=coef(object)))
         if(class(mod) == "try-error"){
            out <- .fn(object)
            out <- rep(NA, length(out)) } else  {out <- .fn(mod)}
     out
     }
    } else {
    boot.f <- function(data, indices, .fn) {
      first <- all(indices == seq(length(indices)))
      res <- residuals(object)
      val <- fitted(object) + res[indices]
      if (!is.null(object$na.action)){
            pad <- object$na.action
            attr(pad, "class") <- "exclude"
            val <- naresid(pad, val)
            }
      assign(".y.boot", val, envir=.carEnv)
      mod <- try(update(object, get(".y.boot", envir=.carEnv) ~ .,
             start=coef(object)))
      if(class(mod) == "try-error"){
            out <- .fn(object)
            out <- rep(NA, length(out)) } else  {out <- .fn(mod)}
      out
      }
  }
  b <- boot::boot(data.frame(update(object, model=TRUE)$model), boot.f, R, .fn=f)
  colnames(b$t) <- labels
  if(exists(".y.boot", envir=.carEnv))
     remove(".y.boot", envir=.carEnv)
  if(exists(".boot.indices", envir=.carEnv))
     remove(".boot.indices", envir=.carEnv)
  options(opt)
  d <- dim(na.omit(b$t))[1]
  if(d != R)
   cat( paste("\n","Number of bootstraps was", d, "out of", R, "attempted", "\n"))
  b
  }

confint.boot <- function(object, parm, level = 0.95,
    type = c("bca", "norm", "basic", "perc", "all"), ...){
  if (!requireNamespace("boot")) "boot package is missing"
  cl <- match.call()
  type <- match.arg(type)
  if(type=="all") stop("Use 'boot.ci' if you want to see 'all' types")
  types <-   c("bca", "norm", "basic", "perc")
  typelab <- c("bca", "normal",   "basic", "percent")[match(type, types)]
  nn <- colnames(object$t)
  names(nn) <- nn
  parm <- if(missing(parm)) which(!is.na(object$t0)) else parm
  out <- list()
  for (j in 1:length(parm)){
   out[[j]] <- boot::boot.ci(object, conf=level, type=type, index=parm[j], ...)
  }
  levs <- unlist(lapply(level, function(x) c( (1-x)/2, 1 - (1-x)/2)))
  ints <- matrix(0, nrow=length(parm), ncol=length(levs))
  rownames(ints) <- nn[parm]
  for (j in 1:length(parm)){
    which <- if(typelab=="normal") 2:3 else 4:5
    ints[j, ] <- as.vector(t(out[[j]][[typelab]][, which]))
  }
  or <- order(levs)
  levs <- levs[or]
  ints <- ints[, or, drop=FALSE]
  colnames(ints) <- paste(round(100*levs, 1), " %",sep="")
  attr(ints,"type") <- typelab
  class(ints) <- c("confint.boot", class(ints))
  ints
}
print.confint.boot <- function(x, ...) {
  cat("Bootstrap quantiles, type = ", attr(x, "type"), "\n\n")
  print(as.data.frame(x), ...)
  }


summary.boot <- function (object, parm, high.moments = FALSE, 
    extremes=FALSE, ...)
{
    cl <- match.call()
    skew1 <- function(x){
       x <- x[!is.na(x)]
       xbar <- mean(x)
       sum((x-xbar)^3)/(length(x) * sd(x)^3)
       }
    kurtosis1 <- function (x) {
       x <- x[!is.na(x)]
       xbar <- mean(x)
       sum((x - xbar)^4)/(length(x) * sd(x)^4) - 3
       }
    not.aliased <-  !is.na(object$t0)
    boots <- object$t[ , not.aliased, drop=FALSE ]
    nc <- if(is.matrix(boots)) ncol(boots) else 1
    stats <- matrix(rep(NA, nc * 10), ncol = 10)
    rownames(stats) <- colnames(boots)
    stats[, 1] <- apply(boots, 2, function(x) sum(!is.na(x))) # num. obs
    stats[, 2] <- object$t0[not.aliased]  # point estimate
    stats[, 3] <- apply(boots, 2, function(x) mean(x, na.rm=TRUE)) - stats[, 2]
    stats[, 5] <- apply(boots, 2, function(x) median(x, na.rm=TRUE))
    stats[, 4] <- apply(boots, 2, function(x) sd(x, na.rm=TRUE))
    stats[, 6] <- apply(boots, 2, function(x) min(x, na.rm=TRUE))
    stats[, 7] <- apply(boots, 2, function(x) max(x, na.rm=TRUE))
    stats[, 8] <- stats[, 7] - stats[, 6]
    stats[, 9] <- apply(boots, 2, skew1)
    stats[, 10] <- apply(boots, 2, kurtosis1)
    colnames(stats) <- c(
      "R", "original", "bootBias", "bootSE", "bootMed", "bootMin",
     "bootMax", "bootRange", "bootSkew", "bootKurtosis")
    stats <- as.data.frame(stats)
    class(stats) <- c("summary.boot", "data.frame")
    use <- rep(TRUE, 10)
    if (high.moments == FALSE) use[9:10] <- FALSE
    if (extremes==FALSE) use[6:8] <- FALSE
    parm <- if(missing(parm)) 1:dim(stats)[1] else parm
    return(stats[parm , use])
}
print.summary.boot <- 
   function(x, digits = max(getOption("digits") - 2, 3), ...)
{
    print.data.frame(x, digits=digits, ...)
}

hist.boot <- function(x, parm, layout=NULL, ask, main="", freq=FALSE,
      estPoint = TRUE, point.col="black", point.lty=2, point.lwd=2,
      estDensity = !freq, den.col="blue", den.lty=1, den.lwd=2,
      estNormal = !freq,  nor.col="red",   nor.lty=2, nor.lwd=2,
      ci=c("bca", "none", "percentile"), level=0.95,
      legend=c("top", "none", "separate"), box=TRUE, ...){
  not.aliased <- which(!is.na(x$t0))
  ci <- match.arg(ci)
  legend <- match.arg(legend)
  pe <- x$t0[not.aliased]
  if(is.null(names(pe))) names(pe) <- colnames(x$t)
  if(missing(parm)) parm <- not.aliased
  nt <- length(parm) + if(legend == "separate") 1 else 0
  if (nt > 1 & (is.null(layout) || is.numeric(layout))) {
    if(is.null(layout)){
         layout <- switch(min(nt, 9), c(1, 1), c(1, 2), c(2, 2), c(2, 2),
                             c(3, 2), c(3, 2), c(3, 3), c(3, 3), c(3, 3))
    }
    ask <- if(missing(ask) || is.null(ask)) prod(layout) < nt else ask
    oma3 <- if(legend == "top") 0.5 + estPoint + estDensity + estNormal
              else 1.5
    op <- par(mfrow=layout, ask=ask, no.readonly=TRUE,
            oma=c(0, 0, oma3, 0), mar=c(5, 4, 1, 2) + .1)
    on.exit(par(op))
    }
  if(ci != "none") clim <- confint(x, method=ci, level=level)
  pn <- colnames(x$t)
  names(pn) <- pn
  what <- c(estNormal & !freq, estDensity & !freq, ci != "none", estPoint)
  for (j in parm){
# determine the range of the y-axis
       z <- na.omit(x$t[, j])
       h <- hist(z, plot=FALSE)
       d <- density(z)
       n <- pnorm(0)/(sd <- sd(z))
       m <- if(freq == FALSE) max(h$density, d$y, n) else max(h$counts)
       plot(h, xlab=pn[j], freq=freq,
            main=if(length(parm)==1) main else "", ylim=c(0, m), ...)
       if(estDensity & !freq){
          lines(d, col=den.col, lty=den.lty, lwd=den.lwd)
          }
       if(estNormal & !freq){
          z <- na.omit(x$t[, j])
          xx <- seq(-4, 4, length=400)
          xbar <- mean(z)
          sd <- sd(z)
          lines( xbar + sd*xx, dnorm(xx)/sd, col=nor.col, lty=nor.lty,
              lwd=nor.lwd)
          }
       if(ci != "none") lines( clim[j ,], c(0, 0), lwd=4)
       if(estPoint) abline(v=pe[j], lty=point.lty, col=point.col, lwd=point.lwd)
       if(box) box()
       if( j == parm[1] & legend == "top" ) { # add legend
		        usr <- par("usr")
		        legend.coords <- list(x=usr[1], y=usr[4] + 1.3 * (1 + sum(what)) *strheight("N"))
            legend( legend.coords,
             c("Normal Density", "Kernel Density",
             paste(ci, " ", round(100*level), "% CI", sep=""),
                   "Obs. Value")[what],
             lty=c(nor.lty, den.lty, 1, point.lty)[what],
             col=c(nor.col, den.col, "black", point.col)[what],
             fill=c(nor.col, den.col, "black", point.col)[what],
             lwd=c(2, 2, 4, 2)[what],
             border=c(nor.col, den.col, "black", point.col)[what],
             bty="n", cex=0.9, xpd=NA)#, #horiz=TRUE, offset=
		}
  }
  mtext(side=3, outer=TRUE, main, cex=1.2)
  if(legend == "separate") {
    plot(0:1, 0:1, xaxt="n", yaxt="n", xlab="", ylab="", type="n")
    use <- (1:4)[c( estNormal, estDensity, TRUE, ci != "none")]
    curves <- c("fitted normal density", "Kernel density est",
              paste(100*level, "% ", ci, " confidence interval", sep=""),
              "Observed value of statistic")
    colors <- c(nor.col, den.col, "black", point.col)
    lines <- c(nor.lty, den.lty, 1, point.lty)
    widths<- c(nor.lwd, den.lwd, 2, point.lty)
    legend("center", curves[use], lty=lines[use], lwd=widths[use],
          col=colors[use], box.col=par()$bg,
          title="Bootstrap histograms")
  }
  invisible(NULL)
  }
  

  
