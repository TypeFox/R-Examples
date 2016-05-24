#
# This is version 1.6 of the reldist functions.
#   by Mark S. Handcock
#
#   July 1, 2004
#
#  See the  README file for details.
#
#  This file contains R functions used to construct figures in
#
#  "Relative Distribution Methods in the Social Sciences" 
#  by Mark S. Handcock and Martina Morris, 
#  Springer-Verlag, 1999, ISBN 0387987789 
#
#  We make no claims that the code given
#  is as efficient as it could be, and we would be interested to learn
#  of more efficient approaches.
#
#  Other software is available from  the Relative Distribution Website:
#
#  http://www.stat.ucla.edu/~handcock/RelDist
#
#  The website contains related software, descriptions of the variables
#  and data file formats.
#
# relative distribution and density functions 
#
"reldist" <- function(y, yo=FALSE, ywgt=FALSE,yowgt=FALSE,
  show="none", decomp="locadd",
  location="median", scale="IQR",
  rpmult=FALSE, 
  z=FALSE, zo=FALSE,
  smooth = 0.35, 
  quiet = TRUE, 
  cdfplot=FALSE,
  ci=FALSE,
  bar="no",
  add=FALSE,
  graph=TRUE, type="l",
  xlab="Reference proportion",ylab="Relative Density",yaxs="r",
  yolabs=pretty(yo), yolabslabs=NULL,
  ylabs=pretty(y), ylabslabs=NULL, 
  yolabsloc=0.6, ylabsloc=1, 
  ylim=NULL, cex=0.8, lty=1,
  binn=100,
  aicc=seq(0.0001, 5, length=30),
  deciles=(0:10)/10,
  discrete=FALSE,
  method="bgk",
  y0=NULL,
  ...) {
#
# missing test
#
   missargs <- NULL
   if(missing(y)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yo)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(ywgt)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yowgt)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(show)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(decomp)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(location)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(scale)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(rpmult)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)} 
   if(missing(z)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(zo)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(smooth)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(quiet)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)} 
   if(missing(cdfplot)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(ci)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(bar)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(add)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(graph)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(type)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(xlab)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(ylab)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yaxs)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yolabs)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(ylabslabs)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yolabslabs)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(ylabs)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yolabsloc)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(ylabsloc)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)} 
   if(missing(ylim)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(cex)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(lty)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(binn)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(aicc)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(deciles)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(method)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   names(missargs) <- c(
   "y",
   "yo",
   "ywgt",
   "yowgt",
   "show",
   "decomp",
   "location",
   "scale",
   "rpmult", 
   "z",
   "zo",
   "smooth",
   "quiet", 
   "cdfplot",
   "ci",
   "bar",
   "add",
   "graph",
   "type",
   "xlab",
   "ylab",
   "yaxs",
   "yolabs",
   "yolabslabs",
   "ylabs",
   "ylabslabs",
   "yolabsloc",
   "ylabsloc", 
   "ylim",
   "cex",
   "lty",
   "binn",
   "aicc",
   "deciles",
   "method")
 #
   if(!missing(y0)){
     warning("You passed y0. Did you mean yo?")
     if(missing(yo)){
       yo <- y0	      
       missargs["yo"] <- FALSE      
       warning("y0 is being used as the sample from the reference distribution (in place of yo).")
     }
   }
   if(discrete){
    out <- rddist(
     y=y,
     yo=yo,
     ywgt=ywgt,
     yowgt=yowgt,
     show=show,
     decomp=decomp,
     location=location, scale=scale,
     rpmult=rpmult, 
     z=z, zo=zo,
     smooth = smooth,
     quiet = quiet, 
     cdfplot=cdfplot,
     ci=ci,
     bar=bar,
     add=add,
     graph=graph, type=type,
     xlab=xlab,ylab=ylab,yaxs=yaxs,
     yolabs=yolabs, yolabslabs=yolabslabs,
     ylabs=yolabs, ylabslabs=ylabslabs,
     yolabsloc=yolabsloc, ylabsloc=ylabsloc, 
     ylim=ylim, cex=cex, lty=lty,
     binn=binn,
     aicc=aicc,
     deciles=deciles,
     method=method, missargs=missargs,
     ...)
   }else{
    out <- rcdist(
     y=y,
     yo=yo,
     ywgt=ywgt,
     yowgt=yowgt,
     show=show,
     decomp=decomp,
     location=location, scale=scale,
     rpmult=rpmult, 
     z=z, zo=zo,
     smooth = smooth,
     quiet = quiet, 
     cdfplot=cdfplot,
     ci=ci,
     bar=bar,
     add=add,
     graph=graph, type=type,
     xlab=xlab,ylab=ylab,yaxs=yaxs,
     yolabs=yolabs, yolabslabs=yolabslabs, ylabs=ylabs,
     yolabsloc=yolabsloc, ylabsloc=ylabsloc, 
     ylim=ylim, cex=cex, lty=lty,
     binn=binn,
     aicc=aicc,
     deciles=deciles,
     method=method, missargs=missargs,
     ...)
   }
   if(!quiet){
    out
   }else{
    invisible(out)
   }
 }
 "rcdist" <- function(y, yo=FALSE, ywgt=FALSE,yowgt=FALSE,
          show="none", decomp="locadd",
          location="median", scale="IQR",
          rpmult=FALSE, 
          z=FALSE, zo=FALSE,
          smooth = 0.35, 
          quiet = TRUE, 
          cdfplot=FALSE,
          ci=FALSE,
          bar="no",
          add=FALSE,
          graph=TRUE, type="l",
          xlab="Reference proportion",ylab="Relative Density",yaxs="r",
          yolabs=pretty(yo), yolabslabs=NULL, ylabs=pretty(y),
          yolabsloc=0.6, ylabsloc=1, 
          ylim=NULL, cex=0.8, lty=1,
          binn=100,
          aicc=seq(0.00001, 5, length=30),
          deciles=(0:10)/10,
          method="quick", missargs,
          ...) {
 #
 #	# INPUT variables
 #	# required variables
 #	y	sample from comparison distribution
 #	yo	sample from reference distribution
 
 #	# smoothness variables
 #       smooth  degree of smoothness required in the fit. Higher values lead
 #               to smoother curves, lower to closer fits to the observed 
 #               data.
 #               If the local-linear binning estimator is used it is the span
 #               of the smoother as a proportion of [0,1].
 #               If the logspline estimator is used it is the number
 #               of knots used in the spline fit. If it is not specified
 #               the  number of knots  that  minimizes the BIC is used.
 #                 
 
 #       # plotting variables
 #	graph	graph the results on the current device
 #	bar	graph the deciles on the current device
 #               "no"   no deciles plotted
 #               "yes"  deciles plotted with the non-paramteric fit
 #               "only" plot the deciles only and not the non-parametric fit.
 #	add	add the density to the current plot?
 #	ylim	plotting limits for the vertical axis
 #	lty	line type to use for the density
 #	xlab	horizontal label
 #	ylab	vertical label
 #	ylabs	locations for labels to be added to the right axis
 #	ylabslabs	labels indicating the original scale
 #	ylabsloc distance of labels to right of axis (in lines)
 #	yolabs	locations for labels to be added to the top axis
 #	yolabslabs	labels indicating the original scale
 #	yolabsloc distance of labels above axis (in lines)
 #	yaxs	style of vertical axis
 
 #	# options
 #	cdfplot	calculate and plot the CDF rather than the density.
 #       ci      include 95% pointwise confidence bands on the plot?
 #	ywgt	weights on the comparison sample
 #	yowgt	weights on the reference sample
 #
 #       show    type of relative distribution to produce
 #                 none       comparison to reference
 #                 residual   location-matched reference to reference
 #                 effect     comparison to location-matched reference
 #       decomp  form of matching to the comparison sample
 #                 locmult     multiplitively scale the reference
 #                 locadd      additively shift the reference
 #                 lsadd       location/scale additive shift 
 #                 covariate   covariate adjust the reference
 #                             (requires z and zo to be specified)
 #      location how to measure location: "mean", "median"
 #       scale   how to measure scale: "standev", "IQR"
 #	rpmult  in calculation of polarization indices ONLY, multiplicatively
 #                scale the reference sample to 
 #                the comparison sample before comparing the two distributions?
 
 #	# debugging options
 #	binn	number of bins used in the smoother
 #	yaxs	style of vertical axis
 #
 #	# OUTPUT list components
 #	x	horizontal ordinates for the density (typically percent)
 #	y	density at x
 #	rp	95% confidence interval for the median relative polarization
 #                as lower bound, estimate, upper bound.
 #	rpl	lower relative polarization
 #	  	95% confidence interval for the lower relative polarization
 #                as lower bound, estimate, upper bound.
 #	rpu	upper relative polarization
 #	  	95% confidence interval for the upper relative polarization
 #                as lower bound, estimate, upper bound.
 #	cdf	list for the CDF
 #                x horizontal ordinates for the CDF (typically percent)
 #                y CDF at x
 #
 #
 #	NOTE: Most of the code is for the plotting and tinkering.
 #	      The guts of the method are forming the relative data
 #	      at the top. The rest is a standard fixed interval
 #	      density estimation with a few bells and whistles.
 #
 # -------------------------------------------------------------------------
 #
 #  if(decomp=="covariate"){
 #    stop("This has not been implemented yet. sorry. MSH")
 #  }
   if(is.logical(bar)){
    if(bar){bar <- "yes"}else{bar <- "no"}
   }
   if(really.missing(yo,missargs)){
    warning("The reference distribution, yo, was not passed. We will presume that y corresponds to a sample from a relative distribution directly.")
    m <- length(y)
    if(really.missing(ywgt,missargs)){ywgt <- rep(1/m,length=m)}
    y <- as.vector(y)
    x <- sort.list(y)
    ywgt <- as.vector(ywgt)[x]
    missargs[ "ywgt"] <- FALSE
    x <- y[x]
    # Force to have support part of [0,1]
    if(min(x,na.rm=TRUE) < 0 | max(x,na.rm=TRUE) > 1){
      x <- (x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
    }
    m <- length(x)
    n <- 0
    rmd <- list(x=x,wgt=ywgt,n=n,m=m)
   }else{
    rmd <- rmdata(y=as.vector(y), yo=as.vector(yo),
                  ywgt=as.vector(ywgt), yowgt=as.vector(yowgt),
                  show=show, decomp=decomp,
                  location=location, scale=scale,
                  z=as.vector(z), zo=as.vector(zo),
                  missargs=missargs
                 )
    y <- rmd$y
    yo <- rmd$yo
    ywgt <- rmd$ywgt
    yowgt <- rmd$yowgt
    missargs[    "y"] <- FALSE
    missargs[   "yo"] <- FALSE
    missargs[ "ywgt"] <- FALSE
    missargs["yowgt"] <- FALSE
    x <- rmd$x
    n <- rmd$n
    m <- rmd$m
   }
   r <- seq(0, 1, length = binn + 1)[-1] - 0.5/binn
 #
 # Now work on the data
 #
 #   In R use binning and local polynomial methods
 #
   if(really.missing(ywgt,missargs)){
 #   xx <- table(cut(x, breaks = seq(0, 1, length = binn + 1),
 #               include.lowest=TRUE,labels=FALSE))/m
     xx <- hist(x,breaks=seq(0, 1, length = binn + 1),
                 plot=FALSE,include.lowest = TRUE)$counts/m
   }else{
     ywgt <- ywgt/sum(ywgt)
     xxx <- cut(x, breaks = seq(0, 1, length = binn + 1),
                include.lowest=TRUE,labels=FALSE)
     xx <- rep(0,length=binn)
     for (i in unique(xxx)){    
       xx[i] <- sum(ywgt[xxx==i],na.rm=TRUE)
     }
   }
   xx <- as.vector(xx)
   gpdf <- xx
 #
   if(method=="loclik"){
 #
 #   Use local-likelihood density estimation 
 #   and AIC to choose the bandwidth
 #   based on a grid search
 #
   locfit <- NULL
   if(requireNamespace(locfit, quietly = TRUE)){
 #
     yl <- locfit::locfit(~ x, weights=ywgt, xlim=c(0,1),
           alpha=c(2*smooth,0.3), flim=c(0,1))
 #    gpdf <- locfit(~ x, weights=ywgt, xlim=c(0,1), renorm=TRUE, 
 #          acri="cp", alpha=smooth, flim=c(0,1),ev="grid",mg=binn)
     gpdf <- predict(yl, newdata=r)
     scalef <- binn/sum(gpdf)
     gpdf <- gpdf * scalef
 #
#    if(ci){
#      warning('Confidence intervals require method="gam"\n')
#    }
   }else{
    method="bgk"
   }}
   if(method=="gam"){
#
#     Use local-likelihood density estimation 
#     and AICC to choose the bandwidth
#     based on a grid search
#
#     requireNamespace(mgcv, quietly = TRUE)
      is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
        abs(x - round(x)) < tol
      }
      family <- ifelse(is.wholenumber(m*xx),"poisson","quasi")
      maxsmooth <- max(aicc)
      if(smooth<=0){
        aicc <- cbind(aicc,aicc)
        dimnames(aicc)[[2]] <- c("smooth","aicc")
        for(i in seq(along=aicc[,1])){
         yl <- mgcv::gam(y ~ s(x, bs="cr"),sp=aicc[i,1],
	           family = family, data=data.frame(x=r,y=m*xx))
         df <- summary(yl)$edf + 2
         if(df >= binn - 2){
           aicc[i,2] <- Inf
         }else{
           aicc[i,2] <- log(yl$deviance) + (1 + df/binn)/(1 - df/binn - 2/binn)
         }
        }
        aicc <- aicc[order(aicc[,2]),]
        smooth <- aicc[1,1]
#       print(aicc)
      }
      if(really.missing(smooth, missargs)){
         yl <- mgcv::gam(y ~ s(x, bs="cr"), scale=-1,
	           family = family, data=data.frame(x=r,y=m*xx))
         smooth <- yl$sp
      }else{
       yl <- mgcv::gam(y ~ s(x, bs="cr"), sp=smooth,
            family = family, data=data.frame(x=r,y=m*xx)) 
      }
      if(smooth > maxsmooth){
        cat(paste("Smoothing the maximum amount\n"))
      }else{
        cat(paste("Smoothing using",format(smooth),"\n"))
      }
      gpdf <- mgcv::predict.gam(yl, type = "response")
#
      scalef <- binn/sum(gpdf)
      gpdf <- gpdf * scalef
      if(ci){
        gpdfse <- mgcv::predict.gam(yl, type = "response", se.fit=TRUE)$se.fit
        gpdfse <- as.vector(gpdfse) * scalef
      }
    }
#
   if(method=="bgk"){
#
#     Use Botev. Z.I., Grotowski J.F and Kroese D. P. (2010)
#
      if(really.missing(smooth, missargs)){
        a=bgk_kde(x,n=2^(ceiling(log(binn)/log(2))),MIN=0,MAX=1)
      }else{
        a=bgk_kde(x,n=2^(ceiling(log(binn)/log(2))),MIN=0,MAX=1,smooth=4*smooth/0.35)
      }
#     gpdf <- approx(x=a[1,],y=a[2,],xout=r,rule=2)$y
#     Use an interpolating cubic spline
      gpdf <- spline(x=a[1,],y=a[2,],xout=r)$y
#     gpdf <- predict(smooth.spline(x=a[1,],y=a[2,],df=45*smooth),x=r)$y
#
      scalef <- binn/sum(gpdf)
      gpdf <- gpdf * scalef
    }
#
   if(method=="quick"){
 #
 #     Anscombe transformation to stabilize variances
 #     not quite as good
 #
 #     requireNamespace(modreg, quietly = TRUE)
 #
       vstxx <- 2*sqrt(m*xx + 3/8)
       yl <- loess(vstxx ~ r, span = smooth, degree = 1)
       gpdf <- yl$fitted * yl$fitted/4 - 3/8
       gpdf[gpdf<0] <- 0
       scalef <- binn/sum(gpdf)
       gpdf <- gpdf * scalef
   }
 #
 # calc the CDF
 #
   if(really.missing(ywgt,missargs)){
     cdfg <- seq(along=x)/m
   }else{
     cdfg <- cumsum(ywgt)
   }
   cdfgr <- cumsum(gpdf)/sum(gpdf)
 #
 # calc deciles
 #
   if(really.missing(ywgt,missargs)){
 #   xx <- table(cut(x, breaks = deciles, include.lowest=TRUE,labels=FALSE))/m
     xx <- hist(x,breaks=deciles,plot=FALSE,include.lowest=TRUE)$counts/m
   }else{
     ywgt <- ywgt/sum(ywgt)
     xxx <- cut(x, breaks = deciles, include.lowest=TRUE,labels=FALSE)
     xx <- rep(0,length=length(deciles))
     for (i in unique(xxx)){    
       xx[i] <- sum(ywgt[xxx==i],na.rm=TRUE)
     }
   }
   cdfgdd <- as.vector(xx)[-length(deciles)]
 #
 # next the old way
 #
 #  cdfgdd <- approx(x=c(0,seq(along=x)/m,1),y=c(0,cdfg,1),
 #                xout=deciles[-c(1,length(deciles))])$y
 #  cdfgdd[is.na(cdfgdd)] <- 1
 #  cdfgdd <- diff(c(0,cdfgdd,1))
 #
 # Calculate the entropy of the relative distribution
 #
   egv <- gpdf[gpdf > 0]
   entropy <- sum(egv*log(egv))/binn
 #
   if(method == "histogram"){
     bar <- "yes"
     type <- "n"
   }
   if(really.missing(ylim,missargs)){
     if(bar!="no"){
       add <- TRUE
 #     barplot(height=rbind((length(deciles)-1)*cdfgdd,0.5),
 #        beside=FALSE,
       barplot(height=(length(deciles)-1)*cdfgdd,
          beside=FALSE,
          space=0, width=diff(deciles), col="grey",
          yaxs=yaxs, xlab = xlab, ylab = ylab,
          cex.names=cex, 
          xaxt="n", yaxt="n", ...)
     }
     if(add){
       if(method != "bgk" & method != "gam" & method != "histogram" & method != "quick"){
         graph <- FALSE}
       if(cdfplot){
         lines(x = r, y = cdfgr, lty=lty)
       }else{if(method != "histogram" & bar != "only"){
         lines(x = r, y = gpdf,lty=lty)
        }
       }
     }else{
       if(graph & !cdfplot & bar!="only"){
         plot(x = r, y = gpdf, type = type, lty=lty,
              cex=cex, yaxs=yaxs,xlab = xlab, ylab = ylab,
              xaxt="n", yaxt="n", ...)
       }
       if(graph & cdfplot){
         plot(x = r, y = cdfgr, type = "l", lty=lty,
              ylim=c(0,1), xlim=c(0,1),cex=cex,
              yaxs=yaxs,xlab = xlab, ylab = ylab, 
              xaxt="n", yaxt="n", ...)
         abline(h = (0:10)/10, lty = 1, lwd = 0.5, col=0.2)
         abline(v = (0:10)/10, lty = 1, lwd = 0.5, col=0.2)
       }
     }
   }else{
     if(bar!="no"){
       add <- TRUE
       barplot(height=(length(deciles)-1)*cdfgdd,space=0,
         width=diff(deciles), ylim=ylim, col="grey", xpd=FALSE,
         cex.names=cex, 
         yaxs=yaxs, xlab = xlab, ylab = ylab, xaxt="n", yaxt="n",
         ...)
     }
     if(add){
       if(method != "bgk" & method != "gam" & method != "histogram" & method != "quick"){
         graph <- FALSE}
       if(cdfplot){
         lines(x = r, y = cdfgr, lty=lty)
       }else{if(method != "histogram" & bar!="only"){
         lines(x = r, y = gpdf,lty=lty)
        }
       }
     }else{
       if(graph & !cdfplot & bar!="only"){
         plot(x = r, y = gpdf, type = type, ylim=ylim,lty=lty,
              cex=cex, yaxs=yaxs,xlab = xlab, ylab = ylab, 
              yaxt="n", xaxt="n", ...)
       }
       if(graph & cdfplot){
         plot(x = r, y = cdfgr, type = "l", ylim=ylim, lty=lty,
              xlim=c(0,1), cex=cex,
              yaxs=yaxs,xlab = xlab, ylab = ylab,
              xaxt="n", yaxt="n", ...)
              abline(h = (0:10)/10, lty = 2, lwd = 0.5)
              abline(v = (0:10)/10, lty = 2, lwd = 0.5)
       }
     }
   }
   if(graph & !cdfplot){abline(h = 1, lty = 2)}
   if(graph & cdfplot){abline(a = 0, b = 1, lty = 2)}
 #
   if(graph){
     if(graph & really.missing(ylabs,missargs)){
       if(cdfplot){
        axis(side=2,cex=cex,at=(0:10)/10,mgp=c(1,ylabsloc,0))
       }else{
        axis(side=2,cex=cex,mgp=c(3,ylabsloc,0))}}
 #     if(cdfplot){axis(2,cex=cex,at=(0:10)/10,mgp=c(2,1,0))}else{
 #     axis(2,cex=cex,mgp=c(2,0.5,0))}}
     else{
       if(cdfplot){axis(side=2,cex=cex,at=(0:10)/10,mgp=c(1,ylabsloc,0))}else{
        axis(side=2,cex=cex,mgp=c(3,ylabsloc,0))}
 #     if(cdfplot){axis(2,cex=cex,at=(0:10)/10,mgp=c(2,1,0))}else{
 #     axis(2,cex=cex,mgp=c(2,0.5,0))}
       if(length(ylabs) == 1 & ylabs[1] == TRUE){ylabs<-pretty(y)}
        yxlabs <- rmdata(ylabs,y)$x
   #
   #   next places it on the top inside (outside is tck0.02)
   #   mgp c(3,2,0) is normal outside c(3,-2,0) is normal inside
   #
 #      axis(side = 4, at = yxlabs, labels = paste(ylabs), tck=-0.02, 
 #           cex=cex, mgp=c(2,0.5,0))
 #          cex=cex, mgp=c(1,.6,0))
       axis(side = 4, at = yxlabs, labels = paste(ylabs), tck=-0.02,
            cex=cex, mgp=c(1,ylabsloc,0))
     }
   #
     if(really.missing(yolabs,missargs)){
       if(cdfplot){axis(1,cex=cex,at=(0:10)/10)}else{
       axis(1,cex=cex,mgp=c(3,yolabsloc,0))}}
 #     if(cdfplot){axis(1,cex=cex,at=(0:10)/10,mgp=c(2,0.5,0))}else{
 #     axis(1,cex=cex,mgp=c(2,0.5,0))}}
     else{
       if(cdfplot){axis(1,cex=cex,at=(0:10)/10)}else{
       axis(1,cex=cex,mgp=c(3,yolabsloc,0))}
 #      if(cdfplot){axis(1,cex=cex,at=(0:10)/10,mgp=c(2,0.5,0))}else{
 #     axis(1,cex=cex,mgp=c(2,0.5,0))}
       if(length(yolabs) == 1 & yolabs[1] == TRUE){yolabs<-pretty(yo)}
       if(really.missing(yo,missargs)){
        yoxlabs <- yolabs
       }else{
        yoxlabs <- rmdata(yolabs,yo)$x
       }
   #
   #   next places it on the top outside
   #
       if(really.missing(yolabslabs,missargs)){yolabslabs <- paste(yolabs)}
        axis(side = 3, at = yoxlabs, labels = yolabslabs, tck=-0.02,
            cex=cex, mgp=c(3,yolabsloc,0))
 #       axis(side = 3, at = yoxlabs, labels = yolabslabs, tck=-0.02,
 #          cex=cex, mgp=c(2,0.5,0))
     }
   }
   gpdf.ci <- NULL
   if(ci){
     if(n > 0){
      if(method=="gam"){
       varg <- gpdfse^2 + (1/(2*sqrt(pi)))*gpdf*gpdf/(n*smooth)
      }else{
       varg <- (1/(2*sqrt(pi)))*(1/n+gpdf/m)*gpdf
      }
     }else{
      if(method=="gam"){
       varg <- gpdfse^2
      }else{
       varg <- (1/(2*sqrt(pi)))*(gpdf/m)*gpdf
      }
     }
     gpdf.ci <- list(x=r,
               l=gpdf-1.96*sqrt(varg),
               u=gpdf+1.96*sqrt(varg))
     if(graph & !cdfplot & bar!="only"){
       points(x=gpdf.ci$x,y=gpdf.ci$l,pch=".",cex=0.3)
       points(x=gpdf.ci$x,y=gpdf.ci$u,pch=".",cex=0.3)
     }
   }
   cdf.out <- list(x=seq(along=x)/m,y=cdfg)
 #
 # Calculate polarizations
 #
   if(!really.missing(yo,missargs)){
    if(!(show=="residual" & decomp=="locmult") & rpmult & !really.missing(yo,missargs)){
     x <- rmdata(y=y, yo=yo, ywgt=ywgt,yowgt=yowgt,
                 show="residual", decomp="locmult",
                 location=location, scale=scale, missargs=missargs)$x
    }else{
     if(decomp=="covariate" | (show!="residual" & !rpmult)){
      x <- rmdata(y=y, yo=yo, ywgt=ywgt,yowgt=yowgt,
                  show="residual", decomp="locadd",
                  location=location, scale=scale, missargs=missargs)$x
     }
    }
   }
   w <- abs(x - 0.5)
   selu <- x  > 0.5
   sell <- x <= 0.5
   if(!really.missing(ywgt,missargs)){
    rpm <- 4 * sum(w*ywgt)/sum(ywgt) - 1
    c1  <- wtd.var(w,weight=ywgt)
    rpu <- 4 * sum(w*ywgt*selu)/sum(ywgt*selu) - 1
    c1u <- wtd.var(w,weight=ywgt*selu)
    rpl <- 4 * sum(w*ywgt*sell)/sum(ywgt*sell) - 1
    c1l <- wtd.var(w,weight=ywgt*sell)
   }else{
    rpm <- 4 * mean(w) - 1
    c1  <- var(w)
    rpu <- 4 * sum(w*selu)/sum(selu) - 1
    c1u <- var(w[selu])
    rpl <- 4 * sum(w*sell)/sum(sell) - 1
    c1l <- var(w[sell])
   }
   if(really.missing(yo,missargs)){
    serp  <- 4*sqrt(c1/m)
    serpu <- 8*sqrt(c1u/m)
    serpl <- 8*sqrt(c1l/m)
   }else{
    if(show=="residual" & decomp=="locmult" | rpmult){
     x <- rmdata(y=yo, yo=y, ywgt=yowgt,yowgt=ywgt,
                 show="residual", decomp="locmult",
                 location=location, scale=scale, missargs=missargs)$x
    }else{
     x <- rmdata(y=yo, yo=y, ywgt=yowgt,yowgt=ywgt,
                 show="residual", decomp="locadd",
                 location=location, scale=scale, missargs=missargs)$x
    }
    w <- abs(x - 0.5)
    selu <- x  > 0.5
    sell <- x <= 0.5
    if(really.missing(yowgt,missargs)){
     c2  <- var(w)
     c2u <- var(w[selu])
     c2l <- var(w[sell])
    }else{
     c2  <- wtd.var(w,weight=yowgt)
     c2u <- wtd.var(w,weight=yowgt*selu)
     c2l <- wtd.var(w,weight=yowgt*sell)
    }
    serp  <- 4*sqrt(c1/m + c2/n)
    serpu <- 8*sqrt(c1u/m + c2u/n)
    serpl <- 8*sqrt(c1l/m + c2l/n)
   }
   rpm <- rpm+1.96*serp *c(-1, 0, 1)
   rpu <- rpu+1.96*serpu*c(-1, 0, 1)
   rpl <- rpl+1.96*serpl*c(-1, 0, 1)
 #
 # Output
 #
   if(!quiet){
       list(x=r,y=gpdf,ci=gpdf.ci,
        rp = rpm, rpl= rpl, rpu= rpu,
        cdf=cdf.out, deciles=cdfgdd, entropy=entropy,
        smooth=smooth,aicc=aicc
       )
   }else{
      invisible(
       list(x=r,y=gpdf,ci=gpdf.ci,
        rp = rpm, rpl= rpl, rpu= rpu,
        cdf=cdf.out, deciles=cdfgdd, entropy=entropy,
        smooth=smooth,aicc=aicc
       )
      )
   }
 }         
 #
 # median polarization index
 #
 "rpy" <- function(y=FALSE, yo=FALSE, ywgt=FALSE,yowgt=FALSE, pvalue=FALSE,
                    z=FALSE, zo=FALSE,
                    show="residual", decomp="locadd",
                    location="median", scale="IQR",
                    rpmult=FALSE
                 ){
 #
 #	pvalue produce two-sided p-value for the polarization index
 #
 #	See "reldist" for a description of the other variables
 #
 #	# OUTPUT list components
 #	  	95% confidence interval for the median relative polarization
 #               as lower bound, estimate, upper bound.
 #
 # -------------------------------------------------------------------------
 #
   if(missing(yo)){
     rdata <- list(x=y, m=length(x), ywgt=rep(1,length(x)))
   }else{
     if(decomp=="covariate" & !missing(z) & !missing(zo)){
      rdata <- rmdata(y=y, yo=yo, ywgt=ywgt,yowgt=yowgt,
                  show=show, decomp=decomp,
                  z=z, zo=zo,
  	         location=location, scale=scale) 
      y <- rdata$y
      yo <- rdata$yo
      ywgt <- rdata$ywgt
      yowgt <- rdata$yowgt
      if(!rpmult){
        decomp <- "locadd"
      }else{
        decomp <- "locmult"
      }
     }
     rdata <- rmdata(y=y, yo=yo, ywgt=ywgt,yowgt=yowgt,
                  show=show, decomp=decomp,
  	         location=location, scale=scale) 
     y <- rdata$y
     yo <- rdata$yo
     yowgt <- rdata$yowgt
     n <- rdata$n
   }
   x <- rdata$x
   m <- rdata$m
   ywgt <- rdata$ywgt
 #
 # Calculate polarizations
 #
   w <- abs(x - 0.5)
   rpm <- 4 * sum(w*ywgt)/sum(ywgt) - 1
   c1  <- wtd.var(w,weight=ywgt)
   x <- rmdata(y=yo, yo=y, ywgt=yowgt,yowgt=ywgt,
               show=show, decomp=decomp,
               location=location, scale=scale)$x
   w <- abs(x - 0.5)
   c2  <- wtd.var(w,weight=yowgt)
   serp  <- 4*sqrt(c1/m + c2/n)
   if(pvalue){
     c(rpm+1.96*serp *c(-1, 0, 1), 1-pnorm(abs(rpm/serp)))
   }else{
     rpm+1.96*serp *c(-1, 0, 1)
   }
 }
 #
 # lower and upper polarization indices
 #
 "rpluy" <- function(y=FALSE, yo=FALSE, ywgt=FALSE,yowgt=FALSE, pvalue=FALSE,
                     upper=FALSE,lower=TRUE,
                     show="residual", decomp="locadd",
                     location="median", scale="IQR",
                     rpmult=FALSE,
                     z=FALSE, zo=FALSE
                    ){
 #
 #	pvalue produce two-sided p-value for the polarization index
 #	upper  create upper polarization index
 #	lower  create lower polarization index
 #
 #	See "reldist" for a description of the other variables
 #
 #	# OUTPUT list components
 #	  	95% confidence interval for the median relative polarization
 #                as lower bound, estimate, upper bound.
 #
 # -------------------------------------------------------------------------
 #
   if(upper){lower<-FALSE}
 #
   if(missing(yo)){
     rdata <- list(x=y, m=length(x), ywgt=rep(1,length(x)))
   }else{
     if(decomp=="covariate" & !missing(z) & !missing(zo)){
      rdata <- rmdata(y=y, yo=yo, ywgt=ywgt,yowgt=yowgt,
                  show=show, decomp=decomp,
                  z=z, zo=zo,
  	         location=location, scale=scale) 
      y <- rdata$y
      yo <- rdata$yo
      ywgt <- rdata$ywgt
      yowgt <- rdata$yowgt
      if(!rpmult){
        decomp <- "locadd"
      }else{
        decomp <- "locmult"
      }
     }
     rdata <- rmdata(y=y, yo=yo, ywgt=ywgt,yowgt=yowgt,
                     show=show, decomp=decomp,
                     z=z, zo=zo,
  	            location=location, scale=scale) 
     y <- rdata$y
     yo <- rdata$yo
     yowgt <- rdata$yowgt
     n <- rdata$n
   }
   x <- rdata$x
   m <- rdata$m
   ywgt <- rdata$ywgt
 #
   w <- abs(x - 0.5)
   if(upper){
     select <- x  > 0.5
   }else{
     select <- x <= 0.5
   }
   rpu <- 4 * sum(w*ywgt*select)/sum(ywgt*select) - 1
   c1u <- wtd.var(w,weight=ywgt*select)
   x <- rmdata(y=yo, yo=y, ywgt=yowgt,yowgt=ywgt,
               show=show, decomp=decomp,
               location=location, scale=scale)$x
   w <- abs(x - 0.5)
   if(upper){
    select <- x  > 0.5
   }else{
    select <- x <= 0.5
   }
   c2u <- wtd.var(w,weight=yowgt*select)
   serp <- 8*sqrt(c1u/m + c2u/n)
   if(pvalue){
     c(rpu+1.96*serp*c(-1, 0, 1), 1-pnorm(abs(rpu/serp)))
   }else{
     rpu+1.96*serp*c(-1, 0, 1)
   }
 }
 "wtd.var" <- function(x, weight=rep(1,length(x))) {
 	sum(x * x * weight)/sum(weight) - (sum(x * weight)/sum(weight))^2
 }
 "wtd.mean" <- function(x, weight=rep(1,length(x))) {
 	(sum(x * weight)/sum(weight))
 }
 "rmdata" <- function(y, yo, ywgt=FALSE,yowgt=FALSE,
                      z=FALSE, zo=FALSE,
                      show="none", decomp="locadd",
                      location="median", scale="IQR",
                      robust=FALSE, missargs=FALSE
                     ){
 #
 #  First find missing values
 #
  if(missing(missargs)){
   missargs <- NULL
   if(missing(y)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yo)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(ywgt)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yowgt)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(show)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(decomp)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(location)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(scale)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(z)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(zo)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(robust)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   names(missargs) <- c(
   "y",
   "yo",
   "ywgt",
   "yowgt",
   "show",
   "decomp",
   "location",
   "scale",
   "z",
   "zo",
   "robust")
  }else{
   if(missing(y)){missargs["y"] <- TRUE}
   if(missing(yo)){missargs["yo"] <- TRUE}
   if(missing(ywgt)){missargs["ywgt"] <- TRUE}
   if(missing(yowgt)){missargs["yowgt"] <- TRUE}
   if(missing(show)){missargs["show"] <- TRUE}
   if(missing(decomp)){missargs["decomp"] <- TRUE}
   if(missing(location)){missargs["location"] <- TRUE}
   if(missing(scale)){missargs["scale"] <- TRUE}
   if(missing(z)){missargs["z"] <- TRUE}
   if(missing(zo)){missargs["zo"] <- TRUE}
  }
 #
 #  First remove NAs
 #
  if(!really.missing(ywgt,missargs)){
    nas <- is.na(y) | is.na(ywgt)
    ywgt <- ywgt[!nas]/sum(ywgt[!nas])
  }else{
    nas <- is.na(y)
    m <- length(y)
    ywgt<-rep(1,length=m)/m
    missargs["ywgt"] <- FALSE
  }
  if(all(nas)){ stop("All the y values are missing!") }
  y <- y[!nas]
  if(!really.missing(z,missargs)){z <- z[!nas]}
  if(!really.missing(yowgt,missargs)){
    nas <- is.na(yo) | is.na(yowgt)
    yowgt <- yowgt[!nas]/sum(yowgt[!nas])
  }else{
    nas <- is.na(yo)
    n <- length(yo)
    yowgt<-rep(1,length=n)/n
    missargs["yowgt"] <- FALSE
  }
  if(all(nas)){ stop("All the yo values are missing!") }
  yo <- yo[!nas]
  if(!really.missing(zo,missargs)){zo <- zo[!nas]}
#
  n <- length(yo)
  m <- length(y)
#
# Do matching
#
# First covariate
#
  if(decomp=="covariate" & !really.missing(z,missargs)){
   if(show=="effect"){
    ywgt <- yowgt*rdsamp(z,zo,ywgt,yowgt)
    ywgt <- ywgt/sum(ywgt)
    y <- yo
    m <- n
   }else{
    yowgt <- yowgt*rdsamp(z,zo,ywgt,yowgt)
    yowgt <- yowgt/sum(yowgt)
   }
  }
#
# Now multiplicative
#
  if(decomp=="locmult"){
   if(show=="residual"){
    if(location=="median"){
     if(!really.missing(yowgt,missargs) ){
      yo <- wtd.quantile(y,weight=ywgt)*yo/wtd.quantile(yo,weight=yowgt)
     }else{
      yo <- median(y)*yo/median(yo)
     }
    }else{
     if(!really.missing(yowgt,missargs) ){
      yo <- wtd.mean(y,weight=ywgt)*yo/wtd.mean(yo,weight=yowgt)
     }else{
      yo <- mean(y)*yo/mean(yo)
     }
    }
   }
   if(show=="effect"){
    if(location=="median"){
     if(!really.missing(ywgt,missargs)){
      y <- wtd.quantile(y,weight=ywgt)*yo/wtd.quantile(yo,weight=yowgt)
      ywgt <- yowgt
      m <- n
     }else{
      y <- median(y)*yo/median(yo)
      m <- n
     }
    }else{
     if(!really.missing(ywgt,missargs)){
      y <- wtd.mean(y,weight=ywgt)*yo/wtd.mean(yo,weight=yowgt)
      ywgt <- yowgt
      m <- n
     }else{
      y <- mean(y)*yo/mean(yo)
      m <- n
     }
    }
   }
  }
#
# Now additive
#
  if(decomp=="locadd"){
   if(show=="residual"){
    if(location=="median"){
     if(!really.missing(yowgt,missargs) ){
      yo <- wtd.quantile(y,weight=ywgt)+yo-wtd.quantile(yo,weight=yowgt)
     }else{
      yo <- median(y)+(yo-median(yo))
     }
    }else{
     if(!really.missing(yowgt,missargs) ){
      yo <- wtd.mean(y,weight=ywgt)+yo-wtd.mean(yo,weight=yowgt)
     }else{
      yo <- mean(y)+(yo-mean(yo))
     }
    }
   }
   if(show=="effect"){
    if(location=="median"){
     if(!really.missing(ywgt,missargs)){
      y <- wtd.quantile(y,weight=ywgt)+yo-wtd.quantile(yo,weight=yowgt)
      ywgt <- yowgt
      m <- n
     }else{
      y <- median(y)+(yo-median(yo))
      m <- n
     }
    }else{
     if(!really.missing(ywgt,missargs)){
      y <- wtd.mean(y,weight=ywgt)+yo-wtd.mean(yo,weight=yowgt)
      ywgt <- yowgt
      m <- n
     }else{
      y <- mean(y)+(yo-mean(yo))
      m <- n
     }
    }
   }
  }
#
# Now location and scale
#
  if(decomp=="lsadd"){
   if(show=="residual"){
    if(location=="median" & scale=="IQR"){
     if(!really.missing(yowgt,missargs)){
      yo <- wtd.quantile(y,weight=ywgt)+wtd.iqr(y,weight=ywgt)*
            (yo-wtd.quantile(yo,weight=yowgt))/wtd.iqr(yo,weight=yowgt)
     }else{
      yo <- median(y)+iqr(y)*(yo-median(yo))/iqr(yo)
     }
    }else{
     if(!really.missing(yowgt,missargs)){
      yo <- wtd.mean(y,weight=ywgt)+sqrt(wtd.var(y,weight=ywgt))*
            (yo-wtd.mean(yo,weight=yowgt))/sqrt(wtd.var(yo,weight=yowgt))
     }else{
      yo <- mean(y)+sqrt(var(y))*(yo-mean(yo))/sqrt(var(yo))
     }
    }
   }
   if(show=="effect"){
    if(location=="median"){
     if(!really.missing(ywgt,missargs)){
      y <- wtd.quantile(y,weight=ywgt)+wtd.iqr(y,weight=ywgt)*
            (yo-wtd.quantile(yo,weight=yowgt))/wtd.iqr(yo,weight=yowgt)
      ywgt <- yowgt
      m <- n
     }else{
      y <- median(y)+iqr(y)*(yo-median(yo))/iqr(yo)
      m <- n
     }
    }else{
     if(!really.missing(ywgt,missargs)){
      y <- wtd.mean(y,weight=ywgt)+sqrt(wtd.var(y,weight=ywgt))*
            (yo-wtd.mean(yo,weight=yowgt))/sqrt(wtd.var(yo,weight=yowgt))
      ywgt <- yowgt
      m <- n
     }else{
      y <- mean(y)+sqrt(var(y))*(yo-mean(yo))/sqrt(var(yo))
      m <- n
     }
    }
   }
  }
#
#  sy is 1:n
#  ys is the ordered y's
#  ry is the ranks of ys in the joint vector of {yo, y}
#  ry - sy is the number of yo's le ys
#
  sy <- seq(along = y)
  ys <- sort.list(y)
  if(!really.missing(ywgt,missargs)){ywgt <- ywgt[ys]}
  ys <- y[ys]
# ry <- sort.list(sort.list(c(yo, ys)))[sy + n]
  ry <- sort.list(sort.list(c(yo, ys)))
  if(robust){
#	y <- y*(1+sqrt(var(y))*(runif(m)-0.5)/100000)}
        jv <- c(yo, ys)
        njv <- !is.na(jv)
        for(i in unique(jv)[duplicated(jv)]) {
          which <- jv == i & njv
          ry[which] <- mean(ry[which])
        }
  } 
  ry <- ry[sy + n]
  if(really.missing(yowgt,missargs)){
    x <- (ry - sy)/n
  }else{
    yos <- sort.list(yo)
    yowgts <- yowgt[yos]
#
#   x is the new relative data, but note it has weight ywgt
#
    x <- c(0,cumsum(yowgts))[ry - sy + 1]
  }
  x[x > 1] <- 1
  x[x < 0] <- 0
# if(really.missing(ywgt,missargs) ){ywgt <- rep(1,length=m)/m}
# if(really.missing(yowgt,missargs)){yowgt <- rep(1,length=n)/n}
#
#  Use x[order(order(y))]
#  to get the relative data back in original order
#
#  returned values are in the original order
#
  rorder <- order(order(y))
  invisible(list(x=x[rorder],wgt=ywgt[rorder],y=y,yo=yo,ywgt=ywgt[rorder],yowgt=yowgt,n=n,m=m))
}
"rddata" <- function(y, yo, ywgt=FALSE,yowgt=FALSE,
                     smooth=1,
                     show="none", decomp="locadd",
                     location="median", scale="IQR",
                     z=FALSE, zo=FALSE,
                     robust=FALSE, missargs=FALSE
                    ){
#
  if(missing(missargs)){
   missargs <- NULL
   if(missing(y)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yo)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(ywgt) | is.logical(ywgt)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yowgt) | is.logical(yowgt)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(show)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(decomp)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(location)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(scale)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(z)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(zo)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(smooth)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
#  if(missing(pool)){
#    missargs <- c(missargs,TRUE)
#  }else{
#    missargs <- c(missargs,FALSE)
#  }
   names(missargs) <- c(
   "y",
   "yo",
   "ywgt",
   "yowgt",
   "show",
   "decomp",
   "location",
   "scale",
   "z",
   "zo",
   "smooth")
#  "pool")
  }else{
   if(missing(y)){missargs["y"] <- TRUE}
   if(missing(yo)){missargs["yo"] <- TRUE}
   if(missing(ywgt)){missargs["ywgt"] <- TRUE}
   if(missing(yowgt)){missargs["yowgt"] <- TRUE}
   if(missing(show)){missargs["show"] <- TRUE}
   if(missing(decomp)){missargs["decomp"] <- TRUE}
   if(missing(location)){missargs["location"] <- TRUE}
   if(missing(scale)){missargs["scale"] <- TRUE}
   if(missing(z)){missargs["z"] <- TRUE}
   if(missing(smooth)){missargs["zo"] <- TRUE}
   mnames <- names(missargs)
#  if(missing(pool)){
#   missargs <- c(missargs,TRUE)
#  }else{
#   missargs <- c(missargs,FALSE)
#  }
#  names(missargs) <- c(mnames,"pool")
  }
# if(    really.missing(pool,missargs) 
#     & !really.missing(smooth,missargs)){ pool <- smooth }
#
#  First remove NAs
#
  if(!really.missing(ywgt,missargs)){
    nas <- is.na(y) | is.na(ywgt)
    ywgt <- ywgt[!nas]/sum(ywgt[!nas])
  }else{
    nas <- is.na(y)
  }
  if(all(nas)){ stop("All the y values are missing!") }
  y <- y[!nas]
  if(!really.missing(yowgt,missargs)){
    nas <- is.na(yo) | is.na(yowgt)
    yowgt <- yowgt[!nas]/sum(yowgt[!nas])
  }else{
    nas <- is.na(yo)
  }
  if(all(nas)){ stop("All the yo values are missing!") }
  yo <- yo[!nas]
#
  n <- length(yo)
  m <- length(y)
#
# Do matching
#
# First covariate
#
  if(decomp=="covariate" & !really.missing(z,missargs) &
!really.missing(zo,missargs)){
   if(really.missing(yowgt,missargs)){yowgt<-rep(1,length=n)/n}
   if(really.missing( ywgt,missargs)){ ywgt<-rep(1,length=m)/m}
   if(show=="effect"){
    ywgt <- yowgt*rdsamp(z,zo,ywgt,yowgt)
    ywgt <- ywgt/sum(ywgt)
    y <- yo
    m <- n
   }else{
    yowgt <- yowgt*rdsamp(z,zo,ywgt,yowgt)
    yowgt <- yowgt/sum(yowgt)
   }
  }
#
# Now multiplicative
#
  if(decomp=="locmult"){
   if(show=="residual"){
    if(location=="median"){
     if(!really.missing(yowgt,missargs)){
      yo <- wtd.quantile(y,weight=ywgt)*yo/wtd.quantile(yo,weight=yowgt)
     }else{
      yo <- median(y)*yo/median(yo)
     }
    }else{
     if(!really.missing(yowgt,missargs)){
      yo <- wtd.mean(y,weight=ywgt)*yo/wtd.mean(yo,weight=yowgt)
     }else{
      yo <- mean(y)*yo/mean(yo)
     }
    }
   }
   if(show=="effect"){
    if(location=="median"){
     if(!really.missing(ywgt,missargs)){
      y <- wtd.quantile(y,weight=ywgt)*yo/wtd.quantile(yo,weight=yowgt)
      ywgt <- yowgt
      m <- n
     }else{
      y <- median(y)*yo/median(yo)
      m <- n
     }
    }else{
     if(!really.missing(ywgt,missargs)){
      y <- wtd.mean(y,weight=ywgt)*yo/wtd.mean(yo,weight=yowgt)
      ywgt <- yowgt
      m <- n
     }else{
      y <- mean(y)*yo/mean(yo)
      m <- n
     }
    }
   }
  }
#
# Now additive
#
  if(decomp=="locadd"){
   if(show=="residual"){
    if(location=="median"){
     if(!really.missing(yowgt,missargs)){
      yo <- wtd.quantile(y,weight=ywgt)+yo-wtd.quantile(yo,weight=yowgt)
     }else{
      yo <- median(y)+(yo-median(yo))
     }
    }else{
     if(!really.missing(yowgt,missargs)){
      yo <- wtd.mean(y,weight=ywgt)+yo-wtd.mean(yo,weight=yowgt)
     }else{
      yo <- mean(y)+(yo-mean(yo))
     }
    }
   }
   if(show=="effect"){
    if(location=="median"){
     if(!really.missing(ywgt,missargs)){
      y <- wtd.quantile(y,weight=ywgt)+yo-wtd.quantile(yo,weight=yowgt)
      ywgt <- yowgt
      m <- n
     }else{
      y <- median(y)+(yo-median(yo))
      m <- n
     }
    }else{
     if(!really.missing(ywgt,missargs)){
      y <- wtd.mean(y,weight=ywgt)+yo-wtd.mean(yo,weight=yowgt)
      ywgt <- yowgt
      m <- n
     }else{
      y <- mean(y)+(yo-mean(yo))
      m <- n
     }
    }
   }
  }
#
# Now location and scale
#
  if(decomp=="lsadd"){
   if(show=="residual"){
    if(location=="median" & scale=="IQR"){
     if(!really.missing(yowgt,missargs)){
      yo <- wtd.quantile(y,weight=ywgt)+wtd.iqr(y,weight=ywgt)*
            (yo-wtd.quantile(yo,weight=yowgt))/wtd.iqr(yo,weight=yowgt)
     }else{
      yo <- median(y)+iqr(y)*(yo-median(yo))/iqr(yo)
     }
    }else{
     if(!really.missing(yowgt,missargs)){
      yo <- wtd.mean(y,weight=ywgt)+sqrt(wtd.var(y,weight=ywgt))*
            (yo-wtd.mean(yo,weight=yowgt))/sqrt(wtd.var(yo,weight=yowgt))
     }else{
      yo <- mean(y)+sqrt(var(y))*(yo-mean(yo))/sqrt(var(yo))
     }
    }
   }
   if(show=="effect"){
    if(location=="median"){
     if(!really.missing(ywgt,missargs)){
      y <- wtd.quantile(y,weight=ywgt)+wtd.iqr(y,weight=ywgt)*
            (yo-wtd.quantile(yo,weight=yowgt))/wtd.iqr(yo,weight=yowgt)
      ywgt <- yowgt
      m <- n
     }else{
      y <- median(y)+iqr(y)*(yo-median(yo))/iqr(yo)
      m <- n
     }
    }else{
     if(!really.missing(ywgt,missargs)){
      y <- wtd.mean(y,weight=ywgt)+sqrt(wtd.var(y,weight=ywgt))*
            (yo-wtd.mean(yo,weight=yowgt))/sqrt(wtd.var(yo,weight=yowgt))
      ywgt <- yowgt
      m <- n
     }else{
      y <- mean(y)+sqrt(var(y))*(yo-mean(yo))/sqrt(var(yo))
      m <- n
     }
    }
   }
  }
#
	sy <- sort(yo)
	sy <- c(min(y,yo),(sy[-1] - 0.5*diff(sy)),max(y,yo))                  
	sy <- unique(sy)     
#
	po <- hist(yo,breaks=sy,plot=FALSE,include.lowest = TRUE)$counts
#
if(really.missing(yowgt,missargs)){
  po <- hist(yo,breaks=sy,plot=FALSE,include.lowest = TRUE)$counts
}else{
  yowgt <- yowgt/mean(yowgt)
  xxx <- cut(yo,breaks=sy,include.lowest = TRUE,labels=FALSE)
  po <- rep(0,length=length(sy)-1)
  for (i in unique(xxx)){    
    po[i] <- sum(yowgt[xxx==i],na.rm=TRUE)
  }
}
#
# remove zero
#
if(sum(po < 1) > 0){
	eee <- seq(along=po)[po < 1] 
	if(eee[length(eee)] == length(po)){eee[length(eee)] <- eee[length(eee)]-1}
#	eee <- eee[eee < length(po)]       
#
	sy <- sy[-(eee+1)]
	if(really.missing(yowgt,missargs)){
	  po <- hist(yo,breaks=sy,plot=FALSE,include.lowest = TRUE)$counts
#         po <- po/n
	}else{
	  yowgt <- yowgt/mean(yowgt)
	  xxx <- cut(yo,breaks=sy,include.lowest = TRUE,labels=FALSE)
	  po <- rep(0,length=length(sy)-1)
	  for (i in unique(xxx)){    
	    po[i] <- sum(yowgt[xxx==i],na.rm=TRUE)
	  }
# 	  po <- po/mean(po)
	}
}
	if(smooth>=5){
#
# remove small cells
#
        ggg <- seq(5,smooth,length=4)
#
#       Pool successively
#
        for ( px in ggg ){

	eee <- seq(along=po)[po <= px] 
#
	if(length(eee)>0 & length(sy) > 2){
	  if(eee[length(eee)] == length(po)){
		eee[length(eee)] <- eee[length(eee)]-1}
	   sy <- sy[-(eee+1)]
	}
	if(really.missing(yowgt,missargs)){
	  po <- hist(yo,breaks=sy,plot=FALSE,include.lowest = TRUE)$counts
	}else{
	  yowgt <- yowgt/mean(yowgt)
	  xxx <- cut(yo,breaks=sy,include.lowest = TRUE,labels=FALSE)
	  po <- rep(0,length=length(sy)-1)
	  for (i in unique(xxx)){    
	    po[i] <- sum(yowgt[xxx==i],na.rm=TRUE)
	  }
#	  po <- po/mean(po)
	}
        }
	}
        po <- po / n
#
	if(really.missing(yowgt,missargs)){
	  pix <- hist(y,breaks=sy,plot=FALSE,include.lowest = TRUE)$counts
	}else{
	  ywgt <- ywgt/mean(ywgt)
	  xxx <- cut(y,breaks=sy,include.lowest = TRUE,labels=FALSE)
	  pix <- rep(0,length=length(sy)-1)
	  for (i in unique(xxx)){    
	    pix[i] <- sum(ywgt[xxx==i],na.rm=TRUE)
	  }
	}
#	pix <- as.vector(pix)+0.5
	pix <- pix/sum(pix)
	xx <- cumsum(po)    
	yy <- pix/po    
	nx <- length(xx)  
#
  xx[xx > 1] <- 1
  xx[xx < 0] <- 0
  if(really.missing(ywgt,missargs)){ywgt <- rep(1,length=m)/m}
#
  invisible(list(x=xx,gx=yy,nx=nx,pix=pix,po=po,y=y,yo=yo,ywgt=ywgt,yowgt=yowgt,
       breaks=sy,n=n,m=m))
}
rddist <- function(y, yo, add=FALSE, ylim=NULL,ywgt=FALSE,yowgt=FALSE,
                   show="none", decomp="locadd",
                   location="median", scale="IQR",
                   rpmult=FALSE, 
		   quiet=TRUE,
                   z=FALSE, zo=FALSE,
                   smooth=1, binn=100,
                   ci=FALSE,graph=TRUE, cex=0.8,yaxs="r",
                   xlab="baseline proportion", ylab="density", 
                   yolabslabs=NULL, ylabslabs=NULL,
                   yolabs=pretty(yo), ylabs=pretty(y),
                   missargs,
		   ...) {
#
  options(warn=-1)
 #
 #  First find missing values
 #
  if(missing(missargs)){
   missargs <- NULL
   if(missing(y)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yo)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(ywgt)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yowgt)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(show)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(decomp)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(location)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(scale)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(z)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(zo)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(ylabs)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yolabs)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(ylabslabs)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(yolabslabs)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   if(missing(ylim)){missargs <- c(missargs,TRUE)}else{missargs <- c(missargs,FALSE)}
   names(missargs) <- c(
   "y",
   "yo",
   "ywgt",
   "yowgt",
   "show",
   "decomp",
   "location",
   "scale",
   "z",
   "zo",
   "ylabs",
   "yolabs",
   "ylabslabs",
   "yolabslabs",
   "ylim")
  }else{
   if(missing(y)){missargs["y"] <- TRUE}
   if(missing(yo)){missargs["yo"] <- TRUE}
   if(missing(ywgt)){missargs["ywgt"] <- TRUE}
   if(missing(yowgt)){missargs["yowgt"] <- TRUE}
   if(missing(show)){missargs["show"] <- TRUE}
   if(missing(decomp)){missargs["decomp"] <- TRUE}
   if(missing(location)){missargs["location"] <- TRUE}
   if(missing(scale)){missargs["scale"] <- TRUE}
   if(missing(z)){missargs["z"] <- TRUE}
   if(missing(zo)){missargs["zo"] <- TRUE}
   if(missing(ylabs)){missargs["ylabs"] <- TRUE}
   if(missing(yolabs)){missargs["yolabs"] <- TRUE}
   if(missing(ylabslabs)){missargs["ylabslabs"] <- TRUE}
   if(missing(yolabslabs)){missargs["yolabslabs"] <- TRUE}
   if(missing(ylim)){missargs["ylim"] <- TRUE}
  }
# if(missing(pool) & !really.missing(smooth,missargs)){ pool <- smooth }
  if(really.missing(yo,missargs)){
   x <- sort(y)
   m <- length(x)
   n <- 0
   wgt <- rep(1/m,length=m)
   gx <- x
   nx <- n
   pix <- wgt
   rmd <- list(x=x,wgt=wgt,n=n,m=m)
  }else{
   rmd <- rddata(y=y, yo=yo, ywgt=ywgt,yowgt=yowgt,
                 show=show, decomp=decomp, smooth=smooth,
                 location=location, scale=scale,
                 missargs)
   gx <- rmd$gx
   nx <- rmd$nx
   pix <- rmd$pix
   po <- rmd$po
   ywgt <- rmd$ywgt
   yowgt <- rmd$yowgt
   missargs[ "ywgt"] <- FALSE
   missargs["yowgt"] <- FALSE
   x <- rmd$x
   n <- rmd$n
   m <- rmd$m
  }
#
  x1 <- rep(c(0,x),rep(2,nx+1))[-1]
  y1 <- rep(gx,c(rep(2,nx-1),3))
  y2 <- rep(gx,c(rep(2,nx-1),4))[-1]             
  x2 <- c(rep(x,rep(2,nx)),1)  
#
# requireNamespace(mgcv, quietly = TRUE)
  r <- seq(0, 1, length = binn + 1)[-1] - 0.5/binn
  pdfgdd <- approx(x=c(0,x),y=c(0,cumsum(pix)), xout=r)
  pdfgdd$y <- diff(c(0,pdfgdd$y))
  pdfgdd$y[is.na(pdfgdd$y)] <- 1
  pdfgdd$y <- binn*mgcv::predict.gam(mgcv::gam(y ~ s(x), data=pdfgdd), type="response")
#
# Calculate the entropy of the relative distribution
#
  dx <- rep(1/binn,binn)
  egv <- pdfgdd$y[pdfgdd$y > 0]
  edx <- dx[pdfgdd$y > 0]
  entropy <- sum(egv*log(egv)*edx)
#
  chisq <- sum((egv-1)^2*edx)
#
  if(graph){
   if(!add){
    if(really.missing(ylim,missargs)){
      plot(x = x1, y = y1, type = "n", xaxt="n", yaxs=yaxs,
        xlab = xlab, ylab = ylab)
    }else{
      plot(x = x1, y = y1, type = "n", ylim=ylim,xaxt="n", yaxs=yaxs,
        xlab = xlab, ylab = ylab)
    }
    lines(x = r, y = pdfgdd$y, lty=1)
    if(really.missing(ylabs,missargs)){
     axis(2)}
    else{
     axis(2)
     if(length(ylabs) == 1 & ylabs[1] == TRUE){ylabs<-sort(unique(y))}
     yxlabs <- rmdata(ylabs,y)$x
     yxlabs <- yxlabs - 0.5*diff(c(0,yxlabs))
#
#    next places it on the top outside
#
     if(really.missing(ylabslabs,missargs)){ylabslabs <- paste(ylabs)}
     axis(side = 4, at = yxlabs, labels = ylabslabs, tck=0.02, 
          cex=cex)
    }
    if(really.missing(yolabs,missargs)){
     axis(1)}
    else{
     axis(1)
     if(length(yolabs) == 1 & yolabs[1] == TRUE){yolabs<-sort(unique(yo))}
     yoxlabs <- rmdata(yolabs,yo)$x
#
#    next places it on the top outside
#
     if(really.missing(yolabslabs,missargs)){yolabslabs <- paste(yolabs)}
     axis(side = 3, at = yoxlabs, labels = yolabslabs, tck=0.02, cex=cex)
    }
    add <- 1
   }
   segments(x1,y1,x2,y2,lty=add)
   abline(h = 1, lty = 2)
  }
#
# Calculate polarization
#
  Q <- length(gx)
  se <- (n+m+1)/(3*n*m)
  cpo <- cumsum(po)
#
  hlf <- match(TRUE,cpo >= 0.5) -1
  g11 <- c(po[1:hlf], 0.5 - cpo[hlf])
  g13 <- gx[1:(hlf + 1)]
  g14 <- c(0, cpo[1:hlf])
  rpl <- 2 * sum(g13 * g11 * (1 - 2 * g14 - g11)) / sum(g11*g13) - 1	
#
  g12 <- c(cpo[hlf + 1] - 0.5, po[(hlf + 2):Q])
  g23 <- gx[(hlf + 1):Q]
  g24 <- c(0.5, cpo[(hlf + 1):(Q - 1)])
  rpu <- -2 * sum(g23 * g12 * (1 - 2 * g24 - g12)) / sum(g12*g23) - 1	
#
  rp <- 0.5 * (rpl + rpu)
#
# CI
#
gx.ci <- NULL
if(ci){
  vargx <- gx*gx*((1-pix)/(m*pix) + (1-po)/(n*po))
  gx.ci <- list(l=gx-1.96*sqrt(vargx), u=gx+1.96*sqrt(vargx))
  if(graph){
   segments(c(0,x),gx.ci$l,c(x,1),gx.ci$l,lty=2)
   segments(c(0,x),gx.ci$u,c(x,1),gx.ci$u,lty=2)
  }
}
if(!quiet){
  list(x=x,y=gx, rp=c(rp-1.96*sqrt(se),rp,rp+1.96*sqrt(se)),
     ci=gx.ci,rpl=rpl,rpu=rpu,breaks=rmd$breaks,chisq=chisq,entropy=entropy)
}else{
  invisible(
    list(x=x,y=gx, rp=c(rp-1.96*sqrt(se),rp,rp+1.96*sqrt(se)),
     ci=gx.ci,rpl=rpl,rpu=rpu,breaks=rmd$breaks,chisq=chisq,entropy=entropy)
  )
}

}         
"entropy" <- function(x,xo=FALSE){
#
# 	Calculate the entropy of the relative distribution
#
	if(missing(xo)){
	  egv <- x$y[x$y > 0]
	  dgv <- diff(c(0,x$x))
	  entropy <- sum(egv*log(egv)*dgv)
	}else{
          xout <- seq(0,1,length=1000)
          x <- approx(x=x$x,y=x$y,xout=xout,rule=2)
          xo <- approx(x=xo$x,y=xo$y,xout=xout,rule=2)
	  egv <- x$y[x$y > 0]
	  egvo <- xo$y[x$y > 0]
	  dgvo <- diff(c(0,xo$x))
	  entropy <- sum(egvo*log(egv)*dgvo)
	}
	entropy
}
really.missing <- function(x,missargs) {
 missargs[as.character(substitute(x))]   
}
#
# a function to change the weights to represent the probs
#
rdsamp <- function(targetvalues, samplevalues, 
                   targetvalueswgt=FALSE, samplevalueswgt=FALSE,
                   target=samplevalues, smooth=1){
#
#	targetvalues is the sample from distribution of the covariate
#	            that we wish to match to
#
#	samplevalues is the sample from distribution of the covariate
#	            that we start with
#
	n <- length(samplevalues)
	m <- length(targetvalues)
#
	sy <- sort(samplevalues)
	sy <- c(min(targetvalues,samplevalues),
                (sy[-1] - 0.5*diff(sy)),
                max(targetvalues,samplevalues))                  
	sy <- unique(sy)     
#
#	spb <- table(cut(samplevalues,breaks=sy,include.lowest = TRUE,labels=FALSE))
	spb <- hist(samplevalues,breaks=sy,plot=FALSE,include.lowest = TRUE)$counts
#
	if(missing(samplevalueswgt)){
#  	  spb <- table(cut(samplevalues,breaks=sy,include.lowest = TRUE,labels=FALSE))
	  spb <- hist(samplevalues,breaks=sy,plot=FALSE,include.lowest = TRUE)$counts
	}else{
  	  samplevalueswgt <- samplevalueswgt/mean(samplevalueswgt)
  	  xxx <- cut(samplevalues,breaks=sy,include.lowest = TRUE,labels=FALSE)
  	  spb <- rep(0,length=length(sy)-1)
  	  for (i in unique(xxx)){    
    	    spb[i] <- sum(samplevalueswgt[xxx==i],na.rm=TRUE)
  	  }
	}
#
# remove zero
#
	eee <- seq(along=spb)[spb < 1] 
	if(length(eee) > 0){
	if(eee[length(eee)] == length(spb)){
	  eee[length(eee)] <- eee[length(eee)]-1}
#	eee <- eee[eee < length(spb)]       
#
	sy <- sy[-(eee+1)]
	if(missing(samplevalueswgt)){
#	  spb <- table(cut(samplevalues,breaks=sy,include.lowest = TRUE,labels=FALSE))
	  spb <- hist(samplevalues,breaks=sy,plot=FALSE,include.lowest = TRUE)$counts
	}else{
	  samplevalueswgt <- samplevalueswgt/mean(samplevalueswgt)
	  xxx <- cut(samplevalues,breaks=sy,include.lowest = TRUE,labels=FALSE)
	  spb <- rep(0,length=length(sy)-1)
	  for (i in unique(xxx)){    
	    spb[i] <- sum(samplevalueswgt[xxx==i],na.rm=TRUE)
	  }
	}
	}
#
# remove small cells
#
        if(smooth > 4){
        ggg <- seq(5,smooth,length=4)
#
#       Pool successively
#
        for ( px in ggg ){

	  eee <- seq(along=spb)[spb <= px] 
#
	  if(length(eee)>0){
	    if(eee[length(eee)] == length(spb)){
		  eee[length(eee)] <- eee[length(eee)]-1}
	     sy <- sy[-(eee+1)]
	  }
	  if(missing(samplevalueswgt)){
#	    spb <- table(cut(samplevalues,breaks=sy,include.lowest = TRUE,labels=FALSE))
	    spb <- hist(samplevalues,breaks=sy,plot=FALSE,include.lowest = TRUE)$counts
	  }else{
	    samplevalueswgt <- samplevalueswgt/mean(samplevalueswgt)
	    xxx <- cut(samplevalues,breaks=sy,include.lowest = TRUE,labels=FALSE)
	    spb <- rep(0,length=length(sy)-1)
	    for (i in unique(xxx)){    
	      spb[i] <- sum(samplevalueswgt[xxx==i],na.rm=TRUE)
	    }
	  }
        }}
#
        spb <- spb / n
#
	if(missing(samplevalueswgt)){
#	  tpb <- table(cut(targetvalues,breaks=sy,include.lowest = TRUE,labels=FALSE))
	  tpb <- hist(targetvalues,breaks=sy,plot=FALSE,include.lowest = TRUE)$counts
	}else{
	  targetvalueswgt <- targetvalueswgt/mean(targetvalueswgt)
	  xxx <- cut(targetvalues,breaks=sy,include.lowest = TRUE,labels=FALSE)
	  tpb <- rep(0,length=length(sy)-1)
	  for (i in unique(xxx)){    
	    tpb[i] <- sum(targetvalueswgt[xxx==i],na.rm=TRUE)
	  }
	}
	tpb <- as.vector(tpb)+0.5
	tpb <- tpb/sum(tpb)
	xx <- cumsum(spb)    
	yy <- tpb/spb    
	nx <- length(xx)  
	x1 <- rep(c(0,xx),rep(2,nx+1))[-1]
	y1 <- rep(yy,c(rep(2,nx-1),3))
	y2 <- rep(yy,c(rep(2,nx-1),4))[-1]             
	x2 <- c(rep(xx,rep(2,nx)),1)  
#
   spb <- as.vector(spb)
#
   relprobs <- tpb/spb
#
#  cbind(sy, spb, tpb, relprobs)
#
#  sy are the category cutpoints
#  spb are the proportions of the sample values in each category.
#  tpb are the proportions of the target values in each category.
#  relprobs is then the density ratio of the targetvalues 
#    to the samplevalues, with values defined by the categories
#
#   relprobs[match(samplevalues, sort(unique(targetvalues)))]/n
#
#   relprobs[cut(samplevalues,breaks=sy,include.lowest = TRUE)]
#
   relprobs[cut(target,breaks=sy,include.lowest = TRUE,labels=FALSE)]/length(target)
#
#  This returns the probabilities for each target value in samplevalues.
#  with the probabilities of the targetvalues.
#
}
rdeciles <- function(y, yo, ywgt=FALSE,yowgt=FALSE,
                     reference=yo,refwgt=yowgt,
                     binn=10,
                     show="none", decomp="locadd",
                     location="median", scale="IQR",
                     z=FALSE, zo=FALSE
                     ){
#
#	# INPUT variables
#	# required variables
#	y	sample from comparison distribution
#	yo	sample from reference distribution

#	# options
#	ywgt	weights on the comparison sample
#	yowgt	weights on the comparison sample
#	binn 	number of categories (default 10)
#	reference sample from an alternative
#                 reference distribution
#                 The reference cutpoints are taken
#		  from this instead of yo
#	refwgt  weights from an alternative
#               reference distribution

#	# OUTPUT list components
#	deciles	list for the CDF
#       x	relative ratio of the deciles
#       refcuts decile cutpoints from the reference 
#       yx	proportion of comparison within
#               each reference cutpoints
#               (If reference=yo this is close to x/binn).
#		This might deviate due to ties in the 
#               reference or comparison distributions.
#       yox	proportion of yo within
#               each reference cutpoints
#               (If reference=yo this should be 
#               close to 1/binn).
#		This might deviate due to ties in the 
#               reference or comparison distributions.
#
#       See reldist for other arguments
#
  rmd <- rmdata(y=y, yo=yo, ywgt=ywgt,yowgt=yowgt,
                show=show, decomp=decomp,
                location=location, scale=scale,
                z=z, zo=zo
               )
  y <- rmd$y
  yo <- rmd$yo
  ywgt <- rmd$ywgt
  yowgt <- rmd$yowgt
  x <- rmd$x
  n <- rmd$n
  m <- rmd$m
#
  if(missing(reference)){reference <- yo}
  if(missing(refwgt)){refwgt <- yowgt}
  sy <- seq(along = reference)
  ys <- sort.list(reference)
  if(!missing(refwgt)){refwgt <- refwgt[ys]}
  reference <- reference[ys]
#
  if(missing(refwgt)){
   cdfref <- sy/m
  }else{
   refwgt <- refwgt/sum(refwgt)
   cdfref <- cumsum(refwgt)
  }
  refcuts <- seq(0, 1, length = binn + 1)[c(-1,-(binn+1))]
  refcuts <- c(min(y,yo,reference)-0.01,
               approx(x=cdfref,y=reference,xout=refcuts)$y,
               max(y,yo,reference)+0.01)
#
  if(missing(ywgt)){
#  yx <- table(cut(y, breaks = refcuts, include.lowest=TRUE,labels=FALSE))/m
   yx <- hist(y,breaks=refcuts,plot=FALSE,include.lowest = TRUE)$counts/m
  }else{
   ywgt <- ywgt/sum(ywgt)
   xxx <- cut(y, breaks = refcuts,include.lowest=TRUE,labels=FALSE)
   yx <- rep(0,length=binn)
   for (i in unique(xxx)){    
    yx[i] <- sum(ywgt[xxx==i],na.rm=TRUE)
   }
  }
#
  if(missing(yowgt)){
#  yox <- table(cut(yo, breaks = refcuts, include.lowest=TRUE,labels=FALSE))/n
   yox <- hist(yo,breaks=refcuts,plot=FALSE,include.lowest = TRUE)$counts/n
  }else{
   yowgt <- yowgt/sum(yowgt)
   xxx <- cut(yo, breaks = refcuts,include.lowest=TRUE,labels=FALSE)
   yox <- rep(0,length=binn)
   for (i in unique(xxx)){    
    yox[i] <- sum(yowgt[xxx==i],na.rm=TRUE)
   }
  }
  ratio <- yx/yox
#
# Output
#
  list(x=as.numeric(ratio),refcuts=refcuts,yx=yx,yox=as.numeric(yox))
}
#"relasp"<- function(y,yo,binn=100,boundary=0.05){
#  y <- sort(y)
#  m <- length(y)
#  y <- c(min(y)-y[1:(boundary*m)],y,2*max(y)-y[((1-boundary)*m):m])
#  m <- length(y)
#  yo <- sort(yo)
#  n <- length(yo)
#  yo <- c(min(yo)-yo[1:(0.1*n)],yo,2*max(yo)-yo[(0.9*n):n])
#  n <- length(yo)
#  write(binn, "cao.in", append = FALSE)
#  write(n+m, "cao.in", append = TRUE)
#  rdata <- cbind(0,c(y,yo),1,rep(c(1,-1),c(n,m)))
#  write(t(rdata), "cao.in", append = TRUE, ncol = 4)
##
## Call converted Pascal program
##
#  aaa <- system("cao < cao.in > cao.out")
#  out <- matrix(as.numeric(scan("cao.out")), ncol = 2,byrow=TRUE)
## out[is.na(out[, 2]), 2] <- 0.1
#  out[is.na(out[, 2]), 2] <- mean(as.numeric(out[!is.na(out[, 2]), 2]))
#  dimnames(out) <- list(NULL,c("x","h(x)"))
#  out
#}
resplot <- function(x, standardize=TRUE,
                    xlab="Gaussian Cumulative Proportion", ...){
  if(standardize){
   stdres <- ( x - mean(x, na.rm=TRUE) ) / sqrt(var(x, na.rm=TRUE)) 
  }else{
   stdres <- x[!is.na(x)]
  }
  xrange <- range(c(0,2,reldist(y=pnorm(q=stdres), graph=FALSE, ...)$y))
  reldist(y=pnorm(q=stdres),
          ylim=xrange,
          xlab=xlab, ...)
}
