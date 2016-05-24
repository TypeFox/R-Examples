shadeDist <- function(xshade=NULL, ddist="dnorm", parm1=NULL, parm2=NULL,
                  lower.tail=TRUE, xlab=NULL, xmin=NULL, xmax=NULL, xtic=TRUE,
                  digits.prob=4, digits.xtic=3, is.discrete=NULL,
                  additional.x.range=NULL, main=NULL, col=c("black","pink"), lwd=5, ...) {
  # Plots a continuous probability density function (pdf) and shades in a region.
  # `xshade' can be a scalar or a vector of size 2.
  # When `xshade' is a scalar and lower.tail is TRUE,
  #    the area from -Inf to xshade is shaded.
  # When `xshade' is a scalar and lower.tail is FALSE,
  #    the area from xshade to +Inf is shaded.
  # When `xshade' is a vector and lower.tail is TRUE,
  #    the tail probabilities are shaded.
  # When `xshade' is a vector and lower.tail is FALSE,
  #    the area from xshade[1] to xshade[2] is shaded.
  # `ddist' is the pdf and must be in quotes.  Some choices are
  #    dnorm, dexp, dcauchy, dlaplace, dt, dt.nc, dchisq, df, df.nc, dunif.
  #    If the random variable is a sample proportion, then set ddist="dprop" and use same parameters from "dbinom".
  # `parm1' and `parm2' are numeric scalars or vectors of the parameters of `distA'
  #    (excluding the first argument), where NULL is the default value.
  # If the distribution needs more than 2 parameters, then assign `parm1'
  #    to be a vector of all (as many as 5) parameters, and keep `parm2' equal to NULL.
  #    For example, to shade a noncentral F pdf with 7 and 9 df and ncp=2.5 in the
  #       right tail beginning at 6, type:
  #       > shadeDist( 6, 'df.nc', c(7,9,2.5), NULL, F, 0, 8 )
  # To avoid ambiguity with discrete distributions, avoid boarder values with `xshade'.
  #    For example, if X ~ Binomial(n=10, p=0.4), then to graph P(X>2) type:
  #    > shadeDist( 2.5, 'dbinom', 10, 0.4, F )
  # `xlab' is the label on the x-axis.  If NULL, `xlab' may be automatically set to `X',
  #    `Z', `T', `F' or `X^2', according to `ddist'.
  # `xmin' and `xmax' represent the limiting values on the x-axis.
  #    Reasonable values of `xmin' and `xmax' are generated when originally set to NULL.
  # `digits.prob' is the number of significant digits listed in the probability.
  # `digits.xtic' is the number of significant digits listed on the x-axis.
  # If `xtic' is TRUE, then the numbers listed on the x-axis include the median and `xshade'.
  # If `xtic' is FALSE, then the default numbers are listed on the x-axis.
  # If `xtic' is a vector, then the numbers listed on the x-axis consist of `xtic'.
  # `is.discrete' is logical, indicating whether or not the distribution is discrete.
  #    If NULL, then the macro attempts to assign the correct logical value.
  # `additional.x.range' is a vector of TWO additional `x' values for evaluating the function.
  #    This is useful when the plot is too sparse.  This is ignored if ddist="dprop".
  # `main' is the main title.  If NULL, the probability is displayed.
  #    If no title is desired, set `main' to an empty string ''.
  # `col' represents the color of the density function as col[1],
  #    and the color of the shading as col[2].
  # `lwd' is line width for discrete densities.
  # ...: optional arguments to `plot'.
 if (ddist!="dprop")    {
  num.grid.points <- 10001;    ylimchisq1 <- 1.2
  options(warn = -1);   
  if (!is.character(ddist))  
     stop("'ddist' must be a probability density function in character format.")
  if (substring(ddist,1,1)!="d") stop("'ddist' must be a probability density function in character format.")
  if (!is.function(try(match.fun(ddist),silent=TRUE))) 
     stop("'ddist' must be a probability density function in character format.")
  if (!is.numeric(xshade) & !is.null(xshade))  stop("'xshade' must be numeric or NULL.")
  if (!is.numeric(parm1) & !is.null(parm1))  stop("'parm1' must be numeric or NULL.")
  if (!is.numeric(parm2) & !is.null(parm2))  stop("'parm2' must be numeric or NULL.")
  if (!is.numeric(xmin) & !is.null(xmin))  stop("'xmin' must be numeric or NULL.")
  if (!is.numeric(xmax) & !is.null(xmax))  stop("'xmax' must be numeric or NULL.")
  if (!is.numeric(xtic) & !is.logical(xtic))  stop("'xtic' must be numeric or logical.")
  if (!is.numeric(additional.x.range) & !is.null(additional.x.range))  stop("'additional.x.range' must be numeric or NULL.")
  if (!is.numeric(digits.prob))  stop("'digits.prob' must be numeric.")
  if (!is.numeric(digits.xtic))  stop("'digits.xtic' must be numeric.")
  if (!is.numeric(lwd))  stop("'lwd' must be numeric.")
  if (!is.logical(is.discrete) & !is.null(is.discrete))  stop("'is.discrete' must be logical or NULL.")
  if (!is.logical(lower.tail))  stop("'lower.tail' must be logical.")
  if (!is.character(xlab) & !is.null(xlab))  stop("'xlab' must be character or NULL.")
  if (!is.character(main) & !is.null(main))  stop("'main' must be character or NULL.")
  all.discrete = c("dbinom", "dgeom", "dhyper", "dnbinom", "dpois")
  if (!is.null(is.discrete)) discrete.pdf=is.discrete
  if (is.null(is.discrete)) discrete.pdf=(ddist %in% all.discrete)
  if (is.null(xlab)) {   xlab="X"
    if (ddist=="dnorm") {    if (is.null(c(parm1,parm2))) xlab="Z"
       if (!is.null(parm1) & is.null(parm2)) { if (parm1==0) xlab="Z" }
       if (!is.null(parm1) & !is.null(parm2)) { if (parm1==0 & parm2==1) xlab="Z" }    }
    if (xlab=="Z") {parm1=0; parm2=1}
    if (ddist=="dnorm" & is.null(parm2)) parm2==1
    temp <- c(parm1, parm2)
    if (ddist=="dt" & length(temp)==1) xlab=paste(c("T_",temp[1]),sep="",collapse="")
    if ( (ddist=="dt.nc" & length(temp)==2) | (ddist=="dt" & length(temp)==2) ) xlab=paste(c("T_",temp[1],",",temp[2]),sep="",collapse="")
    if (ddist=="dchisq" & length(temp)==1) xlab=paste(c("CHI-SQUARE_",temp[1]),sep="",collapse="")
    if (ddist=="dchisq" & length(temp)==2) xlab=paste(c("CHI-SQUARE_",temp[1],",",temp[2]),sep="",collapse="")
    if (ddist=="df" & length(temp)==2) xlab=paste(c("F_",temp[1],",",temp[2]),sep="",collapse="")
    if ( (ddist=="df.nc" & length(temp)==3) | (ddist=="df" & length(temp)==3) )
               xlab=paste(c("F_",temp[1],",",temp[2],",",temp[3]),sep="",collapse="") }
  parm1 <- c(parm1, parm2);  temp.vec <- getMinMax(xmin, xmax, ddist, parm1)
  if (is.null(xshade)) {xshade <- temp.vec$medianA}
  if (length(xshade)!=1 & length(xshade)!=2) stop("xshade must have length 1 or 2.")
  if (length(xshade)==2)
    {if (xshade[1]>xshade[2]) {temp <- xshade[1]; xshade[1] <- xshade[2]; xshade[2] <- temp}}
  if (is.null(xmin) || is.null(xmax))   {
    xmin <- temp.vec$xmin ;    xmax <- temp.vec$xmax
    xmin.temp <- min( xmin, xshade[1]-(xmax-xmin)*0.1 )
    xmax.temp <- max( xmax, xshade[length(xshade)]+(xmax-xmin)*0.1 )
    xmin <- xmin.temp ;   xmax <- xmax.temp  }
  else {  if (xmin>xmax) { temp <- xmin; xmin <- xmax; xmax <- temp }  }
  pdist <- paste("p",substring(ddist,2),sep="")
  if (!is.function(try(match.fun(pdist),silent=TRUE))) 
     stop("'ddist' must have corresponding cumulative distribution function 'pdist'.")
  qdist <- paste("q",substring(ddist,2),sep="")
  if (pdist=="pf.nc") pdist <- "pf" ;     if (pdist=="pt.nc") pdist <- "pt"
  if (!is.function(try(match.fun(qdist),silent=TRUE))) 
     stop("'ddist' must have correspoding quantile function 'qdist'.")
  f <- function(z)   {
                if (length(parm1)==0)
                   {return(match.fun(ddist)(z))}
                if (length(parm1)==1)
                   {return(match.fun(ddist)(z, parm1[1]))}
                if (length(parm1)==2)
                   {return(match.fun(ddist)(z, parm1[1], parm1[2]))}
                if (length(parm1)==3)
                   {return(match.fun(ddist)(z, parm1[1], parm1[2], parm1[3]))}
                if (length(parm1)==4)
                   {return(match.fun(ddist)(z, parm1[1], parm1[2], parm1[3], parm1[4]))}
                if (length(parm1)==5)
                   {return(match.fun(ddist)(z, parm1[1], parm1[2], parm1[3],
                                               parm1[4], parm1[5]))}     }
  g <- function(z, lower.tail=TRUE)   {
                if (length(parm1)==0)
                   {return(match.fun(pdist)(z, lower.tail=lower.tail))}
                if (length(parm1)==1)
                   {return(match.fun(pdist)(z, parm1[1], lower.tail=lower.tail))}
                if (length(parm1)==2)
                   {return(match.fun(pdist)(z, parm1[1], parm1[2], lower.tail=lower.tail))}
                if (length(parm1)==3)
                   {return(match.fun(pdist)(z, parm1[1], parm1[2], parm1[3], lower.tail=lower.tail))}
                if (length(parm1)==4)
                   {return(match.fun(pdist)(z, parm1[1], parm1[2], parm1[3], parm1[4], lower.tail=lower.tail))}
                if (length(parm1)==5)
                   {return(match.fun(pdist)(z, parm1[1], parm1[2], parm1[3],
                                               parm1[4], parm1[5], lower.tail=lower.tail))}     }
  h <- function(z)   {
                if (length(parm1)==0)
                   {return(match.fun(qdist)(z))}
                if (length(parm1)==1)
                   {return(match.fun(qdist)(z, parm1[1]))}
                if (length(parm1)==2)
                   {return(match.fun(qdist)(z, parm1[1], parm1[2]))}
                if (length(parm1)==3)
                   {return(match.fun(qdist)(z, parm1[1], parm1[2], parm1[3]))}
                if (length(parm1)==4)
                   {return(match.fun(qdist)(z, parm1[1], parm1[2], parm1[3], parm1[4]))}
                if (length(parm1)==5)
                   {return(match.fun(qdist)(z, parm1[1], parm1[2], parm1[3],
                                               parm1[4], parm1[5]))}     }
  x <- xmin + (xmax-xmin)*(0:num.grid.points)/num.grid.points; options(warn = -1)
  options(warn = 0)
  if (discrete.pdf) {xmin <- ceiling(xmin); xmax<-floor(xmax); x <- xmin:xmax }
  if (length(xshade)==2 & xshade[1]>xshade[2]) {
     temp <- xshade[1];   xshade[1] <- xshade[2];   xshade[2] <- temp }
  if (xshade[1]<xmin || xshade[length(xshade)]>xmax)
     stop("xshade must be between xmin and xmax.")
  if (lower.tail  & length(xshade)==1)
     {the.prob <- g(xshade);    x1 <- x[ x<=xshade ]     }
  if (!lower.tail & length(xshade)==1)
     {the.prob <- g(xshade, lower.tail=FALSE);  x1 <- x[ x>xshade ]   }
  if (!lower.tail & length(xshade)==2)
     {the.prob <- g(xshade[2])-g(xshade[1], lower.tail=TRUE)
      x1 <- x[ xshade[1]<x & x<=xshade[2] ]   }
  if (lower.tail  & length(xshade)==2)
     {the.prob <- g(xshade[1])+g(xshade[2], lower.tail=FALSE)
      x1 <- x[ x<=xshade[1] ] ;   x2 <- x[ x>xshade[2] ]    }
  maxA <- ifelse( ((ddist %in% c("dchisq", "df")) && parm1[1]==1),
                    min(ylimchisq1, max(f(x),na.rm=TRUE)), max(f(x))  )
  ylim1 <- c(min(f(x),na.rm=TRUE), maxA)
  if (is.null(main))  main=
      paste(c("Probability is", prettyNum(the.prob,digits=digits.prob)), collapse=" ")
  if (discrete.pdf)  {
    plot(x, f(x), xlim=c(xmin, xmax), ylim=ylim1, type="h", xlab=xlab,
         ylab="probability  density  function", xaxt="n", main=main, col=col[1], lwd=lwd, ...)
    curve(f, min(x1), max(x1), (max(x1)-min(x1)+1), add=TRUE, type="h", xaxt="n", col=col[2], lwd=lwd, ...)
    if (lower.tail==TRUE & length(xshade)==2)  {
        curve(f, min(x2), max(x2), (max(x2)-min(x2)+1), add=TRUE, type="h", xaxt="n", col=col[2], lwd=lwd, ...) }  }
  else   {
    plot(x1, f(x1), xlim=c(xmin,xmax), ylim=ylim1, type="h", xlab=xlab,
         ylab="probability  density  function", xaxt="n", main=main, col=col[2], ...)
    if (lower.tail==TRUE & length(xshade)==2)  {
        curve(f, min(x2), max(x2), num.grid.points, add=TRUE, type="h", xaxt="n", col=col[2])  }
    curve(f, xmin, xmax, num.grid.points, add=TRUE, type="p", xaxt="n", pch=20, cex=0.5, col=col[1])  
    if (length(additional.x.range)>=2) 
          curve(f, min(additional.x.range), max(additional.x.range), 
          num.grid.points, add=TRUE, type="p", xaxt="n", pch=20, cex=0.5, col=col[1])  
    }
  if (!is.logical(xtic))    axis(1, at=xtic, labels=xtic)
  if (is.logical(xtic))
     { if (!xtic) axis(1);
       if (xtic) {
               if (length(xshade)==1) xtic0=prettyNum(c(h(0.5), 2*h(0.5)-xshade, xshade),digits=digits.xtic)
               if (length(xshade)==2) xtic0=prettyNum(c(h(0.5), xshade),digits=digits.xtic)
               if (ddist %in% c("dchisq","df")) xtic0=prettyNum(c(h(0.5), xshade, 0),digits=digits.xtic)
               axis(1, at=xtic0, labels=xtic0)    } }
 }   # End of plotting for random variables other than a sample proportion
 if (ddist=="dprop")    
    shadePhat(xshade=xshade, size=parm1, prob=parm2, lower.tail=lower.tail, xmin=xmin, xmax=xmax,
                  xlab=xlab, xtic=xtic, digits.prob=digits.prob, digits.xtic=digits.xtic, 
                  main=main, col=col, lwd=lwd, ...)   
} 


shadePhat <- function(xshade=NULL, size=1, prob=0.5, lower.tail=TRUE, xmin=0, xmax=1,
                  xlab=expression(hat(p)), xtic=TRUE,
                  digits.prob=4, digits.xtic=3, main=NULL, col=c("black","pink"), lwd=5, ...) {
  # Plots a continuous probability density function (pdf) and shades in a region.
  # `xshade' can be a scalar or a vector of size 2.
  # When `xshade' is a scalar and lower.tail is TRUE,
  #    the area from 0 to xshade is shaded.
  # When `xshade' is a scalar and lower.tail is FALSE,
  #    the area from xshade to 1 is shaded.
  # When `xshade' is a vector and lower.tail is TRUE,
  #    the tail probabilities are shaded.
  # When `xshade' is a vector and lower.tail is FALSE,
  #    the area from xshade[1] to xshade[2] is shaded.
  # `size' and `prob' are scalars, corresponding to the parameters in `dbinom'.
  # To avoid ambiguity with discrete distributions, avoid boarder values with `xshade'.
  # `xmin' and `xmax' represent the limiting values on the x-axis.
  # `xlab' is the label on the x-axis.
  # `digits.prob' is the number of significant digits listed in the probability.
  # `digits.xtic' is the number of significant digits listed on the x-axis.
  # If `xtic' is TRUE, then the numbers listed on the x-axis consist of `p' and `xshade'.
  # If `xtic' is FALSE, then the default numbers are listed on the x-axis.
  # If `xtic' is a vector, then the numbers listed on the x-axis consist of `xtic'.
  # `main' is the main title.  If NULL, the probability is displayed.
  #    If no title is desired, set `main' to an empty string ''.
  # `col' represents the color of the density function as col[1],
  #    and the color of the shading as col[2].
  # `lwd' is line width for discrete densities.
  # ...: optional arguments to `plot'.

  if (!is.numeric(xshade) & !is.null(xshade))  stop("'xshade' must be numeric or NULL.")
  if (!is.numeric(size))  stop("'size' must be numeric.")
  if (!is.numeric(prob))  stop("'prob' must be numeric.")
  if (!is.numeric(xmin) & !is.null(xmin))  stop("'xmin' must be numeric or NULL.")
  if (!is.numeric(xmax) & !is.null(xmax))  stop("'xmax' must be numeric or NULL.")
  if (!is.logical(lower.tail))  stop("'lower.tail' must be logical.")
  if (!is.numeric(lwd))  stop("'lwd' must be numeric.")
  if (!is.numeric(digits.prob))  stop("'digits.prob' must be numeric.")
  if (!is.numeric(digits.xtic))  stop("'digits.xtic' must be numeric.")
  if (!is.numeric(xtic) & !is.logical(xtic))  stop("'xtic' must be numeric or logical.")
  if (!is.character(main) & !is.null(main))  stop("'main' must be character or NULL.")
  if (!(is.expression(xlab) | is.character(xlab)) & !is.null(xlab))  stop("'xlab' must be an expression, character, or NULL.")
  if (is.null(xshade)) xshade=0.5
  if (is.null(size))   size=1
  if (is.null(prob))   prob=0.5
  if (is.null(xmin))   xmin=0
  if (is.null(xmax))   xmax=1
  if (is.null(xlab))   xlab=expression(hat(p))
  num.grid.points <- 10001  
  options(warn = -1); n=size; p=prob
  if (is.null(xshade)) xshade <- p
  if (length(xshade)!=1 & length(xshade)!=2) stop("xshade must have length 1 or 2.")
  if (length(xshade)==2)
    {if (xshade[1]>xshade[2]) {temp <- xshade[1]; xshade[1] <- xshade[2]; xshade[2] <- temp}}
  if (xmin>xmax) { temp <- xmin; xmin <- xmax; xmax <- temp }
  x <- xmin + (xmax-xmin)*(0:num.grid.points)/num.grid.points; options(warn = -1)
  options(warn = 0)
  f=function(x){dbinom(x,n,p)} ;   g=function(prop){dbinom(round(prop*n),n,p)}
  xmin1 <- ceiling(xmin*n); xmax1<-floor(xmax*n); x <- xmin1:xmax1
  if (xshade[1]<xmin || xshade[length(xshade)]>xmax)
     stop("xshade must be between xmin and xmax.")
  if (lower.tail & length(xshade)==1)
     {the.prob <- pbinom(xshade*n,n,p);    x1 <- x[ x<=xshade*n ]     }
  if (!lower.tail & length(xshade)==1)
     {the.prob <- 1-pbinom(xshade*n,n,p);  x1 <- x[ x>xshade*n ]   }
  if (!lower.tail & length(xshade)==2)
     {the.prob <- pbinom(xshade[2]*n,n,p)-pbinom(xshade[1]*n,n,p)
      x1 <- x[ xshade[1]*n<x & x<=xshade[2]*n ]   }
  if (lower.tail & length(xshade)==2)
     {the.prob <- 1-pbinom(xshade[2]*n,n,p)+pbinom(xshade[1]*n,n,p)
      x1 <- x[ x<=xshade[1]*n ] ;   x2 <- x[ x>xshade[2]*n ]    }
  ylim1 <- c(0, max(f(x),na.rm=TRUE))
  if (is.null(main))   main=
       paste(c("Probability is", prettyNum(the.prob,digits=digits.prob)), collapse=" ")
  plot(x/n, f(x), xlim=c(xmin, xmax), ylim=ylim1, type="h", xlab=xlab,
       ylab="probability  density  function", xaxt="n", main=main, col=col[1], lwd=lwd, ...)
  curve(g, min(x1/n), max(x1/n), max(x1)-min(x1)+1, add=TRUE, type="h", xaxt="n", col=col[2], lwd=lwd)
  if (lower.tail & length(xshade)==2)  {
      curve(g, min(x2/n), max(x2/n), max(x2)-min(x2)+1, add=TRUE, type="h", xaxt="n", col=col[2], lwd=lwd) }
  if (is.logical(xtic))
     { if (!xtic) axis(1);
       if (xtic) {
           if (length(xshade)==1) xtic=prettyNum(c(xshade, 2*p-xshade, p, xmin, xmax),digits=digits.xtic)
           if (length(xshade)==2) xtic=prettyNum(c(xshade, p, xmin, xmax),digits=digits.xtic)
           axis(1, at=xtic, labels=xtic) } }
  if (!is.logical(xtic))    axis(1, at=xtic, labels=xtic)
} 
                                                                                                                       

plotDist <- function(distA="dnorm", parmA1=NULL, parmA2=NULL,
             distB=NULL, parmB1=NULL, parmB2=NULL, distC=NULL, parmC1=NULL, parmC2=NULL,
             xlab=NULL, xmin=NULL, xmax=NULL, col=c("black","red","darkgreen"), 
             is.discrete=NULL, additional.x.range=NULL, lwd=5, ...) {
  # Plots one, two, or three functions.
  # When plotting the pdf, some choices for `distA', `distB' and `distC' are
  #    dnorm, dexp, dcauchy, dlaplace, dt, dt.nc, dchisq, df, df.nc, dunif, dbinom, 
  #    dgeom, dpois, dhyper, dnbinom.
  # When plotting the cdf, some choices for `distA', `distB' and `distC' are
  #    pnorm, pexp, pcauchy, plaplace, pt, pchisq, pf, punif, pbinom, 
  #    pgeo, ppois, phyper, pnbinom.
  # The values of `distA', `distB' and `distC' must be in quotes if `xmin' or `xmax' is NULL.
  # `parmA1' and `parmA2' are numeric scalars or vectors of the parameters of `distA'
  #    (excluding the first argument), where NULL is the default value;
  #    likewise for the other two functions.
  # `col' may be a scalar or vector, and specifies the colors of the plotted functions.
  #    Type `colors()' for selections.
  # If only one function is to be plotted, then use `distA'.
  # If a distribution, say `distA', needs more than 2 parameters, then assign `parmA1'
  #    to be a vector of all (as many as 5) parameters, and keep `parmA2' equal to NULL.
  #    For example, to plot a noncentral F pdf with 7 and 9 df and ncp=2.5, type:
  #    > plot2.dist( df.nc, c(7,9,2.5) )
  # `xmin' and `xmax' represent the limiting values on the x-axis.
  #    Reasonable values of `xmin' and `xmax' are attempted when originally set to NULL.
  # `xlab' is the label on the x-axis.  If NULL, `xlab' may be automatically set to `X',
  #    `Z', `T', `F' or `X^2', according to `distA', `distB' and `distC'.
  # `is.discrete' is a vector with 1, 2, or 3 logical values, indicating whether or not 
  #    `distA', `distB', and `distC' are discrete.
  #    If NULL, then the macro attempts to assign the correct logical value(s).
  # `additional.x.range' is a vector of TWO additional `x' values for evaluating the function.
  #    This is useful when the plot is too sparse.
  # `lwd' is line width for discrete densities.
  num.grid.points <- 10001;    options(warn= -1);   ylimchisq1 <- 1.2
  if (!is.character(distA))
     stop("'distA' must be a probability density function or cumulative distribution function in character format.")
  if (!(substring(distA,1,1) %in% c("d","p"))) stop("'distA' must be a probability density function or cumulative distribution in character format.")
  if (!is.function(try(match.fun(distA),silent=TRUE)))
     stop("'distA' must be a probability density function or cumulative distribution in character format.")
  if (!is.null(distB)) {
     if (!is.character(distB)) 
        stop("'distB' must be a probability density function or cumulative distribution function in character format or NULL.")
     if (!(substring(distB,1,1) %in% c("d","p"))) stop("'distB' must be a probability density function or cumulative distribution in character format or NULL.")
     if (!is.function(try(match.fun(distB),silent=TRUE)))
        stop("'distB' must be a probability density function or cumulative distribution in character format or NULL.")  }
  if (!is.null(distC)) {
     if (!is.character(distC)) 
        stop("'distC' must be a probability density function or cumulative distribution function in character format or NULL.")
     if (!(substring(distC,1,1) %in% c("d","p"))) stop("'distC' must be a probability density function or cumulative distribution in character format or NULL.")
     if (!is.function(try(match.fun(distC),silent=TRUE)))
        stop("'distC' must be a probability density function or cumulative distribution in character format or NULL.")  }
  if (!is.numeric(parmA1) & !is.null(parmA1))  stop("'parmA1' must be numeric or NULL.")
  if (!is.numeric(parmA2) & !is.null(parmA2))  stop("'parmA2' must be numeric or NULL.")
  if (!is.numeric(parmB1) & !is.null(parmB1))  stop("'parmB1' must be numeric or NULL.")
  if (!is.numeric(parmB2) & !is.null(parmB2))  stop("'parmB2' must be numeric or NULL.")
  if (!is.numeric(parmC1) & !is.null(parmC1))  stop("'parmC1' must be numeric or NULL.")
  if (!is.numeric(parmC2) & !is.null(parmC2))  stop("'parmC2' must be numeric or NULL.")
  if (!is.numeric(xmin) & !is.null(xmin))  stop("'xmin' must be numeric or NULL.")
  if (!is.numeric(xmax) & !is.null(xmax))  stop("'xmax' must be numeric or NULL.")
  if (!is.numeric(additional.x.range) & !is.null(additional.x.range))  stop("'additional.x.range' must be numeric or NULL.")
  if (!is.numeric(lwd))  stop("'lwd' must be numeric.")
  if (!is.character(xlab) & !is.null(xlab))  stop("'xlab' must be character or NULL.")
  all.discrete = c("dbinom", "dgeom", "dhyper", "dpois", "dnbinom", 
                   "pbinom", "pgeom", "phyper", "ppois", "pnbinom")
  if (!is.null(is.discrete)) discrete.pdfA=is.discrete[1]
  discrete.pdfB=FALSE; discrete.pdfC=FALSE
  if (is.null(is.discrete)) discrete.pdfA=(distA %in% all.discrete)
  if (!is.null(distB))
     {  if (!is.null(is.discrete)) discrete.pdfB=is.discrete[2]
        if (is.null(is.discrete)) discrete.pdfB=(distB %in% all.discrete) }
  if (!is.null(distC))
     {  if (!is.null(is.discrete)) discrete.pdfC=is.discrete[3]
        if (is.null(is.discrete)) discrete.pdfC=(distC %in% all.discrete) }
  if (is.null(xlab)) {   xlab="X"
    if ((distA=="dnorm" || distA=="pnorm") & is.null(distB) & is.null(distC)) {
       if (is.null(c(parmA1,parmA2))) xlab="Z"
       if (!is.null(parmA1) & is.null(parmA2)) { if (parmA1==0) xlab="Z" }
       if (!is.null(parmA1) & !is.null(parmA2)) { if (parmA1==0 & parmA2==1) xlab="Z" }    }
    temp = c(parmA1, parmA2)
    if (distA=="dt" || (distA=="pt" & length(temp)==1)) {
       xlab=paste(c("T_",temp[1]),sep="",collapse="")
       if ( !is.null(distB) ) { xlab="X"
           if ( length( union( union(distA,distB), union(distA,distC) ) )==1 ) xlab="T" }  }
    if (distA=="dt.nc" | (distA=="pt" & length(temp)==2) | (distA=="dt" & length(temp)==2)) {
       xlab=paste(c("T_",temp[1],",",temp[2]),sep="",collapse="")
       if ( !is.null(distB) ) { xlab="X"
           if ( length( union( union(distA,distB), union(distA,distC) ) )==1 ) xlab="T" }  }
    if (distA=="dchisq" || distA=="pchisq") {
       xlab=paste(c("CHI-SQUARE_",parmA1[1]),sep="",collapse="")
       if (length(temp)==2) xlab=paste(c("CHI-SQUARE_",temp[1],",",temp[2]),sep="",collapse="")
       if ( !is.null(distB) ) { xlab="X"
           if ( length( union( union(distA,distB), union(distA,distC) ) )==1 ) xlab="CHI-SQUARE" }  }
    if (distA=="df" || (distA=="pf" & length(temp)==2)) {
       xlab=paste(c("F_",temp[1],",",temp[2]),sep="",collapse="")
       if ( !is.null(distB) ) { xlab="X"
           if ( length( union( union(distA,distB), union(distA,distC) ) )==1 ) xlab="F" }  }
    if (distA=="df.nc" || (distA=="pf" & length(temp)==3) || (distA=="df" & length(temp)==3) ) {
       xlab=paste(c("F_",temp[1],",",temp[2],",",temp[3]),sep="",collapse="")
       if ( !is.null(distB) ) { xlab="X"
           if ( length( union( union(distA,distB), union(distA,distC) ) )==1 ) xlab="F" }  }     }
  parmA1 <- c(parmA1, parmA2); parmB1 <- c(parmB1, parmB2); parmC1 <- c(parmC1, parmC2)
  temp.vec <- getMinMax(xmin, xmax, distA, parmA1, NULL, distB, parmB1, NULL, distC, parmC1, NULL )
  xmin <- temp.vec$xmin ;    xmax <- temp.vec$xmax
  if (xmin>xmax) {temp <- xmin; xmin <- xmax; xmax <- temp}
  x <- xmin + (xmax-xmin)*(0:num.grid.points)/num.grid.points
  x0 <- ceiling(xmin):floor(xmax)
  if (length(col)==1) {col <- rep(col, 3)};  if (length(col)==2) {col <- c(col, col[1])}
  fA <- function(z)   {
                if (length(parmA1)==0)
                   {return(match.fun(distA)(z))}
                if (length(parmA1)==1)
                   {return(match.fun(distA)(z, parmA1[1]))}
                if (length(parmA1)==2)
                   {return(match.fun(distA)(z, parmA1[1], parmA1[2]))}
                if (length(parmA1)==3)
                   {return(match.fun(distA)(z, parmA1[1], parmA1[2], parmA1[3]))}
                if (length(parmA1)==4)
                   {return(match.fun(distA)(z, parmA1[1], parmA1[2], parmA1[3], parmA1[4]))}
                if (length(parmA1)==5)
                   {return(match.fun(distA)(z, parmA1[1], parmA1[2], parmA1[3],
                                               parmA1[4], parmA1[5]))}     }
  ylim1 <- NULL
  if (!is.null(distB)) {
     fB <- function(z)   {
                if (length(parmB1)==0)
                   {return(match.fun(distB)(z))}
                if (length(parmB1)==1)
                   {return(match.fun(distB)(z, parmB1[1]))}
                if (length(parmB1)==2)
                   {return(match.fun(distB)(z, parmB1[1], parmB1[2]))}
                if (length(parmB1)==3)
                   {return(match.fun(distB)(z, parmB1[1], parmB1[2], parmB1[3]))}
                if (length(parmB1)==4)
                   {return(match.fun(distB)(z, parmB1[1], parmB1[2], parmB1[3], parmB1[4]))}
                if (length(parmB1)==5)
                   {return(match.fun(distB)(z, parmB1[1], parmB1[2], parmB1[3],
                                                parmB1[4], parmB1[5]))}     }        }
  if (!is.null(distC)) {
     fC <- function(z)   {
                if (length(parmC1)==0)
                   {return(match.fun(distC)(z))}
                if (length(parmC1)==1)
                   {return(match.fun(distC)(z, parmC1[1]))}
                if (length(parmC1)==2)
                   {return(match.fun(distC)(z, parmC1[1], parmC1[2]))}
                if (length(parmC1)==3)
                   {return(match.fun(distC)(z, parmC1[1], parmC1[2], parmC1[3]))}
                if (length(parmC1)==4)
                   {return(match.fun(distC)(z, parmC1[1], parmC1[2], parmC1[3], parmC1[4]))}
                if (length(parmC1)==5)
                   {return(match.fun(distC)(z, parmC1[1], parmC1[2], parmC1[3],
                                                parmC1[4], parmC1[5]))}     }        }
  maxA <- ifelse( ((distA %in% c("dchisq", "df")) && parmA1[1]==1),
                    min(ylimchisq1, max(c(fA(x),fA(x0)),na.rm=TRUE)), max(c(fA(x),fA(x0)))  )
  ylim1 <- c(min(c(fA(x),fA(x0)),na.rm=TRUE), maxA)
  if (!is.null(distB) & is.null(distC))  {
     maxB <- ifelse( ((distB %in% c("dchisq", "df")) && parmB1[1]==1),
                       min(ylimchisq1, max(c(fB(x),fB(x0)),na.rm=TRUE)), max(c(fB(x),fB(x0)))  )
     ylim1 <- c( min(c(fA(x),fB(x),fA(x0),fB(x0)),na.rm=TRUE),  max(maxA, maxB) )   }
  if (is.null(distB) & !is.null(distC))  {
     maxC <- ifelse( ((distC %in% c("dchisq", "df")) && parmC1[1]==1),
                       min(ylimchisq1, max(c(fC(x),fC(x0)),na.rm=TRUE)), max(c(fC(x),fC(x0)))  )
     ylim1 <- c( min(c(fA(x),fC(x),fA(x0),fC(x0)),na.rm=TRUE),  max(maxA, maxC) )   }
  if (!is.null(distB) & !is.null(distC))  {
     maxB <- ifelse( ((distB %in% c("dchisq", "df")) && parmB1[1]==1),
                       min(ylimchisq1, max(c(fB(x),fB(x0)),na.rm=TRUE)), max(c(fB(x),fB(x0)))  )
     maxC <- ifelse( ((distC %in% c("dchisq", "df")) && parmC1[1]==1),
                       min(ylimchisq1, max(c(fC(x),fC(x0)),na.rm=TRUE)), max(c(fC(x),fC(x0)))  )
     ylim1 <- c( min(c(fA(x),fB(x),fC(x),fA(x0),fB(x0),fC(x0)),na.rm=TRUE), max(maxA, maxB, maxC) )   }
  ylab="";   pdfA <- substring(distA, 1, 1)=="d"
  pdfB <- ifelse(is.null(distB), T, (substring(distB, 1, 1)=="d"))
  pdfC <- ifelse(is.null(distC), T, (substring(distC, 1, 1)=="d"))
  if (pdfA & is.null(distB) & is.null(distC)) {ylab <- "probability density function"}
  else { if (pdfA & pdfB & pdfB) ylab <- "probability density functions" }
  cdfA <- substring(distA, 1, 1)=="p"
  cdfB <- ifelse(is.null(distB), T, (substring(distB, 1, 1)=="p"))
  cdfC <- ifelse(is.null(distC), T, (substring(distB, 1, 1)=="p"))
  if (cdfA & is.null(distB) & is.null(distC)) {ylab <- "cumulative distribution function"}
  else { if (cdfA & cdfB & cdfB) ylab <- "cumulative distribution functions" }
  if (discrete.pdfA) { plot(x0, fA(x0), ylim=ylim1, type="h", xlab=xlab, ylab=ylab, col=col[1], lwd=lwd, ...)  }
  else { plot(x, fA(x), ylim=ylim1, type="p", pch=20, cex=0.5, xlab=xlab, ylab=ylab, col=col[1], ...) }
  if (discrete.pdfB) { curve(fB, x0[1], x0[length(x0)], length(x0),
                             add=TRUE, type="h", col=col[2], lwd=lwd)  }
  if (!is.null(distB) & !discrete.pdfB)
            { curve(fB, xmin, xmax, num.grid.points, add=TRUE, type="p",
              pch=20, cex=0.5, col=col[2]) }
  if (discrete.pdfC) { curve(fC, x0[1], x0[length(x0)], length(x0),
                       add=TRUE, type="h", col=col[3], lwd=lwd)  }
  if (!is.null(distC) & !discrete.pdfC)
            { curve(fC, xmin, xmax, num.grid.points, add=TRUE, type="p",
              pch=20, cex=0.5, col=col[3]) }
  if (discrete.pdfA & discrete.pdfB & !discrete.pdfC)   {
     x1 <- x0[ fA(x0) > fB(x0) ];     x2 <- x0[ fA(x0) <= fB(x0) ]
     curve(fA, x1[1], x1[length(x1)], length(x1), add=TRUE, type="h", col=col[1])
     curve(fB, x1[1], x1[length(x1)], length(x1), add=TRUE, type="h", col=col[2])
     curve(fB, x2[1], x2[length(x2)], length(x2), add=TRUE, type="h", col=col[2])
     curve(fA, x2[1], x2[length(x2)], length(x2), add=TRUE, type="h", col=col[1])   }
  if (discrete.pdfA & !discrete.pdfB & discrete.pdfC)   {
     x1 <- x0[ fA(x0) > fC(x0) ];     x2 <- x0[ fA(x0) <= fC(x0) ]
     curve(fA, x1[1], x1[length(x1)], length(x1), add=TRUE, type="h", col=col[1])
     curve(fC, x1[1], x1[length(x1)], length(x1), add=TRUE, type="h", col=col[3])
     curve(fC, x2[1], x2[length(x2)], length(x2), add=TRUE, type="h", col=col[3])
     curve(fA, x2[1], x2[length(x2)], length(x2), add=TRUE, type="h", col=col[1])   }
  if (!discrete.pdfA & discrete.pdfB & discrete.pdfC)   {
     x1 <- x0[ fB(x0) > fC(x0) ];     x2 <- x0[ fB(x0) <= fC(x0) ]
     curve(fB, x1[1], x1[length(x1)], length(x1), add=TRUE, type="h", col=col[2])
     curve(fC, x1[1], x1[length(x1)], length(x1), add=TRUE, type="h", col=col[3])
     curve(fC, x2[1], x2[length(x2)], length(x2), add=TRUE, type="h", col=col[3])
     curve(fB, x2[1], x2[length(x2)], length(x2), add=TRUE, type="h", col=col[2])   }
  if (discrete.pdfA & discrete.pdfB & discrete.pdfC)   {
     x1 <- x0[ fA(x0) >= fB(x0) & fB(x0) >= fC(x0) ]
     x2 <- x0[ fA(x0) >= fC(x0) & fC(x0) >= fB(x0) ]
     x3 <- x0[ fB(x0) >= fA(x0) & fA(x0) >= fC(x0) ]
     x4 <- x0[ fB(x0) >= fC(x0) & fC(x0) >= fA(x0) ]
     x5 <- x0[ fC(x0) >= fA(x0) & fA(x0) >= fB(x0) ]
     x6 <- x0[ fC(x0) >= fB(x0) & fB(x0) >= fA(x0) ]
     curve(fA, x1[1], x1[length(x1)], length(x1), add=TRUE, type="h", col=col[1])
     curve(fB, x1[1], x1[length(x1)], length(x1), add=TRUE, type="h", col=col[2])
     curve(fC, x1[1], x1[length(x1)], length(x1), add=TRUE, type="h", col=col[3])
     curve(fA, x2[1], x2[length(x2)], length(x2), add=TRUE, type="h", col=col[1])
     curve(fC, x2[1], x2[length(x2)], length(x2), add=TRUE, type="h", col=col[3])
     curve(fB, x2[1], x2[length(x2)], length(x2), add=TRUE, type="h", col=col[2])
     curve(fB, x3[1], x3[length(x3)], length(x3), add=TRUE, type="h", col=col[2])
     curve(fA, x3[1], x3[length(x3)], length(x3), add=TRUE, type="h", col=col[1])
     curve(fC, x3[1], x3[length(x3)], length(x3), add=TRUE, type="h", col=col[3])
     curve(fB, x4[1], x4[length(x4)], length(x4), add=TRUE, type="h", col=col[2])
     curve(fC, x4[1], x4[length(x4)], length(x4), add=TRUE, type="h", col=col[3])
     curve(fA, x4[1], x4[length(x4)], length(x4), add=TRUE, type="h", col=col[1])
     curve(fC, x5[1], x5[length(x5)], length(x5), add=TRUE, type="h", col=col[3])
     curve(fA, x5[1], x5[length(x5)], length(x5), add=TRUE, type="h", col=col[1])
     curve(fB, x5[1], x5[length(x5)], length(x5), add=TRUE, type="h", col=col[2])
     curve(fC, x6[1], x6[length(x6)], length(x6), add=TRUE, type="h", col=col[3])
     curve(fB, x6[1], x6[length(x6)], length(x6), add=TRUE, type="h", col=col[2])
     curve(fA, x6[1], x6[length(x6)], length(x6), add=TRUE, type="h", col=col[1])   }
  if (!is.null(distC) & !discrete.pdfC)
     { curve(fC, xmin, xmax, num.grid.points, add=TRUE, type="p",
                  pch=20, cex=0.5, col=col[3]) 
       if (length(additional.x.range)>=2)
          curve(fC, min(additional.x.range), max(additional.x.range), num.grid.points, 
                add=TRUE, type="p", pch=20, cex=0.5, col=col[3])  }
  if (!is.null(distB) & !discrete.pdfB)
     { curve(fB, xmin, xmax, num.grid.points, add=TRUE, type="p",
                  pch=20, cex=0.5, col=col[2]) 
       if (length(additional.x.range)>=2)
          curve(fB, min(additional.x.range), max(additional.x.range), num.grid.points, 
                add=TRUE, type="p", pch=20, cex=0.5, col=col[2])  }
  if (!is.null(distA) & !discrete.pdfA)
     { curve(fA, xmin, xmax, num.grid.points, add=TRUE, type="p",
                  pch=20, cex=0.5, col=col[1]) 
       if (length(additional.x.range)>=2)
          curve(fA, min(additional.x.range), max(additional.x.range), num.grid.points, 
                add=TRUE, type="p", pch=20, cex=0.5, col=col[1])  }
} 


getMinMax <- function(xmin=NULL, xmax=NULL, distA, parmA1=NULL, parmA2=NULL,
             distB=NULL, parmB1=NULL, parmB2=NULL, distC=NULL, parmC1=NULL, parmC2=NULL) {
  # Outputs reasonable minimum and maximum boundaries for plotting functions.
  # Also, outputs the medians.
  # `xmin' and `xmax' represent the limiting values on the x-axis.
  #    Reasonable values of `xmin' and `xmax' are attempted when originally set to NULL.
  # When plotting the pdf, some choices for `distA', `distB' and `distC' are
  #    dnorm, dexp, dcauchy, dlaplace, dt, dt.nc, dchisq, df, df.nc, dunif, dbinom, 
  #    dgeom, dhyper, dpois, dnbinom.
  # When plotting the cdf, some choices for `distA', `distB' and `distC' are
  #    pnorm, pexp, pcauchy, plaplace, pt, pchisq, pf, punif, pbinom, 
  #    pgeom, phyper, ppois, pnbinom.
  # The values of `distA', `distB' and `distC' must be in quotes if `xmin' or `xmax' is NULL.
  # `parmA1' and `parmA2' are numeric scalars or vectors of the parameters of `distA' 
  #    (excluding the first argument), where NULL is the default value; 
  #    likewise for the other two functions.
  # If only one function is to be plotted, then use `distA'.
  # If a distribution, say `distA', needs more than 2 parameters, then assign `parmA1'
  #    to be a vector of all (as many as 5) parameters, and keep `parmA2' equal to NULL.
  #    For example, to plot a noncentral F pdf with 7 and 9 df and ncp=2.5, type:
  #    > plot2.dist( 0, 8, "df.nc", c(7,9,2.5) )
  if (!is.character(distA))
     stop("'distA' must be a probability density function or cumulative distribution function in character format.")
  if (!(substring(distA,1,1) %in% c("d","p"))) stop("'distA' must be a probability density function or cumulative distribution in character format.")
  if (!is.function(try(match.fun(distA),silent=TRUE)))
     stop("'distA' must be a probability density function or cumulative distribution in character format.")
  if (!is.null(distB)) {
     if (!is.character(distB)) 
        stop("'distB' must be a probability density function or cumulative distribution function in character format or NULL.")
     if (!(substring(distB,1,1) %in% c("d","p"))) stop("'distB' must be a probability density function or cumulative distribution in character format or NULL.")
     if (!is.function(try(match.fun(distB),silent=TRUE)))
        stop("'distB' must be a probability density function or cumulative distribution in character format or NULL.")  }
  if (!is.null(distC)) {
     if (!is.character(distC)) 
        stop("'distC' must be a probability density function or cumulative distribution function in character format or NULL.")
     if (!(substring(distC,1,1) %in% c("d","p"))) stop("'distC' must be a probability density function or cumulative distribution in character format or NULL.")
     if (!is.function(try(match.fun(distC),silent=TRUE)))
        stop("'distC' must be a probability density function or cumulative distribution in character format or NULL.")  }
  if (!is.numeric(parmA1) & !is.null(parmA1))  stop("'parmA1' must be numeric or NULL.")
  if (!is.numeric(parmA2) & !is.null(parmA2))  stop("'parmA2' must be numeric or NULL.")
  if (!is.numeric(parmB1) & !is.null(parmB1))  stop("'parmB1' must be numeric or NULL.")
  if (!is.numeric(parmB2) & !is.null(parmB2))  stop("'parmB2' must be numeric or NULL.")
  if (!is.numeric(parmC1) & !is.null(parmC1))  stop("'parmC1' must be numeric or NULL.")
  if (!is.numeric(parmC2) & !is.null(parmC2))  stop("'parmC2' must be numeric or NULL.")
  if (!is.numeric(xmin) & !is.null(xmin))  stop("'xmin' must be numeric or NULL.")
  if (!is.numeric(xmax) & !is.null(xmax))  stop("'xmax' must be numeric or NULL.")
  pmin <- 0.125 ; amountover <- 1 ; buffer <- 0.15
  return.now <- FALSE ;  medianA <- medianB <- medianC <- NULL
  parmA1 <- c(parmA1, parmA2); parmB1 <- c(parmB1, parmB2); parmC1 <- c(parmC1, parmC2)
  if (!is.null(xmin) & !is.null(xmax)) medianA <- medianB <- medianC <- mean(xmin,xmax)
  if ("dt.nc" %in% c(distA, distB, distC)) return.now <- TRUE
  if ("df.nc" %in% c(distA, distB, distC)) return.now <- TRUE
  # if (distA=="pt" & length(parmA1)>1) return.now <- TRUE
  # if (!is.null(distB)) {if (distB=="pt" & length(parmB1)>1) return.now <- TRUE}
  # if (!is.null(distC)) {if (distC=="pt" & length(parmC1)>1) return.now <- TRUE}
  # if (distA=="pf" & length(parmA1)>2) return.now <- TRUE
  # if (!is.null(distB)) {if (distB=="pf" & length(parmB1)>2) return.now <- TRUE}
  # if (!is.null(distC)) {if (distC=="pf" & length(parmC1)>2) return.now <- TRUE}
  if (return.now & (is.null(xmin) || is.null(xmax))) 
    stop("Neither `xmin' nor `xmax' can be NULL for this distribution.")
  if (!return.now) {
    if (!is.null(distA))   qdistA <- paste("q",substring(distA,2),sep="")
    if (!is.null(distB))   qdistB <- paste("q",substring(distB,2),sep="")
    if (!is.null(distC))   qdistC <- paste("q",substring(distC,2),sep="")    
    if (!is.null(distA)) {
       qfA <- function(z)   {
                  if (length(parmA1)==0)
                     {return(match.fun(qdistA)(z))}
                  if (length(parmA1)==1)
                     {return(match.fun(qdistA)(z, parmA1[1]))}
                  if (length(parmA1)==2)
                     {return(match.fun(qdistA)(z, parmA1[1], parmA1[2]))}
                  if (length(parmA1)==3)
                     {return(match.fun(qdistA)(z, parmA1[1], parmA1[2], parmA1[3]))}
                  if (length(parmA1)==4)
                     {return(match.fun(qdistA)(z, parmA1[1], parmA1[2], parmA1[3], parmA1[4]))}
                  if (length(parmA1)==5)
                     {return(match.fun(qdistA)(z, parmA1[1], parmA1[2], parmA1[3],
                                                 parmA1[4], parmA1[5]))}     }     }
    if (!is.null(distB)) {
       qfB <- function(z)   {
                  if (length(parmB1)==0)
                     {return(match.fun(qdistB)(z))}
                  if (length(parmB1)==1)
                     {return(match.fun(qdistB)(z, parmB1[1]))}
                  if (length(parmB1)==2)
                     {return(match.fun(qdistB)(z, parmB1[1], parmB1[2]))}
                  if (length(parmB1)==3)
                     {return(match.fun(qdistB)(z, parmB1[1], parmB1[2], parmB1[3]))}
                  if (length(parmB1)==4)
                     {return(match.fun(qdistB)(z, parmB1[1], parmB1[2], parmB1[3], parmB1[4]))}
                  if (length(parmB1)==5)
                     {return(match.fun(qdistB)(z, parmB1[1], parmB1[2], parmB1[3],
                                                  parmB1[4], parmB1[5]))}     }    }
    if (!is.null(distC)) {
       qfC <- function(z)   {
                  if (length(parmC1)==0)
                     {return(match.fun(qdistC)(z))}
                  if (length(parmC1)==1)
                     {return(match.fun(qdistC)(z, parmC1[1]))}
                  if (length(parmC1)==2)
                     {return(match.fun(qdistC)(z, parmC1[1], parmC1[2]))}
                  if (length(parmC1)==3)
                     {return(match.fun(qdistC)(z, parmC1[1], parmC1[2], parmC1[3]))}
                  if (length(parmC1)==4)
                     {return(match.fun(qdistC)(z, parmC1[1], parmC1[2], parmC1[3], parmC1[4]))}
                  if (length(parmC1)==5)
                     {return(match.fun(qdistC)(z, parmC1[1], parmC1[2], parmC1[3],
                                                  parmC1[4], parmC1[5]))}     }    }
    xminA.temp <- NULL; xmaxA.temp <- NULL;  xminB.temp <- NULL; xmaxB.temp <- NULL
    xminC.temp <- NULL; xmaxC.temp <- NULL
    xmin2A.temp <- NULL; xmax2A.temp <- NULL;  xmin2B.temp <- NULL; xmax2B.temp <- NULL
    xmin2C.temp <- NULL; xmax2C.temp <- NULL
    if (is.null(xmin)) {
      if (!is.null(distA))   xminA.temp <- qfA(pmin)
      if (!is.null(distB))   xminB.temp <- qfB(pmin)
      if (!is.null(distC))   xminC.temp <- qfC(pmin)    }
    else  { xminA.temp <- xmin; xminB.temp <- xmin; xminC.temp <- xmin }
    if (is.null(xmax)) {
      if (!is.null(distA))   xmaxA.temp <- qfA(1-pmin)
      if (!is.null(distB))   xmaxB.temp <- qfB(1-pmin)
      if (!is.null(distC))   xmaxC.temp <- qfC(1-pmin)  }
    else  { xmaxA.temp <- xmax; xmaxB.temp <- xmax; xmaxC.temp <- xmax }
    if (!is.null(distA)) {
      xmin2A.temp <- xminA.temp - (xmaxA.temp-xminA.temp)*amountover
      xmax2A.temp <- xmaxA.temp + (xmaxA.temp-xminA.temp)*amountover
      qfA0 <- ifelse( is.na(qfA(0)), -Inf, qfA(0) )
      qfA1 <- ifelse( is.na(qfA(1)),  Inf, qfA(1) )
      left.bounded  <- ( xmin2A.temp < qfA0 )
      right.bounded <- ( xmax2A.temp > qfA1 )
      xmin2A.temp <- max( xmin2A.temp, qfA0 )
      xmax2A.temp <- min( xmax2A.temp, qfA1 )
      if (left.bounded)  xmin3A.temp <- xmin2A.temp - (xmax2A.temp - xmin2A.temp) * buffer
      if (right.bounded) xmax3A.temp <- xmax2A.temp + (xmax2A.temp - xmin2A.temp) * buffer
      if (left.bounded)  xmin2A.temp <- xmin3A.temp
      if (right.bounded) xmax2A.temp <- xmax3A.temp     }
    if (!is.null(distB)) {
      xmin2B.temp <- xminB.temp - (xmaxB.temp-xminB.temp)*amountover
      xmax2B.temp <- xmaxB.temp + (xmaxB.temp-xminB.temp)*amountover
      qfB0 <- ifelse( is.na(qfB(0)), -Inf, qfB(0) )
      qfB1 <- ifelse( is.na(qfB(1)),  Inf, qfB(1) )
      left.bounded  <- ( xmin2B.temp < qfB0 )
      right.bounded <- ( xmax2B.temp > qfB1 )
      xmin2B.temp <- max( xmin2B.temp, qfB0 )
      xmax2B.temp <- min( xmax2B.temp, qfB1 )
      if (left.bounded)  xmin3B.temp <- xmin2B.temp - (xmax2B.temp - xmin2B.temp) * buffer
      if (right.bounded) xmax3B.temp <- xmax2B.temp + (xmax2B.temp - xmin2B.temp) * buffer
      if (left.bounded)  xmin2B.temp <- xmin3B.temp
      if (right.bounded) xmax2B.temp <- xmax3B.temp     }
    if (!is.null(distC)) {
      xmin2C.temp <- xminC.temp - (xmaxC.temp-xminC.temp)*amountover
      xmax2C.temp <- xmaxC.temp + (xmaxC.temp-xminC.temp)*amountover
      qfC0 <- ifelse( is.na(qfC(0)), -Inf, qfC(0) )
      qfC1 <- ifelse( is.na(qfC(1)),  Inf, qfC(1) )
      left.bounded  <- ( xmin2C.temp < qfC0 )
      right.bounded <- ( xmax2C.temp > qfC1 )
      xmin2C.temp <- max( xmin2C.temp, qfC0 )
      xmax2C.temp <- min( xmax2C.temp, qfC1 )
      if (left.bounded)  xmin3C.temp <- xmin2C.temp - (xmax2C.temp - xmin2C.temp) * buffer
      if (right.bounded) xmax3C.temp <- xmax2C.temp + (xmax2C.temp - xmin2C.temp) * buffer
      if (left.bounded)  xmin2C.temp <- xmin3C.temp
      if (right.bounded) xmax2C.temp <- xmax3C.temp     }
    if (is.null(xmin))   xmin <- min(xmin2A.temp, xmin2B.temp, xmin2C.temp)
    if (is.null(xmax))   xmax <- max(xmax2A.temp, xmax2B.temp, xmax2C.temp)
    if (xmin>xmax) {temp <- xmin; xmin <- xmax; xmax <- temp}
    medianA <- ifelse( is.null(distA), NA, qfA(0.5) )
    medianB <- ifelse( is.null(distB), NA, qfB(0.5) )
    medianC <- ifelse( is.null(distC), NA, qfC(0.5) )
  }  
  structure(list(xmin=xmin, xmax=xmax, medianA=medianA, medianB=medianB, medianC=medianC))
} 


plotLine <- function(x, y=NULL, data=NULL, xlab=NULL, ylab=NULL, pch=19, col=c("black", "red"),
                       digits.intercept=NULL, digits.slope=NULL, ...) {
  # Constructs a scatterplot with the fitted regression line.
  # `x' is the vector of independent observations.
  # `y' is the vector of dependent observations.
  # `xlab' is the label for the x-axis.
  # `ylab' is the label for the y-axis.
  # `pch' is the plotting character, usually 19 to 25; 21 is a hollow circle.
  #    Type `?points' for more details.
  # `col' may be a scalar or vector of length 2, and specifies the color.
  #    col[1] is the color of the points, and col[2] is the color of the plotted line.
  # `digits.intercept' is the desired number of digits after the decimal point for the intercept.
  # `digits.slope' is the desired number of digits after the decimal point for the slope.
  # ...: optional arguments to `plot'.
  if (!is.numeric(digits.intercept) & !is.null(digits.intercept))  stop("'digits.intercept' must be numeric or NULL.")
  if (!is.numeric(digits.slope) & !is.null(digits.slope))  stop("'digits.slope' must be numeric or NULL.")
  if (length(col)==1) col <- rep(col,2)
  if  (is.null(y))   form = x
  if (!is.null(y))   form = y~x
  if (is.null(xlab))  xlab = colnames(model.frame(lm(form, data=data)))[2]
  if (is.null(ylab))  ylab = colnames(model.frame(lm(form, data=data)))[1]
  z <- lm(form, data=data)
  intercept <- as.numeric(z$coefficients[1]); slope <- as.numeric(z$coefficients[2])
  main1 <- c("Simple Linear Regression Line",  paste( ylab, "=",
              format(intercept, digits=digits.intercept),
              ifelse(slope>=0,"+","-"),
              format(abs(slope), digits=digits.slope), xlab )  )
  ylim <- c( min(c(model.frame(z)[,1],fitted(z))), max(c(model.frame(z)[,1],fitted(z))) )
  plot(form, ylim=ylim, data=data,  type="p", xlab=xlab, ylab=ylab, main=main1, col=col[1], pch=pch, ...)
  abline(a=intercept, b=slope, col=col[2], lwd=2)
}

