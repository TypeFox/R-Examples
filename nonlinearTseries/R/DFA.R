################################################################################
#' Detrended Fluctuation Analysis
#' @description
#' Functions for performing Detrended Fluctuation Analysis (DFA), a widely used
#' technique for detecting long range correlations in time series. These functions
#' are able to estimate several scaling exponents from the time series being analyzed. 
#' These scaling exponents  characterize short or long-term fluctuations, depending of
#' the range used for regression (see  details).
#' @details The Detrended Fluctuation Analysis (DFA) has become a widely used
#' technique for detecting long range correlations in time series. The DFA procedure
#' may be summarized as follows:
#' \enumerate{
#'  \item Integrate the time series to be analyzed. The time series resulting from the
#'  integration will be referred to as the profile.
#'  \item Divide the profile into N non-overlapping segments.
#'  \item  Calculate the local trend for each of the segments using least-square regression.
#'  Compute the total error for each ofi the segments.
#'  \item Compute the average of the total error over all segments and take its root square. By repeating 
#'  the previous steps for several segment sizes (let's denote it by t), we obtain the
#'  so-called Fluctuation function \eqn{F(t)}.
#'  \item  If the data presents long-range power law correlations:  \eqn{F(t) \sim t^\alpha}{F(t) proportional t^alpha} and
#'  we may estimate using regression.
#'  \item  Usually, when plotting \eqn{\log(F(t))\;Vs\;log(t)}{log(F(t)) Vs log(t)} we may distinguish two linear regions.
#'  By regression them separately, we obtain two scaling exponents, \emph{\eqn{\alpha_1}{alpha1}}  
#'  (characterizing short-term fluctuations) and \emph{\eqn{\alpha_2}{alpha2}} (characterizing long-term fluctuations). 
#' }
#' Steps 1-4 are performed using the \emph{dfa} function. In order to obtain a estimate 
#' of some scaling exponent, the user must use the  \emph{estimate} function specifying
#' the regression range (window sizes used to detrend the series).
#' @param time.series The original time series to be analyzed.
#' @param npoints The number of different window sizes that will be used to estimate 
#' the Fluctuation function in each zone.
#' @param window.size.range Range of values for the windows size that will be used
#' to estimate the fluctuation function. Default: c(10,300).
#' @param do.plot logical value. If TRUE (default value), a plot of the Fluctuation function is shown.
#' @param ... Additional graphical parameters.
#' @return A \emph{dfa}  object.
#' @author Constantino A. Garcia
#' @examples
#' \dontrun{
#'  white.noise = rnorm(5000)
#' dfa.analysis = dfa(time.series = white.noise, npoints = 10, 
#'                  window.size.range=c(10,1000), do.plot=FALSE)
#' white.estimation = estimate(dfa.analysis,do.plot=TRUE)
#' cat("Theorical: 0.5---Estimated: ",white.estimation ,"\n")
#' 
#' 
#' library(fArma)
#' fgn = as.numeric(fArma::fgnSim(n = 2000, H = 0.75))
#' dfa.analysis = dfa(time.series = fgn, npoints = 30, 
#'                    window.size.range=c(10,1000), 
#'                    do.plot=FALSE)
#' fgn.estimation = estimate(dfa.analysis, do.plot = TRUE,
#'                           fit.col="blue",fit.lwd=2,fit.lty=2,
#'                           main="Fitting DFA to fGn")
#' cat("Theorical: 0.75---Estimated: ",fgn.estimation ,"\n")
#' 
#' 
#' fbm = as.numeric(fArma::fbmSim(n = 2000, H = 0.25))
#' dfa.analysis = dfa(time.series = fbm, npoints = 50, 
#'                    window.size.range=c(10,300), 
#'                    do.plot=FALSE)
#' fbm.estimation = estimate(dfa.analysis,do.plot = TRUE,
#'                           add.legend=F, main="DFA of fBm")
#' cat("Theorical: 1.25 ---Estimated: ",fbm.estimation ,"\n")
#' 
#' }
#' @rdname dfa
#' @export dfa
#' @exportClass dfa
dfa=function(time.series, window.size.range=c(10,300), npoints=20,do.plot=TRUE,
             ...){
  f="y~x"
  window.sizes=c()
  fluctuation.function=c()
  
  if ( (window.size.range[[2]] < window.size.range[[1]]) || 
         any(window.size.range < 1)){
    stop("Invalid window size range")
  }
  # compute windows' sizes... 
  if (npoints >= (window.size.range[[2]] - window.size.range[[1]] + 1)){
    # compute Fn for all possible windows
    window.sizes = window.size.range[[1]]:window.size.range[[2]]
  }else{
    # Try equally spaced windows in log scale
    log.window.sizes = seq(log10(window.size.range[[1]]),
                           log10(window.size.range[[2]]),
                           len=npoints)
    window.sizes = unique(round(10^log.window.sizes))
    # if there are repeated elements we select missing windows until
    # reaching npoints
    if (length(window.sizes) < npoints){ 
      missing = setdiff(window.size.range[[1]]:window.size.range[[2]],
                        window.sizes)
      points.needed = npoints - length(window.sizes) 
      add.points = missing[seq(1,points.needed,len=points.needed)]
      window.sizes = sort(c(window.sizes,add.points))
    }
  }
  npoints = length(window.sizes)
  #integrate the time series time.series: this is called the "profile" in Penzel et al.
  Y=cumsum( time.series - mean(time.series) )
  
  #Calculate the local trend for
  #each of the segments by a least-square fit
  #of the time.series and compute the error
  for (i in 1:npoints){
    fluctuation.function[[i]]=calculateFluctuationFunction(Y,f,window.sizes[[i]])
  }
  
  #create dfa object
  dfa.structure = list(window.sizes=window.sizes,fluctuation.function=fluctuation.function)
  class(dfa.structure) = "dfa"
  # add attributes
  id=deparse(substitute(time.series))
  attr(dfa.structure,"id") = id
  
    
  #plot
  if (do.plot){
    plot(dfa.structure,...)
  }
  
  return(dfa.structure)
}



#' Returns the window sizes used for DFA in a \emph{dfa} object.
#' @param x A \emph{dfa} object.
#' @return The \emph{windowSizes} function returns the windows sizes used
#' to detrend the time series in the DFA. 
#' @seealso \code{\link{dfa}}
#' @export windowSizes
windowSizes = function(x){
  UseMethod("windowSizes")
}

#' @return The \emph{windowSizes} function returns the windows sizes used
#' to detrend the time series. 
#' @rdname dfa
#' @method windowSizes dfa
#' @export
windowSizes.dfa = function(x){
  return (x$window.sizes)
}


#' Returns the fluctuation function obtained in a DFA and represented by a
#' \emph{dfa} object.
#' @param x A \emph{dfa} object.
#' @return The \emph{fluctuationFunction} function returns the fluctuation function used
#' obtained in the DFA. 
#' @seealso \code{\link{dfa}}
#' @export fluctuationFunction
fluctuationFunction = function(x){
  UseMethod("fluctuationFunction")
}


#' @return The \emph{fluctuationFunction} function returns the fluctuation function
#' obtained in the DFA represented by the \emph{dfa} object.
#' @rdname dfa
#' @method fluctuationFunction dfa
#' @export
fluctuationFunction.dfa = function(x){
  return (x$fluctuation.function)
}

#' @rdname dfa
#' @method plot dfa
#' @param main A title for the plot.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param log A character string which contains "x" if the x axis is to be 
#' logarithmic, "y" if the y axis is to be logarithmic and "xy" or "yx" if 
#' both axes are to be logarithmic.
#' @export
#' 
plot.dfa = function(x, main="Detrended fluctuation analysis\n",
                    xlab="Window size: t", ylab="Fluctuation function: F(t)",
                    log="xy",...){
  plot(x$window.sizes, x$fluctuation.function,
       main=main,xlab=xlab,ylab=ylab,log=log,...)
}
#' @rdname dfa
#' @method estimate dfa
#' @param x A \emph{dfa} object.
#' @param regression.range Vector with 2 components denoting the range where the function will perform linear regression.
#' @param fit.col A colors to plot the regression line.
#' @param fit.lty The type of line to plot the regression line.
#' @param fit.lwd The width of the line for the regression line.
#' @param add.legend add a legend with the resulting estmation to the plot?
#' @export
#' 
estimate.dfa = function(x, regression.range = NULL,do.plot=FALSE,
                        fit.col=2,fit.lty=1,fit.lwd=1,add.legend=T,...){
  log.fluctuation.function = log10(x$fluctuation.function)
  log.window.sizes = log10(x$window.sizes)
  if (is.null(regression.range)){ regression.range = range(x$window.sizes)}
  #prepare interpolations to calculate the exponents
  index.alpha = which(x$window.sizes>=regression.range[[1]]&x$window.sizes<=regression.range[[2]])
  #get alpha
  interp1=lm(log.fluctuation.function[index.alpha]~log.window.sizes[index.alpha])
  alpha=interp1$coefficients[[2]]
  if (do.plot){
    plot(x,...)
    lines(x$window.sizes[index.alpha],10^interp1$fitted.values,
          col=fit.col,lty=fit.lty,lwd=fit.lwd)
    if(add.legend){
    legend("bottomright",col=fit.col,lty=fit.lty,lwd=fit.lwd,bty="n",
           legend=paste("Scaling exponent estimate =", sprintf("%.2f",alpha)))
    }
  }
  return(alpha)
}

# private
# Calculate the local trend for
# each of the segments by a least-square fit
# of the time.series and compute the error
# Y: the profile
# fint: interpolation formula
# l: size of the window
calculateFluctuationFunction=function(Y,fint,l){
  #preliminary definitions
  ldata=length(Y)
  ngroups=floor(ldata/l)
  x=1:l
  
  # Divide the profile in ngroups non-overlaping groups. Interpolate and calculate the residual error
  F2=c()
  for (i in 1:ngroups){
    beg=(i-1)*l+1; en=i*l;
    y=Y[beg:en]
    #interpolate
    fit=lm(as.formula(fint))
    #get error
    F2[i]=sum(fit$residuals^2)/l
  }
  return (sqrt(sum(F2)/ngroups))
}

