#' Goodness of fit tests for distance sampling models
#'
#' Generic function that computes chi-square goodness of fit test for detection function models with binned data and Cramer-von Mises and Kolmogorov-Smirnov tests for exact distance data. By default a Q-Q plot is generated for exact data (and can be surpressed using the \code{qq=FALSE} argument.
#'
#' @aliases ddf.gof gof.io gof.io.fi gof.trial gof.trial.fi gof.rem
#'   gof.rem.fi
#' @export
#' @param model model object
#' @param breaks Cutpoints to use for binning data
#' @param nc Number of distance classes
#' @param qq Flag to indicate whether quantile-quantile plot is desired
#' @param \dots Graphics parameters to pass into qqplot function
#' @return List of class \code{ddf.gof} containing \item{chi-square }{Goodness of fit test statistic} \item{df}{Degrees of freedom associated with test statistic} \item{p-value }{Significance level of test statistic}
#' @author Jeff Laake
#' @seealso \code{\link{qqplot.ddf}}
#' @keywords utility
ddf.gof <- function(model,breaks=NULL,nc=NULL,qq=TRUE,...){
  # Functions Used: gof.ds, gof.io, gof.io.fi, gof.trial,
  #                 gof.trial.fi, qqplot.df

  if(!is.null(breaks)){
    breaks <- test.breaks(breaks,model$meta.data$left,model$meta.data$width)
    nc <- length(breaks)-1
  }else if(model$meta.data$binned){
    breaks <- model$meta.data$breaks
  }

  # call method specific function
  result <- switch(model$method,
                   ds       = gof.ds(model,breaks,nc),
                   io       = gof.io(model,breaks,nc),
                   io.fi    = gof.io.fi(model,breaks,nc),
                   trial    = gof.trial(model,breaks,nc),
                   trial.fi = gof.trial.fi(model,breaks,nc),
                   rem      = gof.rem(model,breaks,nc),
                   rem.fi   = gof.rem.fi(model,breaks,nc))

  if(!model$meta.data$binned){
    result <- list(chisquare=result,dsgof=qqplot.ddf(model, plot=qq, ...))
  }else{
    result <- list(chisquare=result)
  }

  class(result)=c("ddf.gof")
  return(result)
}
