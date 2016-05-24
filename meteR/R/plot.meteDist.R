#' @title Plot METE distributions and associated data
#'
#' @description
#' \code{plot.meteDist} plots both the theoretical prediction and data for a 
#' \code{meteDist} object using either a rank or cumulative distribution plot
#'
#' @details
#' \code{plot.meteDist} automatically extracts the prediction and data (if used 
#' in \code{meteESF}) from the \code{meteDist} object. Additional plotting 
#' arguments can be passed to \code{...}.
#' 
#' @param x a \code{meteDist} object
#' @param ptype type of plot; either "cdf" or "rad"
#' @param th.col line color of theoretical prediction 
#' @param lower.tail logical; choose TRUE to highlight differences between data and theory at low abundance; choose FALSE to highlight differences at high abundance.
#' @param add.legend logical; add a legend
#' @param ... arguments to be passed to \code{plot}
#' @param add.line add the curve for a fitted model to the existing plot
# @keywords manip
#' @export
#' @importFrom graphics curve legend plot points
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                abund=arth$count,
#'                power=arth$mass^(.75),
#'                minE=min(arth$mass^(.75)))
#' ipd1 <- ipd(esf1)
#' plot(ipd1)
#' plot(ipd1, ptype='rad')
#' 
# @return list
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso sad, ipd, ssad, sipd, print.meteDist
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.


plot.meteDist <- function(x, ptype=c("cdf","rad"), th.col="red", 
                          lower.tail=TRUE, add.legend=TRUE, add.line=FALSE, ...) {
  ptype <- match.arg(ptype,c("cdf", "rad"))
  
  plot.par <- list(...)
  
  if(!('ylab' %in% names(plot.par))) {
    ylab <- ifelse(ptype=='cdf', 'Cumulative probability', '%s')
    ylab <- sprintf(ylab, switch(x$type,
                                 'sad' = 'Abundance',
                                 'ipd' = 'Metabolic rate',
                                 'sipd' = 'Metabolic rate'))
    plot.par$ylab <- ylab
  }
  
  if(!('xlab' %in% names(plot.par))) {
    xlab <- ifelse(ptype=='cdf', '%s', 'Rank')
    xlab <- sprintf(xlab, switch(x$type,
                                 'sad' = 'Abundance',
                                 'ipd' = 'Metabolic rate',
                                 'sipd' = 'Metabolic rate'))
    plot.par$xlab <- xlab
  }
  
  if(ptype=="cdf") {
    this.curve <- x$p
    ## if no data, don't plot it, just plot the curve
    if(!is.null(x$data)) {
      xmax <- max(x$data)
      X <- .ecdf(x$data, !lower.tail)
    } else {
      xmax <- ifelse(is.finite(max(plot.par$xlim)), 
                     max(plot.par$xlim), 
                     x$state.var['N0']/x$state.var['S0'])
      X <- cbind(c(1, floor(xmax)), this.curve(c(1, floor(xmax))))
      plot.par$type <- 'n'
    }
    
    do.call(plot, append(list(x=X),plot.par))
    
    if(x$type %in% c("gsd", "sad")) {
      this.supp <- 1:xmax
      points(this.supp, this.curve(this.supp,lower.tail=lower.tail),
             type="l", col=th.col)
    } else {
      curve(this.curve(x,lower.tail=lower.tail), add=TRUE, col=th.col)
    }
  } else {
    ## if no data, don't plot it, just plot the rank fun
    if(is.null(x$data)) {
      X <- meteDist2Rank(x)
      plot.par$type <- 'n'
    } else {
      X <- x$data
    }
    
    ## if ylim not already specified make sure both data and theory fit
    if(!('ylim' %in% names(plot.par))) {
      plot.par$ylim <- range(X, x$rankFun)
    }
    
    ## do plotting
    do.call(plot, append(list(x=X),plot.par))
    points(meteDist2Rank(x), type="l", col=th.col)
  }
  
  if(add.legend) legend(ifelse(ptype=='cdf', 'bottomright', 'topright'), 
                        legend=c('data', 'METE'), col=c('black', 'red'),
                        lty=c(NA, 1), pch=c(21, NA), bty='n') 
}

