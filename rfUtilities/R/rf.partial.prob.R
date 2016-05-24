#' @title Random Forest probability scaled partial dependency plots
#' @description Produces partial dependency plots with probability distribution based on scaled margin distances.
#'       
#' @param x Object of class randomForest 
#' @param pred.data Training data.frame used for constructing the plot, 
#' @param xname Name of the variable for calculating partial dependence 
#' @param which.class The class to focus on
#' @param w Weights to be used in averaging (if not supplied, mean is not weighted) 
#' @param prob Scale distances to probabilities
#' @param plot Plot results (TRUE/FALSE)
#' @param smooth Apply spline.smooth to y
#' @param raw Plot unsmoothed values
#' @param rug Draw hash marks on plot representing deciles of x
#' @param n.pt Number of points on the grid for evaluating partial dependence.
#' @param xlab x-axis plot label
#' @param ylab y-axis plot label
#' @param main Plot label for main
#' @param ... Additional graphical parameters passed to plot
#
#' @return A list class object with fit x,y. If smooth=TRUE y represents smoothed scaled margin distance values 
#'
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references
#' Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species distribution and change using Random Forests CH.8 in Predictive Modeling in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer 
#' @references
#' Baruch-Mordo, S., J.S. Evans, J. Severson, J. D. Naugle, J. Kiesecker, J. Maestas, & M.J. Falkowski (2013) Saving sage-grouse from the trees: A proactive solution to reducing a key threat to a candidate species Biological Conservation 167:233-241  
#'
#' @examples 
#' require(randomForest)
#'   data(iris)
#'     iris.rf <- randomForest(iris[,1:4], iris[,5])
#'       par(mfrow=c(2,2))
#'         for(i in names(iris)[1:4]) {     
#'           rf.partial.prob(iris.rf, iris, i, "setosa", smooth=TRUE, raw=TRUE, rug=FALSE)
#'          }
#'
#' @export
rf.partial.prob <- function(x, pred.data, xname, which.class, w, prob=TRUE, plot=TRUE,
                            smooth=FALSE, raw=FALSE, rug=FALSE, n.pt, xlab, ylab, main, ...) {  
    if (!inherits(x, "randomForest")) stop("x is not randomForest class object")
	if (is.null(x$forest)) stop("Object does not contain an ensemble \n")
	if(!x$type != "regression")	stop("Regression not supported \n")	   
	if (missing(which.class)) stop("Class name missing \n")
	if (missing(xname)) stop("X Variable name missing \n")
	if (missing(x)) stop("randomForest object missing \n")
    if (missing(pred.data)) stop("New data missing \n")
	
	focus <- charmatch(which.class, colnames(x$votes))
    if (is.na(focus)) stop(which.class, "is not one of the class labels")
    xv <- pred.data[, xname]
	n <- nrow(pred.data)
	  if(missing(n.pt)) n.pt <- min(length(unique(pred.data[, xname])), 51)  
	  if (missing(w)) w <- rep(1, n)
	scale.dist <- function(d) { return( (exp(d) - min(exp(d))) / (max(exp(d)) - min(exp(d))) ) }
	
    if(is.factor(xv) && !is.ordered(xv)) {
      x.pt <- levels(xv)
      y.pt <- numeric(length(x.pt))
        for (i in seq(along = x.pt)) {
          x.data <- pred.data
          x.data[, xname] <- factor(rep(x.pt[i], n), levels = x.pt)
          pr <- stats::predict(x, x.data, type = "prob")
          y.pt[i] <- stats::weighted.mean(log(ifelse(pr[, focus] > 0, pr[, focus], .Machine$double.eps)) -
                       rowMeans(log(ifelse(pr > 0, pr, .Machine$double.eps))), w, na.rm=TRUE)
      }  
      if( prob == TRUE) { y.pt <- scale.dist(y.pt) } 
      if (plot) {
        if(missing(xlab)) xlab=xname
	    if(missing(ylab)) ylab=which.class
	    if(missing(main)) main="Partial Dependency Plot"	
	      graphics::barplot(y.pt, width=rep(1,length(y.pt)), col="blue",
                            xlab=xlab, ylab=ylab, main=main,
                            names.arg=x.pt, ...)
      }				  
    }	
	
    if (is.numeric(xv)) {
      x.pt <- seq(min(xv), max(xv), length=n.pt)
      y.pt <- numeric(length(x.pt))
        for (i in seq(along=x.pt)) {
          x.data <- pred.data
          x.data[, xname] <- rep(x.pt[i], n)
          pr <- stats::predict(x, x.data, type = "prob")
          y.pt[i] <- stats::weighted.mean(log(ifelse(pr[, focus] == 0, .Machine$double.eps, 
               pr[, focus])) - rowMeans(log(ifelse(pr == 0, .Machine$double.eps, pr))),
               w, na.rm=TRUE)
        }
    if( prob == TRUE) { y.pt <- scale.dist(y.pt) }       
	  if (plot) {
        if(missing(xlab)) xlab=xname
	      if(missing(ylab)) ylab="probability"
	        if(missing(main)) main=paste("Partial Dependency Plot for Class", which.class, sep=" - ") 			  
    if(smooth)  
      graphics::plot(x.pt, stats::smooth.spline(x.pt, y.pt)$y, type = "l", xlab=xlab, 
		             ylab=ylab, main=main, ...)
	if(raw == TRUE) graphics::lines( y=y.pt, x=x.pt, col="grey", lty=3) 
	  } else {
      graphics::plot(x.pt, y.pt, type = "l", xlab=xlab, ylab=ylab, main=main, ...)		
      }	
    if (rug) {
      if (n.pt > 10) {
        graphics::rug(stats::quantile(xv, seq(0.1, 0.9, by=0.1)), side = 1)
          } else {
        graphics::rug(unique(xv, side = 1))
        }
	  }
    }
  invisible(list(x=x.pt, y=y.pt))
}
