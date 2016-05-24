#' Plot massively parallel semiparametric models
#' 
#' Given a massively parallel smoothing object produced by
#' \code{\link{semipar.mp}}, the function plots the fitted smooth(s) for a
#' given point (e.g., at a given voxel).
#' 
#' 
#' @param x an object of class "\code{\link{semipar.mp}}".
#' @param Y an \eqn{n \times V} outcome matrix.
#' @param arr.ind a 3-element vector specifying the element of the
#' 3-dimensional array of locations (e.g., voxels) for which plotting is
#' desired. If \code{NULL}, \code{which.vox} must be specified.
#' @param which.vox the index of the voxel to be plotted. If \code{NULL},
#' \code{arr.ind} must be specified.
#' @param which.smooth the index of the smooth term of which the confidence
#' interval plot is to be displayed. The default value is \code{NULL} which
#' refers to displaying the plots for all the smooth terms in the model.
#' @param coverage the confidence level of the pointwise confidence intervals
#' in the plot.
#' @param length.new length of the vector of ordered variables with which to
#' predict.
#' @param ylim,ylab, arguments to be passed to \code{\link[graphics]{plot}}.
#' @param ... arguments to be passed to \code{\link[graphics]{plot}}.
#' @author Yin-Hsiu Chen \email{enjoychen0701@@gmail.com}, Philip Reiss
#' \email{phil.reiss@@nyumc.org}, and Lan Huo
#' @examples
#' 
#' n<-32
#' Ys <- matrix(0, n, 5)
#' for(i in 1:n) Ys[i,]<--2:2+rnorm(5, i^2, i^0.5)+sin(i)
#' x1 <- rnorm(n,0,5)
#' x2 <- 1:n+runif(n, 1, 20)
#' semipar.obj <- semipar.mp(~x1+sf(x2,k=10),Y=Ys,lsp=seq(5,50,,30))
#' plot(semipar.obj, Y=Ys, which.vox=2)
#' @export
plot.semipar.mp <-		
function(x, Y, arr.ind = NULL, which.vox = NULL, which.smooth = NULL, coverage = 0.95, length.new = 100, ylim = NULL, ylab = NULL,  ...) {		
	if (is.null(arr.ind) & is.null(which.vox)) stop("Must specify either 'arr.ind' or 'which.vox'")	
		
    tf <- terms.formula(x$formula, specials = "sf")		
    trmstrings <- attr(tf, "term.labels")		
    terms <- rep(NA, length(x$where.sf)); smooth.terms <- terms		
    for (i in 1:length(x$where.sf)) {		
	    formula.term <- sub("sf\\(", "", sub("\\)", "", trmstrings[x$where.sf[i]]))	
        split.term <- as.vector(strsplit(formula.term, ",")[[1]])		
        terms[i] <- split.term[1]		
        smooth.terms[i] <- paste("sf(", split.term[1], ")", sep="")		
	}	
		
    alpha = 1-coverage	
    intercept <- attr(terms.formula(x$formula, specials="sf"), "intercept")		
		
    ################################################################################		
 		
    B.list <- plot.list <- vector("list", length(x$where.sf))		
    for  (i in 1:length(x$where.sf)) {		
	 xarg <- x$list.all[[x$where.sf[i]]]$argvals	
       x.new <- seq(range(xarg)[1], range(xarg)[2], length.out = length.new)		
	 plot.list[[i]]$x <- xarg	
    	 plot.list[[i]]$x.new <- x.new	
       B.list[[i]] = eval.basis(plot.list[[i]]$x.new, x$list.all[[x$where.sf[i]]]$basis)           		
    }		
     pred.ind <- rep(NA, length(x$where.sf)+length(x$where.nsf))		
     for(i in (x$where.sf)) pred.ind[i] <- x$list.all[[i]]$k		
     for(i in (x$where.nsf)) pred.ind[i] <- ncol(x$list.all[[i]]$modmat)		
     start.ind <- c(1, cumsum(pred.ind)[-length(pred.ind)]+1) + intercept		
     end.ind <- start.ind + pred.ind - 1   		
     modmat.tmp <- matrix(NA, length.new, NCOL(x$modmat))		
     if(intercept) modmat.tmp[, 1] <- 1		
     if(length(x$where.nsf)>0) {		
     for(i in 1:length(x$where.nsf)) {		
    		modmat.tmp[, (start.ind[x$where.nsf[i]]):(end.ind[x$where.nsf[i]])] <- matrix(apply(as.matrix(x$modmat[, (start.ind[x$where.nsf[i]]):(end.ind[x$where.nsf[i]])]), 2, mean), length.new, byrow = TRUE)
     }}		
     for(i in 1:length(x$where.sf)) {		
    		modmat.tmp[, (start.ind[x$where.sf[i]]):(end.ind[x$where.sf[i]])] <- 0
     }		
     Xarray <- array(modmat.tmp, dim = c(dim(modmat.tmp), length(x$where.sf))) 		
     for(i in 1:length(x$where.sf)) {		
        Xarray[,(start.ind[x$where.sf[i]]):(end.ind[x$where.sf[i]]),i] <- B.list[[i]]		
    } 		
    ################################################################################		
       coef <- x$coef[, which.vox]		
          for  (i in 1:length(x$where.sf)) {		
		coef.seq <- (start.ind[x$where.sf[i]]):(end.ind[x$where.sf[i]])
            plot.list[[i]]$y.hat <- Xarray[,coef.seq,i]%*%coef[coef.seq]		
          }		
     if (is.null(which.smooth)) plot.ind <- 1:length(x$where.sf)		
	else plot.ind <- which.smooth	
        par(mfrow=c(1, length(plot.ind)))		
        for (i in plot.ind) {    		
		B <- Xarray[,,i]%*%x$ttu      
		pwvar <- x$sigma2[which.vox] * (B %*% x$RinvU)^2 %*% (1 / (1 + exp(x$pwlsp[which.vox]) * x$tau))
		
		ord <- order(plot.list[[i]]$x.new)
		mul = qnorm(1-alpha/2)	
            lower <- min(plot.list[[i]]$y.hat-mul*sqrt(pwvar))		
            upper <- max(plot.list[[i]]$y.hat+mul*sqrt(pwvar))		
            range <- upper-lower	
            if (is.null(ylim)) ylim = c(lower-range*0.1, upper+range*0.1)
            if (is.null(ylab)) ylab = smooth.terms[i]
            plot(plot.list[[i]]$x.new[ord], plot.list[[i]]$y.hat[ord], xlab = terms[i], ylab = ylab, type="l", ylim= ylim, ...)		
            rug(plot.list[[i]]$x)		
            lines(plot.list[[i]]$x.new[ord], plot.list[[i]]$y.hat[ord]-mul*sqrt(pwvar[ord]), lty=2)		
            lines(plot.list[[i]]$x.new[ord], plot.list[[i]]$y.hat[ord]+mul*sqrt(pwvar[ord]), lty=2)    		
        } 		
}		
		
