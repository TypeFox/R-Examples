#' @title Calculate marginal correlations with response, 
#' from a funeigen object
#' @description A function for internal use.  Its main job is to be called by \code{MarginalCor},
#' and do the technical work for calculating estimated marginal correlations. It 
#' uses R. A. Fisher's classic r-to-z transform to create confidence
#' intervals for the correlations.  This process is explained
#' in easy-to-follow detail by David Shen and Zaizai Lu in a technical report.
#' @param object An object of type \code{funeigen}.
#' @param id The vector of subject id's. These tell which responses in \code{response} 
#' correspond to which curves in \code{object}.
#' @param response The vector of responses
#' @param alpha The alpha level for confidence intervals (one minus the two-sided coverage)
#' @note  The confidence intervals are simply based on Fisher's
#'  r-to-z transform and do not 
#'  take into account the uncertainty in estimating the 
#'  smoothed value of x(t). 
#' @references Shen, D., and Lu, Z. (2006). Computation of correlation coefficient
#' and its confidence interval in SAS (r). SUGI 31 (March 26-29, 2006), paper 170-31. Available online at 
#' \url{http://www2.sas.com/proceedings/sugi31/170-31.pdf}.
#' @return Returns a \code{data.frame} with four columns. 
#'  The first, \code{time}, is the time index of the rows.  
#'  That is, it is a grid of points t along the time axis and
#'  these points correspond to the rows.  The next three are the 
#'  lower bound, best estimate, and upper bound, of the 
#'  correlation between the smoothed value of the covariate x(t)
#'  and the response y at each of the time points t. 
#'  We refer to the correlation function estimated here as marginal because
#' it ignores any other functional covariates (rather than 
#' trying to adjust or control for them).
#'@export
marginal.cor.funeigen <- function(object, 
                                  id,  
                                  response, 
                                  alpha=.05
                          ) {
    if (!(class(object)=="funeigen")) {
        stop("object needs to be a funeigen object.");
    }
    if (!(identical(sort(as.vector(unique(id))),sort(as.vector(object$id))))) {
        stop("The id variables may be mismatched.");
    }
    xhat <- fitted(object);
    short.id <- rownames(xhat);
    short.response <- rep(NA,length(short.id));
    for (i in short.id) {
        short.response[which(short.id==i)] <- response[which(id==i)][1];
    }
    time.point <- as.numeric(colnames(xhat));
    nsub <- length(short.id);
    correlation.estimate <- as.vector(cor(xhat,short.response));
    z.Fisher <- .5*log((1+correlation.estimate)/(1-correlation.estimate+1e-30));
    z.crit <- qnorm(1-alpha/2);
    csi.lower <- z.Fisher - z.crit*sqrt(1/(nsub-3));
    csi.upper <- z.Fisher + z.crit*sqrt(1/(nsub-3));
    correlation.lower <- (exp(2*csi.lower)-1)/(exp(2*csi.lower)+1);
    correlation.upper <- (exp(2*csi.upper)-1)/(exp(2*csi.upper)+1);
    answer <- data.frame(time=time.point,
                         r.lower.bound=correlation.lower,
                         r.estimate=correlation.estimate,
                         r.upper.bound=correlation.upper);
    return(answer);
}