#' @name Inference
#' @rdname Inference
#' 
#' @title Inference summaries
#' @description Credibility intervals, estimates
#' @note \code{Inference} is a generic name for the functions documented. 
#' 
#' @param a,b,c,d Prior parameters
#' @param S,T sample sizes
#' @param x,y Observed counts
#' @param level confidence level
#' @param intervals a character vector, the intervals to be returned
#' @param parameter parameter of interest \code{"phi"} or \code{"VE"} (\code{=1-phi})
#' @param ... arguments passed to \link{IntrinsicInference} and \link{Intrinsic2Inference}
#' @return A list of confidence intervals (\code{brr_intervals}) or estimates (\code{brr_estimates})
#' @seealso \code{\link{confint.brr}}
#' 
#' @examples 
#' brr_intervals(x=4, y=5, S=10, T=10, a=0.5, b=0, c=0.5, d=0)
#' brr_intervals(x=4, y=5, S=10, T=10, a=0.5, b=0, c=0.5, d=0, intervals=c("left","equi-tailed"))
#' brr_estimates(x=4, y=5, S=10, T=10, a=0.5, b=0, c=0.5, d=0)
NULL

#' @rdname Inference
#' @importFrom TeachingDemos hpd
#' @importFrom stats setNames
#' @export
brr_intervals <- function(x, y, S, T, a=0.5, b=0, c=0.5, d=0, level=.95, intervals="equi-tailed*", ...){
  post.icdf <- function(q){
    qpost_phi(q, a, b, c, d, S, T, x, y)
  }
  hpd2 <- function(x, y, S, T, a, b, c, d, level){
    if(c+x<1){
      bounds <- c(0, post.icdf(level))
    }else{
      bounds <- hpd(post.icdf, conf=level)
    }
    return(bounds)
  }
  bounds <- sapply(intervals, function(interval) setNames(
    switch(interval, 
           left=c(0, post.icdf(level)), 
           right=c(post.icdf(1-level), Inf),
           "right*"=c(sign(x)*post.icdf(1-level), Inf),
           "equi-tailed"=post.icdf(c((1-level)/2, (1+level)/2)),
           "equi-tailed*"=c(sign(x)*post.icdf((1-level)/2),post.icdf((1+level)/2)),
           hpd=hpd2(x, y, S, T, a, b, c, d, level), 
           intrinsic=intrinsic_bounds(x, y, S, T, a, b, c, d, level, ...),
           intrinsic2=if(a==0.5 && b==0) c(NA,NA) else intrinsic2_bounds(x, y, S, T, a, b, c, d, level, ...)
    ), c("lwr", "upr")
  ), simplify=FALSE)
  if(is.null(bounds)) stop("invalid interval name")
  return(bounds)
}
#'
#' @rdname Inference
#' @export
brr_estimates <- function(x, y, S, T, a=0.5, b=0, c=0.5, d=0, parameter="phi", ...){
  if(!parameter %in% c("phi","VE")) stop("parameter must be 'phi' or 'VE'")
  estimates <- spost_phi(a, b, c, d, S, T, x, y)[c("mode","mean","Q2")]
  names(estimates)[3] <- "median"
  intrinsic <- intrinsic_estimate(x, y, S, T, a, b, c, d, ...)
  intrinsic2 <- ifelse(a==0.5 && b==0, NA, intrinsic2_estimate(x, y, S, T, a, b, c, d, ...))
  out <- c(estimates, intrinsic=intrinsic, intrinsic2=intrinsic2)
  if(parameter=="VE") out <- lapply(out, function(x) 1-x)
  return(out)
}
