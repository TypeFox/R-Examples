innovationVariance <- function(z, method=c("AR", "Kolmogoroff"), ...){
  wm <- match.arg(method, c("AR", "Kolmogoroff"))
  if (length(z)<10) 
    stop("Series length should be at least 10!")
  if (wm=="AR"){
    sigmaSq <- ar(z, aic=TRUE, method="burg")$var.pred
  } else {
    plotQ <- ifelse(length(grep("plot", c(substitute(do.call(...)))))>0, TRUE
                    , FALSE)
    if(plotQ) {
      sdf <- spec.pgram(z, ...)$spec
    } else {
      sdf <- spec.pgram(z, plot=FALSE, ...)$spec
    }
    sigmaSq <- exp(2*sum(log(2*pi*sdf))/length(z))/(2*pi)
  }
  sigmaSq
}