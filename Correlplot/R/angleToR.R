angleToR <-
function(x,ifun="cos") {
#
# from angles (in radians) to a correlation matrix
# all angles w.r.t. positive x-axis.
#
   theta <- x
   ncor <- length(theta)
   CC <- matrix(rep(theta,ncor),nrow=ncor)
   DD <- CC-t(CC) # differences in angle, 0's on diagonal.
   R <- switch(ifun, cos = cos(DD) , lincos = lincos(DD), stop("angleToR: invalid interpretation function"))
   return(R)
}
