# Rotate the data and then smooth the diagonal elements. We use quadratic terms on either direction, rather than only orthogonal to the diagonal.
# xout: a matrix of two columns containing the diagonal elements.

RotateLwls2DV2 <- function(bw, kern='epan', xin, yin, win=NULL, xout) {
  
  if (length(bw) == 1){
    bw <- c(bw, bw)
  }  

  if (missing(win) || is.null(win)){
    win <- rep(1, length(xin))
  }

  if (  is.vector(xout)){
   xout =  matrix(c(xout,xout),ncol= 2)
  }
    
  fit <- Rrotatedmullwlsk(bw, kern, t(xin), yin, win, t(xout), npoly=1, bwCheck=FALSE)
  
  if (any(is.nan(fit)))
    stop('Something wrong with the rotate smoothed results')
    
  return(fit)
}
