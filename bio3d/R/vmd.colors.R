vmd.colors <- function(n=33, picker=FALSE, ...){
  
  ## RGB numbers
  red <- c(0, 1, 0.35, 1, 1, 0.5, 0.6, 0, 1, 1,
           0.25, 0.65, 0.5, 0.9, 0.5, 0.5, 0, 0.88, 0.55, 0,
           0, 0, 0, 0.02, 0.01, 0.27, 0.45, 0.9, 1, 0.98,
           0.81, 0.89, 0.96)

  green <- c(0, 0, 0.35, 0.50, 1, 0.5, 0.6, 1, 1, 0.6,
             0.75, 0, 0.9, 0.4, 0.3, 0.5, 0, 0.97, 0.9, 0.9,
             0.9, 0.88, 0.76, 0.38, 0.04, 0, 0, 0, 0, 0,
             0, 0.35, 0.72)

  blue <- c(1, 0, 0.35, 0, 0, 0.2, 0.6, 0, 1, 0.6,
            0.75, 0.65, 0.4, 0.7, 0, 0.75, 0, 0.02, 0.02, 0.04,
            0.5, 1, 1, 0.67, 0.93, 0.98, 0.9, 0.9, 0.66, 0.23,
            0, 0, 0)

  ## Setup color indices
  max.col <- length(red)
  inds <- (1:n)
  if( n > max.col ) {
    inds <- inds %% max.col
    inds[inds==0] <- max.col
    warning( paste("Colors will be recycled: input 'n' >", max.col) )
  }

  cols <- rgb(red[inds], green[inds], blue[inds], ...)
  names(cols) <- c(1:n)
  
  if(picker) {
    ## Draw a pie chart to help with color choice
    if(n > 50) { warning("Chart will likely be crowded, set n=33 to see all colors") }
    pie(rep(1, length(cols)), labels=paste(inds, cols), col=cols, cex=0.75)
  }
  return(cols)
}

