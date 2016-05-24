plot.bounds <- function(x, row, column, labels = TRUE, order =
                        NULL, intersection = TRUE, xlab, ylab, col = par("fg"),
                        lty = par("lty"), lwd = par("lwd"), ...){ 
  
  if(class(x) != "bounds"){
    stop("'x' must be output from 'bounds'")
  }

  bounds <- x$bounds
  idx <- paste(row, ".", column, sep="")
  if (!(idx %in% names(bounds))){
    stop("selected row or column bounds not in 'x' - please choose a different row or column")
  }

  if(all(is.na(bounds[[idx]]))){
    stop("selected row or column bounds not in 'x' - no precincts satisfy threshold")
  }

  "%wo%" <- function(x,y){x[!x %in% y]}
  threshold <- 100*x$threshold

  if(is.null(order)){
    order <- (1:nrow(bounds[[idx]]))/(nrow(bounds[[idx]])+1)
    xl <- 0:1
    axes <- FALSE
  }
  else {
    xl <- range(order)
    axes <- TRUE
  }
  if (missing(xlab)) {
    xlab <- paste("Precincts with at least", threshold, "% ", row)
  }
  if (missing(ylab)) {
    ylab <- paste("Proportion ", column, sep="")
  }
  plot(xl, 0:1, type = "n", xlab = xlab, ylab = ylab,
       axes = axes, ...)
  axis(side = 2, labels=TRUE)  

  segments(order, bounds[[idx]][,"lower"], order, bounds[[idx]][,"upper"], col = col,
           lty = lty, lwd = lwd)

  if(labels){
    text(order, bounds[[idx]][,"upper"]+.02,
         rownames(bounds[[idx]]), cex=0.4)
  }
  
  if(intersection){
    if(!all(is.na(x$intersection[[row]]))){
      abline(h=c(x$intersection[[row]][1],
               x$intersection[[row]][2]), col="grey80")
    }
  }
}
