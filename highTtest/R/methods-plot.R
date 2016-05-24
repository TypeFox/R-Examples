setMethod("plot",
  signature=c(x="highTtest"),
  definition = function(x, ...){
    if(length(x@gammas) > 1){
      plotI_Lines(x,...)
    } else {
      plotI_Points(x,...)
    }
  }
)

plotI_Lines <- function(x,...){

  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

  xaxis <- x@gammas

  args <- list(...)

  args$x <- xaxis
  args$y <- colSums(x@CK)
  args$type <- "l"
  args$lwd <- 2
  if(is.null(args$ylim))args$ylim <- c(0, max(colSums(x@CK)))
  if(is.null(args$xlab)) args$xlab <- "FDR control level"
  if(is.null(args$ylab)) args$ylab <- "number significant"
  if(is.null(args$cex.lab)) args$cex.lab <- 1
  if(is.null(args$cex.axis)) args$cex.axis <- 1

  do.call(plot,args)

  leg <- "CK"
  j <- 1

  opts <- c("BH","ST")

  for(i in 1:2){
    temp <- slot(x,opts[i])

    if(is.null(temp)) next

    leg <- c(leg,opts[i])
    j <- j + 1
    lines(xaxis, colSums(temp), 
          col = j, 
          lty = j, 
          lwd = 2)
  }

  legend("topright",inset=c(-0.2,0),  
        lty = 1:j,  
        col = 1:j,  
        legend= leg,  
        lwd = 2)

}

plotI_Points <- function(x,...){

  xaxis <- x@gammas

  args <- list(...)

  args$x <- xaxis
  args$y <- colSums(x@CK)
  args$type <- "p"
  args$lwd <- 2
  if(is.null(args$ylim))args$ylim <- c(0, max(colSums(x@CK)))
  if(is.null(args$xlab)) args$xlab <- "FDR control level"
  if(is.null(args$ylab)) args$ylab <- "number significant"
  if(is.null(args$cex.lab)) args$cex.lab <- 1
  if(is.null(args$cex.axis)) args$cex.axis <- 1

  do.call(plot,args)

  leg <- "CK"
  j <- 1

  opts <- c("BH","ST")

  for(i in 1:2){
    temp <- slot(x,opts[i])

    if(is.null(temp)) next

    leg <- c(leg,opts[i])
    j <- j + 1
    points(xaxis, colSums(temp), 
           col = j, 
           pch = j, 
           lwd = 2)
  }

  legend("topright",inset=c(-0.2,0),  
        pch = 1:j,  
        col = 1:j,  
        legend= leg,  
        lwd = 2)

}
