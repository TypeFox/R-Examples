#-------------------------------------------------------------------#
#                                                                   #
#                     PARETO CHART                                  #
#                                                                   #
#-------------------------------------------------------------------#

pareto.chart <- function(x, ylab = "Frequency", ylab2 = "Cumulative Percentage", xlab, cumperc = seq(0, 100, by = 25), ylim, main, col = heat.colors(length(x)), plot = TRUE, ...)
{
  call <- match.call(expand.dots = TRUE)
  varname <- deparse(substitute(x))
  x <- as.table(x)
  if (length(dim(x))>1) 
     stop("only one-dimensional object (table, vector, etc.) may be provided")
  # 
  x <- sort(x, decreasing = TRUE, na.last = TRUE)
  cumsum.x <- cumsum(x)
  cumperc <- cumperc[cumperc >= 0 & cumperc <= 100]
  q <- quantile(seq(0, max(cumsum.x, na.rm = TRUE), 
                      by = max(cumsum.x, na.rm = TRUE) / 100), 
                cumperc/100)
  tab <- cbind(x, cumsum.x, 
               x/max(cumsum.x, na.rm = TRUE)*100, 
               cumsum.x/max(cumsum.x, na.rm = TRUE)*100) 
  colnames(tab) <- c("Frequency", "Cum.Freq.", 
                     "Percentage", "Cum.Percent.")
  names(dimnames(tab)) <- c("", paste("\nPareto chart analysis for", varname))
  # tab <- as.table(tab)
  
  if(plot)
  { if(missing(xlab)) xlab <- ""
    if(missing(ylim)) ylim <- c(0, max(cumsum.x, na.rm = TRUE)*1.05)
    if(missing(main)) main <- paste("Pareto Chart for", varname)
    if(missing(col))  col <- heat.colors(length(x))
    # set las and mar if not provided by user
    w <- max(sapply(names(x), nchar))
    if(is.null(call$las)) las <- 3 else las <- call$las
    if(is.null(call$mar))
      { if (las==1) mar <- c(1,1,0,2)  
        else        mar <- c(log(max(w),2),0,0,2) }
    else mar <- call$mar
    oldpar <- par(mar = pmax(par("mar")+mar,c(4.1,4.1,3.1,4.1)), 
                  las = las, 
                  cex = qcc.options("cex"),
                  no.readonly = TRUE)
    on.exit(par(oldpar))
    pc <- barplot(x, width = 1, space = 0.2, main = main, 
                  ylim = ylim, ylab = ylab, xlab = xlab, yaxt = "n", col = col, ...)
    # adding line for percentage level overwrite bars...
    abline(h = q, col = "lightgrey", lty = 3)
    # ... so we redraw bars (not nice but works!)
    rect(pc-0.5, rep(0,length(x)), pc+0.5, x, col = col)
    lines(pc, cumsum.x, type = "b", cex = 0.7, pch = 19)
    box()
    axis(2, las = 3)
    axis(4, at = q, las = 3, labels = paste(cumperc, "%", sep = ""))
    mtext(ylab2, 4, line = 2.5, las = 3) 
  }    
  
  return(tab)
}
