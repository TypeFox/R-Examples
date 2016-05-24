laterhist <- function (data, catch="Food", hand="Hand", col = 1:nlevels(data[[hand]]), ylim = NULL
                       , ylab = "Number of grips", main="Type of grips regarding to the performed task"
                       , legend.text = FALSE, beside = TRUE, legendlocation=TRUE, cex=1, pt.cex=2, pch=15)
{
  graphlaterhist<-barplot(x<-table(data[[hand]],data[[catch]]), names.arg = levels(data[[catch]]), beside = beside, ylab=ylab, main=main, legend.text = legend.text, col=col, ylim=ylim, xpd=FALSE)

  #Legend
  if (legendlocation == TRUE) {
      message("Click where you want to place the legend")
      legendplace <- locator(1)
      legend(legendplace$x,legendplace$y,as.vector(levels(data[[hand]])),col=col,bty="n",pch=pch, cex=cex, pt.cex=pt.cex)
  } else {
    }
}

