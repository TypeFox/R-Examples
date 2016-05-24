### iplotLegend.R

iplotLegend <- function(colramp = IDPcolorRamp,
                        ncol = NULL, cex.axis = par("cex.axis"),
                        border = FALSE, mar = c(0,0,0,3), las = 1, ...)
  ## internal functions for ipairs and ilagplot
  ## Version  2009-03-26
{
      ##-------------------
      ## drawing legend
      ## calculate colcut only when common ncol is defined
      if(is.null(ncol)) {
          ## warning(c("Separate color coding for each image\n"),
          ##         call. = FALSE)
          border <- FALSE
          mycol <- c(par("bg"),colramp(100))
      }else {
          ncol <- ceiling(ncol)
          mycol <- c(par("bg"),colramp(ncol))
      }

      lev <- 0:length(mycol)
      par(mar=mar, las=las, cex.axis=cex.axis, ...)
      plot.new()
      plot.window(xlim=c(0, 1), ylim=range(lev,na.rm=TRUE),
                  xaxs="i", yaxs="i")

      if(border)
          rect(0, lev[-length(lev)], 1, lev[-1], col=mycol)
      else
          rect(0, lev[-length(lev)], 1, lev[-1],
               col=mycol, border=mycol)
      box()
      if(is.null(ncol)) {
          axis(side=4, at=c(0.5,length(mycol)-0.5), labels=c("0","max"))
      }else {
          ap <- pretty(lev)
          axis(side=4, at=ap+0.5, labels=paste(ap))
      }

} ## PlotLegend

