plotTD <- function(bHat, sdVec, times) {

  nTimes <- length(times)

  intv <- 1.96*sdVec

  maxy <- max(bHat + intv)
  miny <- min(bHat - intv)

  col <- 2L
  graphics::plot(x = times, 
                 y = bHat[,1L], 
                 col = col, lty = 1, lwd = 2, 
                 xlim = c(times[1]*0.9, times[nTimes]*1.1), 
                 ylim = c(miny*0.9,maxy*1.1), 
                 xlab = "time", ylab = "beta(t)",
                 main = "Parameter Estimates",
                 graphics::arrows(times, y0 = bHat[,1L]-intv[,1L], 
                                  times, y1 = bHat[,1L]+intv[,1L],
                                  code = 3, angle = 90, 
                                  length= 0.1, col = col))

  graphics::lines(times, bHat[,1L], type = "l", col = col)
  j <- 2L
  col <- col + 1L
  while( j <= ncol(bHat) ) {
    graphics::points(x = times, 
                     y = bHat[,j], 
                     col = col,
                     graphics::arrows(times, y0 = bHat[,j]-intv[,j], 
                                      times, y1 = bHat[,j]+intv[,j],
                                      code = 3, angle = 90,  
                                      length= 0.1, col=col))
    graphics::lines(times, bHat[,j], type = "l", col = col)
    j <- j + 1L
    col <- col + 1L
  }

  graphics::legend('topright', 
                   legend = paste("beta",0L:(ncol(bHat)-1L),sep=""),
                   lty = 1,
                   col = 2L:(col-1L), bty = 'n', cex = 0.75)

  return()

}
