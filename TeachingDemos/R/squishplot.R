squishplot <- function(xlim,ylim,asp=1, newplot=TRUE){
  if(length(xlim) < 2) stop('xlim must be a vector of length 2')
  if(length(ylim) < 2) stop('ylim must be a vector of length 2')

  if(newplot) plot.new()

  tmp <- par(c('plt','pin','xaxs','yaxs'))

  if( tmp$xaxs == 'i' ){ # not extended axis range
	xlim <- range(xlim, na.rm=TRUE)
  } else { # extended range
	tmp.r <- diff(range(xlim, na.rm=TRUE))
	xlim <- range(xlim, na.rm=TRUE) + c(-1,1)*0.04*tmp.r
  }

  if( tmp$yaxs == 'i' ){ # not extended axis range
	ylim <- range(ylim, na.rm=TRUE)
  } else { # extended range
	tmp.r <- diff(range(ylim, na.rm=TRUE))
	ylim <- range(ylim, na.rm=TRUE) + c(-1,1)*0.04*tmp.r
  }


  tmp2 <- (ylim[2]-ylim[1])/(xlim[2]-xlim[1])

  tmp.y <- tmp$pin[1] * tmp2 * asp

  if(tmp.y < tmp$pin[2]){ # squish vertically
	par(pin=c(tmp$pin[1], tmp.y))
	par(plt=c(tmp$plt[1:2], par('plt')[3:4]))
  } else { # squish horizontally
	tmp.x <- tmp$pin[2]/tmp2/asp
	par(pin=c(tmp.x, tmp$pin[2]))
	par(plt=c(par('plt')[1:2], tmp$plt[3:4]))
  }

  return(invisible(tmp['plt']))
}

