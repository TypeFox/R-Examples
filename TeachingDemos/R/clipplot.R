clipplot <- function(fun, xlim=par('usr')[1:2],
	ylim=par('usr')[3:4] ){
  old.par <- par(c('plt','xpd'))


  if( length(xlim) < 2 ) stop('xlim must be a vector with at least 2 elements')
  if( length(ylim) < 2 ) stop('ylim must be a vector with at least 2 elements')

  xl <- range(xlim)
  yl <- range(ylim)

  pc <- cnvrt.coords(xl,yl)$fig

  box(col='#00000000') # works better with this, don't know why
  par(plt=c(pc$x,pc$y),xpd=FALSE)
  box(col='#00000000') # same

  fun

  par(old.par)
  box(col='#00000000') # need to plot something to reset
  invisible(NULL)
}

