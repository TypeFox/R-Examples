map.legend <-
function(x, y=NULL, group=NULL, data=NULL, col='black', symbol='circle', alpha=0.2, lwd=NA, lty=1, cex=1,
            title='', text=NULL, box='gray', inset=0.02, ...){

  if( is.numeric(x) ){ # if x, y given then assume proportion of the plot wndow in usr coordinates
    x <- par('usr')[1] + x*(par('usr')[2]-par('usr')[1])
    y <- par('usr')[3] + y*(par('usr')[4]-par('usr')[3])
  }
  
  if( is.integer(symbol) ){ # sort out symbols, circle as default
    pch <- symbol
  } else if( is.character(symbol) ){
    pch <- rep(0, length(symbol))
    for(i in 1:length(symbol) )
      pch[i] <- switch(symbol[i], 'cross'=3, 'star'=8, 'square'=15, 'circle'=16, 'diamond'=18, 16)
  } 
  
  if( !is.null(text) )
    gp.levels <- text
  else if( !is.null(group) ){ # now group labels
    if( is.numeric(group) )
      gp.data <- as.factor((data[, group]@data[ , 1]))
    else if( is.character(group) )
      gp.data <- as.factor(data[, which(names(data)==group)]@data[ , 1])
    gp.levels <- levels(gp.data)
  }
  if( title=='') # if no title given use grouping variable
    title <- group
  
  bgcol <- paste('#FFFFFF', as.hexmode(floor(alpha*255L)), sep='')
  
  legend(x=x, y=y, legend=gp.levels, title=title, col=col , pch=pch, bg=bgcol, box.col=box, 
         lwd=lwd, lty=lty, cex=cex, inset=inset, ...)

}
