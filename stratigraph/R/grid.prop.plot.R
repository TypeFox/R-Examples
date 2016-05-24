grid.prop.plot <- function(x, y = NULL,
                           increasing.down = FALSE){

  if(is.null(y)) y <- 1:nrow(x)
  if(increasing.down){
    bot <- y[length(y)]
    top <- y[1]
  }else{
    bot <- y[1]
    top <- y[length(y)]
  }
  
  #cat(bot, top)
  xpc <- x / rowSums(x) * 100
  grid.newpage()
  pushViewport(viewport())
  pushViewport(viewport(width = unit(1, 'npc') - unit(6, 'lines'),
                        height = unit(1, 'npc') - unit(9, 'lines'),
                        xscale = c(0,100), yscale = c(bot, top)))
  grid.rect()
  for(i in 1:ncol(x)){
    if(i == 1){
      grid.polygon(c(0, xpc[,i], 0),
                   c(y[1], y, y[length(y)]),
                   default.units = 'native',
                   gp = gpar(col = 'black', fill = 'black'))
    }else if(i == ncol(x)){
      prev <- rowSums(as.matrix(xpc[,-i]))
      grid.polygon(c(100, prev, 100),
                   c(y[1], y, y[length(y)]),
                   default.units = 'native',
                   gp = gpar(col = 'black', fill = 'white'))
    }else{
      prev <- rowSums(as.matrix(xpc[,1:(i-1)]))
      this <- rowSums(as.matrix(xpc[,1:i]))
      grid.polygon(c(this, prev[length(prev):1]),
                   c(y, y[length(y):1]),
                   default.units = 'native',
                   gp = gpar(col = i, fill = i))
    }         
  }
  grid.xaxis()
  grid.yaxis()
}