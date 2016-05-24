map.density2 <-
  function(data, grid.size=300){
    
    #melt <- NULL
# depends on packages reshape and ks    
    
    # extract locations over which kernel to be fitted
    xmat <- as.matrix(cbind(coordinates(data)[ , 1], coordinates(data)[ , 2]))
    
    # do the smoothing
    x.H <- Hpi(xmat[runif(1000, 1, length(xmat[ , 1])), ])
    x.sm <- kde(xmat, H=x.H, gridsize=grid.size)
    
    # setup the topology for mapping
    ep <- x.sm$eval.points # for convenience
    cell.offset <- c(min(ep[[1]]), min(ep[[2]])) # the start of the grid
    cell.size <- c(ep[[1]][2]-ep[[1]][1], ep[[2]][2]-ep[[2]][1])
    ncells <- c(length(ep[[1]]), length(ep[[2]]))
    gt <- GridTopology(cell.offset, cell.size, ncells)
    
    f1=function(m){ # rotates the smoothed grid to match the Spatial way of doing things
      m=t(m)[ncol(m):1,]
      m=t(m)[ncol(m):1,] # this duplication seems to be needed...
      apply(t(m),1,rev)
    }
    
    # create the SpatialGrid
    x.obs <- kde(x=x.sm$x, H=x.sm$H, eval.points=x.sm$x)$estimate
    maxht <- ifelse(max(x.obs,na.rm=T)>quantile(x.obs,prob=1,na.rm=T), max(x.obs), quantile(x.obs,prob=1,na.rm=T)) + 0.1
    hts <- c(0,quantile(x.obs,probs=seq(from=0.05, to=0.95, length.out=9),na.rm=T),maxht)
    # this next line should make sure that hts is strictly increasing!
    hts <- hts+seq(from=1e-14, to=1e-14*length(hts),length.out=length(hts))
        
    # generate the colours
    sm.col <- c('transparent',rainbow(9, end=1/6, alpha=0.5))
    
    # plot the densties 
    .filled.contour(x=x.sm$eval.points[[1]], y=x.sm$eval.points[[2]], z=x.sm$estimate,
                    levels=hts,col=sm.col)  
    
  }
