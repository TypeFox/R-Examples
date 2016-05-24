make.grid <-
function(x,y,z,byx,byy,xlim,ylim,
                      fun=function(x) sum(x,na.rm=T)){
  # the origin of the grid depends on the xlim[1] and ylim[1]
  
  # X and Y midpoints 
  X <- seq(xlim[1],xlim[2],by=byx)
  Y <- seq(ylim[1],ylim[2],by=byy)
  
  # find which grid cells each datapoint falls into
  i <- which(x>=xlim[1]-0.5*byx & x<xlim[2]+0.5*byx & y>=ylim[1]-0.5*byy & y<ylim[2]+0.5*byy)
  xi <- findInterval(x[i],X-0.5*byx)
  yi <- findInterval(y[i],Y-0.5*byy)
  
  # aggregate the data
  grid0 <- tapply(z[i],list(X[xi],Y[yi]),fun)
  
  # expand the grid to include null values
  Xi <- match(X,rownames(grid0))
  Yi <- match(Y,colnames(grid0))
  grid1 <- grid0[Xi,Yi]
  dimnames(grid1) <- list(X,Y)
  
  return(grid1)
}

