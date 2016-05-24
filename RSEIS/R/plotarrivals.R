plotarrivals<-function(x, THEORY, add=FALSE)
  {
    if(missing(add)) { add=FALSE }
    if(!add){
      plot(range(c(0,x)), c(0, 2.4) , type='n')


    }
    klay = dim(THEORY$trefrac)[1]
    for(n in 1:klay)
      {
        lines(x, -THEORY$trefrac[n,], col=n)
      }

    klay = dim(THEORY$treflex)[1]
    
    for(n in 1:klay)
      {
        lines(x, -THEORY$treflex[n,], col=n)
      }
    lines(x, -THEORY$tair, col=rgb(.5, .5, 1) )

  }


