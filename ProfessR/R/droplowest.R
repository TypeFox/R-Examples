droplowest<-function(z)
  {
##########   z is a matrix of midterm grades

    ### index of the minimum grade in original matrix

    minind =  apply(z, 1, which.min)
    ###  sort grades

    
    z3 = t( apply(z, 1, sort))
  ##  cbind(z3, minind, z)

    ##########   best grades remain, but order is re-arranged
    best = z3[, -1]

    midgrade = apply(best, 1, mean, na.rm=TRUE)

    return(list(minind=minind , best=best, midgrade=midgrade ))

  }
