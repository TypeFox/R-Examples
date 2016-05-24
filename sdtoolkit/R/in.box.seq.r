`in.box.seq` <-
function(x, y, box.seq)
{
  m <- box.seq$num.class
  d <- ncol(x)
  n <- nrow(x)
  
  x.ind <- rep(TRUE, n)
  xy.list <- list()

  for (k in 1:m)
  {
    x.ind.curr <- x.ind    
    box.curr <- box.seq$box[[k]]
    
    for (j in 1:d)
      x.ind.curr <- x.ind.curr & (x[,j]>= box.curr[1,j]) & (x[,j] <= box.curr[2,j])
    
    x.curr <- x[x.ind.curr & x.ind,]
    box.mass.curr <- sum(x.ind.curr)/n
  
    xy.list$x[[k]] <- x.curr
    if (!missing(y))
    {
      y.curr <-  y[x.ind.curr & x.ind]
      y.mean.curr <- mean(y.curr)
      xy.list$y[[k]] <- y.curr
      xy.list$y.mean[[k]] <- y.mean.curr
    }
    xy.list$box[[k]] <- box.curr
    xy.list$mass[[k]] <- box.mass.curr
   

    x.ind <- x.ind & !x.ind.curr
  }
  return (xy.list)
}

