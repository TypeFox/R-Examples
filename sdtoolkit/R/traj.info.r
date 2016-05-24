#traj.info - based of of doung's... box.info?  This function is modified to work
#on a trajectory and thus be list based, rather than on a single box.

#takes a trajectory and calculates statistics on it, such as mass, mean,
#relative and marginal coverage, etc.

`traj.info` <-
function(x, y, box.seq, npts=NA,ninter=NA)
{
  m <- box.seq$num.class
  d <- ncol(x)
  n <- nrow(x)

  x.ind <- rep(TRUE, n)
  xy.list <- list()

  for (k in 1:m){
    x.ind.curr <- x.ind
    box.curr <- box.seq$box[[k]]

		#browser()
    for (j in 1:d){
      x.ind.curr <- x.ind.curr & (x[,j]>= box.curr[1,j]) & (x[,j] <= box.curr[2,j])
		}
		
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
    xy.list$dimlist[[k]] <- dimchecker(x,box.curr)


#    x.ind <- x.ind & !x.ind.curr

  }
  
  precov <- xy.list$mass * xy.list$y.mean
  xy.list$relcoverage  <- precov*nrow(x) / sum(y) #relative coverage
  xy.list$marcoverage  <- precov*nrow(x) / ninter #coverage relative to original
#  xy.list$abscoverage  <- precov*npts/ninter
  return (xy.list)
}

