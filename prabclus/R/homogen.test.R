"homogen.test" <-
function(distmat,ne=ncol(distmat),testdist="erdos"){
  nc <- ncol(distmat)
  vdist <- distmat[upper.tri(distmat)]
  sdist <- sort(vdist)
  distcut <- sdist[ne]
  iv <- nc
  for (i in 1:nc)
    iv <- iv - (min(distmat[-i,i])<=distcut)
  if(testdist=="erdos"){
    lambda <- exp(-2*ne/(nc-1) + log(nc))
    p <- 1-ppois(iv-1,lambda)
    if (iv>=lambda)
      psymm <- 2*pnorm((lambda-iv)/sqrt(lambda))
    else
      psymm <- 2*pnorm((iv-lambda)/sqrt(lambda))
    out <- list(p=p, p.twoside=psymm, iv=iv, lambda=lambda, distcut=distcut, ne=ne)
  } 
  if(testdist=="ling"){
    print("Computation of Ling-probabilities...")
    prob <- rep(0,ne*(nc-iv))
    dim(prob) <- c(ne,(nc-iv))
#    prob <- rep(0,ne*nc)
#    dim(prob) <- c(ne,nc)
    prob[1,2] <- 1
    n2 <- choose(nc,2)
    for (i in 2:ne)
      for(j in 3:(nc-iv))
#      for(j in 3:nc)
        prob[i,j] <- (choose(nc-j+2,2)*prob[i-1,j-2]+
          (j-1)*(nc-j+1)*prob[i-1,j-1]+(choose(j,2)-(i-1))*prob[i-1,j])/
          (n2-(i-1))
#    cat(prob)
    p2 <- p <- 0
    for(j in 2:(nc-iv))
      p <- p+prob[ne,j]
    print("finished.")
    out <- list(p=p, iv=iv, distcut=distcut, ne=ne)
  }
  out
}

