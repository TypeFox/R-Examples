"coin.spin" <-
function( para, r=0, n=0, n.points=199)
{   # compute posterior density for n coin spins of which r heads, with prior in para
  para <- as.matrix(para)
  para[,1] <- para[,1]/sum(para[,1])
  k <- nrow(para)
  w <- rep(0,k)
  for (i in 1:k) 
	w[i] <- (para[i, 1] * beta(para[i, 2] + r, para[i, 3] + n - r))/
        beta(para[i, 2], para[i, 3])
  para[,1] <- w/sum(w)
  x <- c(1:n.points)/(n.points +1)
  y <- rep(0,n.points)
  for (i in 1:k)  
    y <- y + para[i,1]*dbeta(x, shape1=para[i,2]+r, shape2=para[i,3]+n-r )
  list( x=x, y=y )
}

