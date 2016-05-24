eigenCAR1 <- function(nr,nc)
{

eigenrow <- eigenRW1(nr)
eigencol <- eigenRW1(nc)

list( values = kronecker( rep(1,nc), eigenrow$values) +
        kronecker(eigencol$values, rep(1,nr)) ,
      vectors = kronecker( eigencol$vectors, eigenrow$vectors) )

}

