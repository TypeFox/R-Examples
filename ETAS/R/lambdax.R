
lambdax <- function(rt, rx, ry, theta, revents)
{
  if (rt < revents[1,1])
   return(0)
  theta <- sqrt(theta)
  storage.mode(revents) <- "double"
  .Call("lambdax", as.double(rt), as.double(rx), as.double(ry), 
        as.double(theta), revents, PACKAGE="ETAS")[[1]]
}



