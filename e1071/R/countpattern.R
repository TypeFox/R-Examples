countpattern <- function(x, matching=FALSE)
{
     nvar <- dim(x)[2]
     n <- dim(x)[1]

     ## build matrix of all possible binary vectors
     b <- matrix(0, 2^nvar, nvar)
     for (i in 1:nvar)
         b[, nvar+1-i] <- rep(rep(c(0,1),c(2^(i-1),2^(i-1))),2^(nvar-i))

     namespat <- b[,1]
     for (i in 2:nvar)
         namespat <- paste(namespat, b[,i], sep="")

     xpat <- x[,1]
     for (i in 2:nvar)
         xpat <- 2*xpat+x[,i]
     xpat <- xpat+1

     pat <- tabulate(xpat, nbins=2^nvar)
     names(pat) <- namespat

     if (matching)
         return(list(pat=pat, matching=xpat))
     else
         return(pat)
 }
