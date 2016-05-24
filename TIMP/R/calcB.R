"calcB" <-
function (k) 
{
    lenk <- length(k) 
    m <- rep(0, lenk * lenk)
    m <- as.matrix(.C("calcB", m = as.double(m), as.double(k), 
                      as.integer(lenk), PACKAGE="TIMP")$m)
    dim(m) <- c(lenk, lenk)
    m 
  }

