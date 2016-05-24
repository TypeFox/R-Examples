"intmew" <-
function(iopt,tl,tu,til=1e-4) {
b1 <- 0
z  <- .Fortran("srintmw",iwgt=as.integer(iopt),tl=as.double(tl),tu=as.double(tu),
      b1=as.double(b1),til=as.double(til),sum=double(1))
z$sum}

