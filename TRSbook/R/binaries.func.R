bin2dec <- function(x) {

  bin2dec.int <- function(binaryvector) {   
    sum(2^(which(rev(binaryvector)==TRUE)-1)) 
  }

  bin2dec.frac <- function(binaryvector) {
    sum(2^-(which(binaryvector==TRUE))) 
  }
  
  x <- as.character(x)
  res <- strsplit(x,".",fixed=TRUE)[[1]]
  a <- as.numeric(strsplit(res[1],"",fixed=TRUE)[[1]])
  b <- as.numeric(strsplit(res[2],"",fixed=TRUE)[[1]])
  return(bin2dec.int(a)+bin2dec.frac(b))
}

dec2bin <- function(x,prec=52) {

  dec2bin.int <- function(x) {
    as.integer(paste(rev(as.integer(intToBits(x))), collapse=""))
  }

  dec2bin.frac <- function(x,prec=52) {
    res <- rep(NA,prec)
    for (i in 1:prec) {
      res[i] <- as.integer(x*2)
      x <- (x*2) %% 1
    }
    return(paste(res,collapse=""))
  }
  
  x <- as.character(x)
  res <- strsplit(x,".",fixed=TRUE)[[1]]
  int <- res[1]
  fract <- res[2]
  if (is.na(fract)) {
    fract <- 0 
    prec <- 0
  }
  left <- dec2bin.int(as.numeric(int))
  right <- dec2bin.frac(as.numeric(paste("0.",fract,sep="")),prec)
  if (right == "0") {return(left)} else {return(paste(left,right,sep="."))}
}
