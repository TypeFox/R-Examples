matricize <- function (a, b, bins, radius, aspect, type=NULL) {
  if(type == "hist") {
    x.mat <- a 
    cols <- bins
    names.b <- "frequency"
    names.a <- b
    } else {
      x <- hexbin (a[ ,1] ~ b[ ,1], xbins=bins, shape=aspect)
      x.vec <- vector (length = x@dimen[1] * x@dimen[2])
      x.vec[x@cell] <- x@count
      x.mat <- data.frame(Xcoord = rep(seq(1:(x@dimen[2]))-1,x@dimen[1]),
                          Ycoord = rev(rep(c(0:(x@dimen[1]-1)),each=x@dimen[2])),
                          V3     = x.vec)
      cols <- x@dimen[1]
      names.b <- names(b)
      names.a <- names(a)
    }
  
  #now calculate the number of empty bins around each bin
  x.asmat <- matrix (x.mat$V3, ncol = cols)
  x.asmat <- x.asmat[ ,ncol (x.asmat):1] ##reverse-order
  x.dist  <- apply ( (x.mat[ ,1:2]+1), 1, function(y){
    x.cig <- max (x.mat$Xcoord + 1)
    y.cig <- max (x.mat$Ycoord + 1)
    length (which (x.asmat[ (max ( (y[1] - radius),1):min (x.cig,(y[1] + radius) ) ),
                            (max ( (y[2] - radius),1):min (y.cig,(y[2] + radius) ) )]==0) )
  })
  x.mat[,4] <- x.dist
  names(x.mat)[3] <- paste(names.b,names.a,sep="_")
  names(x.mat)[4] <- paste(names.b,names.a,"dist",sep="_")
  x.mat
}
