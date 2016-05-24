tpe1d <- function(d,verbose=TRUE) {
  hc <- hclust(d,"single")
  n <- length(hc$height)
  d <- as.matrix(d)
  ind <- list()
  x <- list()
  for (i in 1:n) {
    if (verbose) {
      message("Iteration ",i," of ",n)
    }
    dmin <- hc$height[i]
    pair <- hc$merge[i,]
    if (max(pair)<0) {
      ind[[i]] <- -pair
      x[[i]] <- rbind(0,dmin)
    } else if (min(pair)<0) {
      c1 <- -min(pair)
      c2 <- ind[[max(pair)]]
      ind[[i]] <- c(c1,c2)
      x[[i]] <- align1d(d[c1,c2,drop=FALSE],matrix(0,1,1),x[[max(pair)]],dmin)
    } else {
      c1 <- ind[[pair[1]]]
      c2 <- ind[[pair[2]]]
      ind[[i]] <- c(c1,c2)
      x[[i]] <- align1d(d[c1,c2,drop=FALSE],x[[pair[1]]],x[[pair[2]]],dmin)
    }
  }
  x <- x[[n]][match(1:(n+1),ind[[n]]),,drop=FALSE]
  rownames(x) <- hc$labels
  x
}

align1d <- function(d,x1,x2,dmin) {
  n1 <- dim(x1)[1]
  n2 <- dim(x2)[1]
  obj <- function(x) {
    td <- as.matrix(dist(x))[1:n1,(n1+1):(n1+n2),drop=FALSE]
    sum((td-d)^2)
  }
  x1r <- reflect1d(x1)
  x2r <- reflect1d(x2)
  t1 <- min(x2)-max(x1)-dmin
  t2 <- max(x2)-min(x1)+dmin
  tx1 <- rbind(x1+t1,x2)
  tx2 <- rbind(x1+t2,x2)
  tx3 <- rbind(x1r+t1,x2)
  tx4 <- rbind(x1r+t2,x2)
  tx5 <- rbind(x1+t1,x2r)
  tx6 <- rbind(x1+t2,x2r)
  tx7 <- rbind(x1r+t1,x2r)
  tx8 <- rbind(x1r+t2,x2r)
  val <- c(obj(tx1),obj(tx2),obj(tx3),obj(tx4),obj(tx5),obj(tx6),obj(tx7),obj(tx8))
  switch(which.min(val),tx1,tx2,tx3,tx4,tx5,tx6,tx7,tx8)
}

reflect1d <- function(x) {
  max(x)-x+min(x)
}
