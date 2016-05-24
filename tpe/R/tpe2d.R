tpe2d <- function(d,verbose=TRUE) {
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
      x[[i]] <- rbind(c(0,0),c(dmin,0))
    } else if (min(pair)<0) {
      c1 <- -min(pair)
      c2 <- ind[[max(pair)]]
      ind[[i]] <- c(c1,c2)
      if (length(c2)==2) {
        x[[i]] <- cmdscale(d[ind[[i]],ind[[i]]])
      } else {
        x[[i]] <- align2d(d[c1,c2,drop=FALSE],matrix(0,1,2),x[[max(pair)]],dmin)
      }
    } else {
      c1 <- ind[[pair[1]]]
      c2 <- ind[[pair[2]]]
      ind[[i]] <- c(c1,c2)
      x[[i]] <- align2d(d[c1,c2,drop=FALSE],x[[pair[1]]],x[[pair[2]]],dmin)
    }
  }
  x <- x[[n]][match(1:(n+1),ind[[n]]),]
  rownames(x) <- hc$labels
  x
}

align2d <- function(d,x1,x2,dmin) {
  n1 <- dim(x1)[1]
  n2 <- dim(x2)[1]
  x1 <- scale(x1,scale=FALSE)
  x2 <- scale(x2,scale=FALSE)
  obj <- function(par) {
    td <- as.matrix(dist(rigid2d(x1,x2,par)))[1:n1,(n1+1):(n1+n2)]
    sum((td-d)^2)
  }
  pen <- function(par) {
    tdmin <- min(as.matrix(dist(rigid2d(x1,x2,par)))[1:n1,(n1+1):(n1+n2)])
    if (dmin>tdmin) {
      Inf      
    } else {
      (tdmin-dmin)^2
    }
  }
  par0 <- c(0,dmin+2*max(dist(x1),dist(x2)),0)
  par <- sumt(obj,pen,par0)
  rigid2d(x1,x2,par)
}

rot2d <- function(theta) {
  rbind(c(cos(theta),sin(theta)),c(-sin(theta),cos(theta)))
}

rigid2d <- function(x1,x2,par) {
  rbind(x1,sweep(x2%*%rot2d(par[1]),2,par[-1]))
}

sumt <- function(obj,pen,par0) {
  penobj <- function(lambda) {
    function(par) {
      obj(par)+lambda*pen(par)
    }
  }
  lambda <- 1
  repeat {
    par <- optim(par0,penobj(lambda))$par
    if (max(abs(par-par0))<sqrt(.Machine$double.eps)) {
      break
    }
    par0 <- par
    lambda <- 10*lambda
  }
  par
}
