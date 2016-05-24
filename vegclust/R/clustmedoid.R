clustmedoid <-
function(x, y, m=1) {
  if(inherits(y,"vegclust") || inherits(y, "vegclass")) {
    u = as.matrix(y$memb)^y$m
    if(y$method=="NC"||y$method=="HNC"||y$method=="NCdd"||y$method=="HNCdd") {
      u = u[,-ncol(u)]
    }
  } else if(is.null(dim(y))) {
    u = as.memb(y)
    colnames(u) = levels(as.factor(y))
  } else {
    u = as.matrix(y)^m
  }
  c = ncol(u)
  med = numeric(c)
  if(!inherits(x,"dist")) {
    d = as.matrix(dist(x))    
  }
	for(k in 1:c) {
	  med[k] = which.min((u[,k]^m)%*%d)
	}   
  names(med) = row.names(as.matrix(d))[med]
	return(med)
}

