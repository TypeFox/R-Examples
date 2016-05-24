clustcentroid <-
function(x,y, m=1) {
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
  s = t(as.matrix(x))%*%(u)
  centers = sweep(t(s),1,colSums(u),"/")
  colnames(centers) = names(x)
  rownames(centers) = colnames(u)   
	return(centers)
}

