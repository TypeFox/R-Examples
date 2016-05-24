`compHclust` <-
function(x,xhc) {
  if ((!is.matrix(x))|(!is.numeric(x))) {
    stop("'x' must be a numeric matrix")
  }
  if (any(is.na(as.vector(x)))) {
    stop("'x' has missing values, please impute missing values")
  }
  if (class(xhc)!="hclust") {
    stop("'xhc' must be of class 'hclust'")
  }
  if (ncol(x)!=length(xhc$order)) {
    stop("'xhc' must be a clustering of the columns of 'x'")
  }
  x.hat <- as.matrix(x)%*%A.chc(xhc)
  return(list("x.prime"=x-x.hat,"gene.imp"=RGI.chc(x,x.hat)))
}

