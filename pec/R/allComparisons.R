allComparisons <- function(x,sep=" <= "){
  n <- length(x)
  pro <- rep(x,(n-1):0)
  contra <- x[unlist(lapply(2:n,function(i)i:n))]
  out <- lapply(1:length(pro),function(i)c(pro[i],contra[i]))
  names(out) <- sapply(out,paste,collapse=sep)
  out
}
