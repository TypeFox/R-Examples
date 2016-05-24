diffexprm <-
function(x, batch, y, method=c("fabatch", "combat", "sva", "meancenter", "standardize", "ratioa", "ratiog", "none")) {

  if(!(method %in% c("fabatch", "combat", "sva", "meancenter",
    "standardize", "ratioa", "ratiog", "none")))
    stop("Input parameter 'method' has to be one of the following:\n'fabatch',  'combat', 'fsva', 'meancenter', 'standardize', 'ratioa', 'ratiog', 'none'.")

  batchessuited <- apply(table(batch, y), 1, function(x) min(x)>1)
  xbrlist <- lapply(levels(batch)[batchessuited], function(x2) 
    ba(x=x[batch!=x2,], y=y[batch!=x2], batch=factor(as.numeric(factor(as.numeric(batch[batch!=x2])))), method=method)$xadj)

  diffexprmAfterBR(x, xbrlist, y, batch, batchessuited)

}
