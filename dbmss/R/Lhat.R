Lhat <-
function(X, r = NULL, ReferenceType = "", NeighborType = "", CheckArguments = TRUE) {
  if (CheckArguments) {
    CheckdbmssArguments()
  }
  
  K <- Khat(X, r, ReferenceType, NeighborType)
  Columns <- names(K)[-1]
  for(i in Columns) {
    K[[i]] <- sqrt(K[[i]]/pi)-K$r
  }
  attr(K, "ylab") <- attr(K, "yexp") <- quote(L(r))
  attr(K, "fname") <- "L"
  return (K)
}
