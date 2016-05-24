#' @rdname csSAM-internal
make.monotone <-
function(x) {
  n=length(x)
  for(j in 2:n){ x[j]=min(x[j-1],x[j])}
  return(x)
}

