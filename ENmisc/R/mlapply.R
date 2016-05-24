mlapply <- function(lol,FUN,...){
  llol <- sapply(lol, length)
  nrows <- llol[1]
  if (any(llol != nrows)) stop("lists not of same length")
  
  arglists <- lapply(as.list(1:length(lol[[1]])),
                     function(i) lapply(lol, `[[`, i))
  names(arglists) <- names(lol[[1]])
  lapply(arglists, function(x)do.call(FUN,c(x,...)))
}
