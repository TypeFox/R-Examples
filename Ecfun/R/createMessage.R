createMessage <- function(x, width.cutoff=45, default='x', collapse='; ', 
                          endchars='...'){
  if(sum(nchar(x))<1){
    x <- default 
  } 
  x. <- paste(x, collapse=collapse)
  nchx <- nchar(x.)
  maxch <- (width.cutoff-nchar(endchars))
  if(nchx>maxch){
    x2 <- substring(x., 1, maxch)
    x. <- paste0(x2, endchars)
  }
  x. 
}
