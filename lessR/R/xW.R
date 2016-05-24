xW <- function(x, w=90, indent=5) {

  n.lines <- (nchar(x) + (w-1)) %/% w

  spc <- paste(rep(" ", indent), collapse = "")

  b <- gregexpr(" ", x, fixed=TRUE)[[1]]

  y <- ""
  stop <- 0
  for (i in 1:n.lines) {
    start <- stop + 1
    if (i < n.lines) {
      k <- 1
      while (!any(b >= i*w-k)) k <- k+1 
      stop <- b[which(b>=(i*w-k))[1]]
    }
    else
      stop <- nchar(x)
    y <- paste(y, substr(x, start, stop), sep="")
    if (i < n.lines) y <- paste(y, "\n", spc, sep="") 
  }

  return(y)
}
