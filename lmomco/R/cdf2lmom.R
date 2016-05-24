"cdf2lmom" <-
function(r, para, fdepth=0, silent=TRUE, ...) {
   if(r < 1) { warning("r < 1, returning NA"); return(NA) }
   if(! check.fs(fdepth)) return(NA)
   lower <- par2qua(  fdepth, para)
   upper <- par2qua(1-fdepth, para)
   if(r == 1) {
      tmp <- NULL
      try(tmp <- integrate(par2qua, 0, 1, para, ...), silent=silent)
      tmp <- ifelse(is.null(tmp), NA, tmp$value)
      return(tmp)
   }
   "sfunc" <- function(j) {
      tmpA <- (-1)^j * exp(lchoose(r-2,j) + lchoose(r, j+1))
      "afunc" <- function(x, j) { Fx <- par2cdf(x, para, ...);
                                  return( Fx^(r-j-1) * (1-Fx)^(j+1) ) }
      tmpB <- NULL
      try(tmpB <- integrate(afunc, lower, upper, j=j), silent=silent)
      tmpB <- ifelse(is.null(tmpB), 0, tmpB$value)
      return(tmpA*tmpB)
   }
   tmp <- sum(sapply(0:(r-2), sfunc)) / r
   return(tmp)
}

