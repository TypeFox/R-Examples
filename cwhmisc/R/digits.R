allDigits <- function(str) {
  k <- length(str)
  result <- logical(k)
  for(i in 1:k) {
    st <- str[i]
    ls <- nchar(st)
    ex <- strsplit(st, NULL)
    result[i] <- all(match(ex,c('0','1','2','3','4','5','6','7','8','9'),nomatch=0)>0)
  }
  result
}   ## end allDigits

isNumeric <- function (str) {
    oldop <- options(warn = -1)
    on.exit(options(oldop))
    !is.na(as.numeric(str))
}   ## end isNumeric

str2dig <- function( str ) {
    as.numeric( unlist(strsplit(str,NULL)) )
}   ## end str2dig

xToBase <- function( x, base=2 ) {
  e <- trunc(log(abs(x),base=base));  a <- x/base^e
  if ( abs(a) < 1.0 ) { e <- e-1.0;   a <- a*base
  } else { if (abs(a) > base ){ e <- e+1.0; a <- a/base }
  }
  return(list(a=a,e=e))
}   ## end xToBase
