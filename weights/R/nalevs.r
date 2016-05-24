nalevs <- function(x, naset=NULL, setmid=NULL, set1=NULL, set0=NULL, setmean=NULL, weight=NULL){
  x[x %in% naset] <- NA
  q <- (x %in% setmid)
  r <- (x %in% set1)
  s <- (x %in% set0)
  t <- (x %in% setmean)
  x[q] <- NA
  x[r] <- NA
  x[s] <- NA
  x[t] <- NA
  x <- as.numeric(x)
  x <- (x-range(x, na.rm=TRUE)[1])/range((x-range(x, na.rm=TRUE)[1]), na.rm=TRUE)[2]
  x[q] <- .5
  x[r] <- 1
  x[s] <- 0
  if(!is.null(weight))
    x[t] <- wtd.mean(x, weight, na.rm=TRUE)
  else
    x[t] <- mean(x, na.rm=TRUE)
  x
}
