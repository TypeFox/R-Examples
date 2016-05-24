na.levs <- function(x, naset=NULL, setmid=NULL, set1=NULL, set0=NULL, setmean=NULL, setmedian=NULL, setmode=NULL, weight=NULL)
  UseMethod("na.levs")



na.levs.factor <- function(x, naset=NULL, setmedian=NULL, setmode=NULL, weight=NULL){
  x[x %in% naset] <- NA
  q <- (x %in% setmedian)
  r <- (x %in% setmode)
  x[q] <- NA
  x[r] <- NA
  x <- drop.levels(x, reorder=FALSE)
  if(!is.null(weight))
    x[q] <- names(table(x))[eval(sapply(1:length(table(x)), function(g) sum(wtd.table(x, weights)$sum.of.weights[1:g])/sum(wtd.table(x, weights)$sum.of.weights))>.5)][1]
  if(is.null(weight))
    x[q] <- names(table(x))[eval(sapply(1:length(table(x)), function(g) sum(table(x)[1:g])/sum(table(x)))>.5)][1]
  x[r] <- names(sort(table(x), decreasing=TRUE))[1]
  x
}

na.levs.numeric <- function(x, naset=NULL, setmid=NULL, set1=NULL, set0=NULL, setmean=NULL, weight=NULL){
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

na.levs.default <- function(x, naset=NULL, setmid=NULL, set1=NULL, set0=NULL, setmean=NULL, weight=NULL){
  x <- as.numeric(x)
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
