median.factor <-
function (x,na.rm=TRUE) {
  if (!is.ordered(x)) {stop(paste0(deparse(substitute(x))," is not ordered"))}
  med.num <- median(as.numeric(x),na.rm=na.rm)
  if ((med.num/0.5)%%2==1) {med.num <- c(floor(med.num),ceiling(med.num))}
  med.lev <- paste(levels(x)[med.num],collapse=",")
  med.lev
}
