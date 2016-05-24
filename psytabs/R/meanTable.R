meanTable <-
  function (value, ver.factor, hor.factor = NA, variation="se") {
  #Todo:
  # - mean Table for unequal N. Seite 23.
  if(is.na(hor.factor)) {
    if(variation=="sd") { 
      tab <- stats::aggregate(value ~ ver.factor, FUN=function(x) c(mn=mean(x), sd=stats::sd(x)))
      myFun <- function(x) {paste(round(x$mn, 2), round(x$sd, 2), sep="+-")}
    } else {
      tab <- stats::aggregate(value ~ ver.factor, FUN=function(x) c(mn=mean(x), se=standarderror(x)))
      myFun <- function(x) {paste(round(x$mn, 2), round(x$se, 2), sep="+-")}  
    }
    
    myFun2 <- function(x) {data.frame(x[,2])}  
    tab.splitted <- split(tab, tab[,1])
    tab.selected <- lapply(tab.splitted, myFun2)
    tab <- data.frame(t(data.frame(lapply(tab.selected, myFun))))
    rownames(tab) <- levels(ver.factor)
    colnames(tab) <- "Value"
    return(tab)
  }
    
    
  if(variation=="sd") { 
    tab <- stats::aggregate(value ~ ver.factor + hor.factor, FUN=function(x) c(mn=mean(x), sd=stats::sd(x)))
    myFun <- function(x) {paste(round(x$mn, 2), round(x$sd, 2), sep="+-")}
  } else {
    tab <- stats::aggregate(value ~ ver.factor + hor.factor, FUN=function(x) c(mn=mean(x), se=standarderror(x)))
    myFun <- function(x) {paste(round(x$mn, 2), round(x$se, 2), sep="+-")}  
  }
  
  myFun2 <- function(x) {data.frame(x[,3])}  
  tab.splitted <- split(tab, tab[,2])
  tab.selected <- lapply(tab.splitted, myFun2)
  tab <- data.frame(lapply(tab.selected, myFun))
  rownames(tab) <- levels(ver.factor)
  tab
}