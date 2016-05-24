`m4plNoMoreSummary` <-
function(x){
 stat1 <- sapply(x, mean, na.rm=TRUE) ##
 stat2 <- sapply(x, sd,   na.rm=TRUE) ##
 res <- data.frame(t(data.frame(mean=stat1,sd=stat2)))
 return(res)
 }
