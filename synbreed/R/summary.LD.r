# summary for LD objects

summary.LDdf <- function(object,...){
     LDdf <- object
     # help functions
     noM <- function(x) length(unique(x$marker1))+1
     avgr2 <- function(x) mean(x$r2,na.rm=TRUE)
     minr2 <- function(x) min(x$r2,na.rm=TRUE)
     maxr2 <- function(x) max(x$r2,na.rm=TRUE)
     pr <- function(x) mean(x$r2>0.2,na.rm=TRUE)
     maxd <- function(x) max(x$dist,na.rm=TRUE)

     ret <- data.frame(noM=unlist(lapply(LDdf,noM)),avgr2= unlist(lapply(LDdf,avgr2)),minr2= unlist(lapply(LDdf,minr2)),maxr2= unlist(lapply(LDdf,maxr2)),Pr02=unlist(lapply(LDdf,pr)),maxDist = unlist(lapply(LDdf,maxd)))
     return(ret)   
}

summary.LDmat <- function(object,...){
     LDmat <- object
     # help functions
     avgr2 <- function(x) mean(x[upper.tri(x)],na.rm=TRUE)
     minr2 <- function(x) min(x[upper.tri(x)],na.rm=TRUE)
     maxr2 <- function(x) max(x[upper.tri(x)],na.rm=TRUE)
     pr <- function(x) mean(x[upper.tri(x)]>0.2,na.rm=TRUE)

     
     ret <- data.frame(noM=unlist(lapply(LDmat$LD,ncol)),avgr2= unlist(lapply(LDmat$LD,avgr2)),minr2= unlist(lapply(LDmat$LD,minr2)),maxr2= unlist(lapply(LDmat$LD,maxr2)),Pr02=unlist(lapply(LDmat$LD,pr)),maxDist = unlist(lapply(LDmat$distance,avgr2)))
     return(ret)   
}

