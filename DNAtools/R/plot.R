plot.dbcompare <- function(x,log="y",las=3,xlab="Match/Partial",ylab="Counts",...){
  nl <- attributes(x)$call$loci
  levs <- dbCats(nl,vector=TRUE)
  if(is.matrix(x$m)) mvec <- t(x$m)[up.tri(x$m)]
  else mvec <- x$m
  mvec[mvec==0] <- NA
  if(attributes(x)$call$collapse){
    mcol <- ifelse(x$m==0,NA,x$m)
    if(xlab=="Match/Partial") xlab <- "Total number of matching alleles"
    plot(0:(length(mcol)-1),mcol,log=log,xlab=xlab,ylab=ylab,...)
  }
  else{
    plot(1:length(levs),mvec,axes=FALSE,xlab=xlab,ylab=ylab,log=log,...)
    axis(1,at=1:length(levs),labels=levs,las=las)
    axis(2)
    box()
  }
}
