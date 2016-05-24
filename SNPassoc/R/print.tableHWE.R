`print.tableHWE` <-
function(x, digits=4, sig=0.05, na="-", ...) 
 {
  x<-round(x,digits)
  x<-data.frame(x)
  if (ncol(x)<3)
   {
    names(x)[1]<-"HWE (p value)" 
    x$flag<- apply(x,1,function(x)ifelse(any(x<sig) & !is.na(x),"<-",""))
   }
  print(as.matrix(x), na.print=na, quote=FALSE, ...)
}

