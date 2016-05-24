## print.Monthmean.R
## Prints basic results from monthmean
## Oct 2009
print.Monthmean<-function(x, digits=1, ...){
## Check
  if (class(x) != "Monthmean"){
    stop("Object must be of class 'Monthmean'")
  } 
## Print
  toprint<-as.data.frame(cbind(month.name,round(x$mean,digits)))
  names(toprint)<-c('Month','Mean')
  print(toprint,row.names=F, ...)
} # end of function
