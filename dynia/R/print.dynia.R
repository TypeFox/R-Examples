print.dynia<-function(x,...){
  printout<-list("delta"=x$delta,"Int.Mod"=x$Int.Mod,"p-value"=x$pvalue)
  print(printout)
}