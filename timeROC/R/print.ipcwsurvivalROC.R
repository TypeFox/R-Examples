print.ipcwsurvivalROC <- function(x,No.lines=5,digits=2,...){ 
  # {{{ generate table to print
  if(x$iid==TRUE){
    tab_ou_print<-round(cbind(x$Stats,x$AUC*100,x$inference$vect_sd_1*100),digits=digits)
    colnames(tab_ou_print)<-c("Cases","Survivors","Censored","AUC (%)","se")
  }
  else{
    tab_ou_print<-round(cbind(x$Stats,x$AUC*100),digits=digits)
    colnames(tab_ou_print)<-c("Cases","Survivors","Censored","AUC (%)")
  }
  # }}}
  cat(paste("Time-dependent-Roc curve estimated using IPCW  (n=",x$n, ", without competing risks). \n",sep=""))
  # {{{ print the table or the quantile rows of the table
  l<-length(x$times)
  if(l<=No.lines){  print(tab_ou_print) }
  else{print(tab_ou_print[unique(round(quantile(1:length(x$times),probs=seq(0,1,length.out=No.lines)),0)),])}
  # }}}
  cat("\n")
  cat("Method used for estimating IPCW:")
  cat(paste(x$weights$method,"\n"))
  cat("\n")
  cat("Total computation time :",round(x$computation_time,2)," secs.")
  cat("\n")
}
