print.ipcwcompetingrisksSeSpPPVNPV <- function(x,No.lines=5,digits=2,...){
  # {{{ generate table to print
  ## browser()
  if(x$iid==TRUE){
    tab_ou_print<-round(cbind(x$Stats[,1:4],
                              x$TP*100,x$inference$vect_se_Se*100,
                              (1-x$FP_1)*100,x$inference$vect_se_Sp1*100,
                              (1-x$FP_2)*100,x$inference$vect_se_Sp2*100,
                              x$PPV*100,x$inference$vect_se_PPV*100,
                              x$NPV_1*100,x$inference$vect_se_NPV1*100,
                              x$NPV_2*100,x$inference$vect_se_NPV2*100
                              ),digits=digits)
    colnames(tab_ou_print)<-c("Cases","Survivors","Other events","Censored",
                              "Se","se_Se","Sp_1","se_Sp1","Sp_2","se_Sp_2",
                              "PPV (%)","se_PPV","NPV_1","se_NPV_1","NPV_2","se_NPV_2")
  }
  else{
    tab_ou_print<-round(cbind(x$Stats[,1:4],
                              x$TP*100,
                              (1-x$FP_1)*100,
                              (1-x$FP_2)*100,
                              x$PPV*100,
                              x$NPV_1*100,
                              x$NPV_2*100
                              ),digits=digits)
    colnames(tab_ou_print)<-c("Cases","Survivors","Other events","Censored",
                              "Se","Sp_1","Sp_2",
                              "PPV (%)","NPV_1","NPV_2"
                              )
  }
  # }}}
  cat(paste("Predictive accuracy measures at cutpoint c=",x$cutpoint," estimated using IPCW (n=",x$n, ", with competing risks). \n",sep=""))
  cat(paste("No. of positive (X>c) =",x$Stats[1,5],", \n",sep=""))
  cat(paste("No. of negative (X<=c) =",x$Stats[1,6],". \n",sep=""))
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
