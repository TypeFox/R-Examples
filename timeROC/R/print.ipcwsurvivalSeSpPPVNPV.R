print.ipcwsurvivalSeSpPPVNPV <- function(x,No.lines=5,digits=2,...){
  # {{{ generate table to print
  ## browser()
  if(x$iid==TRUE){
    tab_ou_print <- round(cbind(x$Stats[,1:3],
                                x$TP*100,x$inference$vect_se_Se*100,
                                (1-x$FP)*100,x$inference$vect_se_Sp1*100,
                                x$PPV*100,x$inference$vect_se_PPV*100,
                                x$NPV*100,x$inference$vect_se_NPV1*100
                                ),digits=digits)
    colnames(tab_ou_print) <- c("Cases","Survivors","Censored",
                                "Se (%)","se_Se","Sp (%)","se_Sp",
                                "PPV (%)","se_PPV","NPV (%)","se_NPV")
  }
  else{
    tab_ou_print <- round(cbind(x$Stats[,1:3],
                                x$TP*100,
                                (1-x$FP)*100,
                                x$PPV*100,
                                x$NPV*100
                                ),digits=digits)
    colnames(tab_ou_print) <- c("Cases","Survivors","Censored",
                                "Se (%)","Sp (%)",
                                "PPV (%)","NPV (%)")
  }
  # }}}
  cat(paste("Predictive accuracy measures at cutpoint c=",x$cutpoint," estimated using IPCW (n=",x$n, ", with competing risks). \n",sep=""))
  cat(paste("No. of positive (X>c) =",x$Stats[1,4],", No. of negative (X<=c) =",x$Stats[1,5],". \n",sep=""))
  cat("\n")
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
