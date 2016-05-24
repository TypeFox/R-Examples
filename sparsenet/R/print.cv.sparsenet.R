print.cv.sparsenet=function(x,digits = max(3, getOption("digits") - 3),...){
    cat("\nCall: ", deparse(x$call), "\n\n")
    whichs=rbind(x$which.min,x$which.1se)[,c(2,1)]
    dimnames(whichs)=list(NULL,c("i","j"))
outtab=data.frame(whichs,Gamma=c(x$parms.min[1],x$parms.1se[1]),Lambda=c(x$parms.min[2],x$parms.1se[2]),"CV Error"=c(x$cvm[x$which.min[1],x$which.min[2]],x$cvm[x$which.1se[1],x$which.1se[2]]),"CV up"=c(x$cvup[x$which.min[1],x$which.min[2]],NA))
    row.names(outtab)=c("Min","1SE")
 print(outtab,na.print="" )
  }
