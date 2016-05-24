##' @export 
print.summary.prodlim <- function(x,digits=ifelse(x$percent,1,3),...){
    model <- x$model
    cotype <- x$cotype
    sumtable <- x$table
    if (x$asMatrix==TRUE){
        print(sumtable,digits=digits,quote=FALSE,...)
    }
    else{
        if (model=="survival"){
            if (cotype==1){
                print(sumtable,digits=digits,quote=FALSE,...)
            } else{
                  print.listof(sumtable,digits=digits,quote=FALSE,...)
              }
        } else{
              if (model=="competing.risks"){
                  for (cc in 1:length(sumtable)){
                      cat("\n\n----------> Cause: ",names(sumtable)[cc],"\n\n")
                      if (cotype==1){
                          print(sumtable[[cc]],digits=digits,quote=FALSE,...)
                      }
                      else{
                          print.listof(sumtable[[cc]],digits=digits,quote=FALSE,...)
                      }
                  }
              }
          }
    }
}
