### print.reclassification.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct  3 2015 (16:26) 
## Version: 
## last-updated: Oct  3 2015 (16:26) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
print.riskReclassification <- function(x,percent=TRUE,digits=ifelse(percent,1,2),...){
    cat("Observed overall re-classification table:\n\n")
    print(x$reclassification)
    cat("\nExpected re-classification probabilities (%) among subjects with event until time ",x$time,"\n\n",sep="")
    fmt <- paste0("%1.", digits[[1]], "f")
    dnames <- dimnames(x$reclassification)
    dim <- dim(x$reclassification)
    if (percent==TRUE){
        rlist <- lapply(x$event.reclassification,function(x){
                            matrix(sprintf(fmt=fmt,100*c(x)),nrow=dim[1],ncol=dim[2],dimnames=dnames)
                        })
    }else{
         rlist <- lapply(x$event.reclassification,function(x){
                             matrix(sprintf(fmt=fmt,c(x)),nrow=dim[1],ncol=dim[2],dimnames=dnames)
                         })
     }
    if (x$model=="competing.risks"){
        print.listof(rlist[-length(rlist)],quote=FALSE)
    } else{
          print.listof(rlist[1],quote=FALSE)
      }
    cat("\nExpected re-classification probabilities (%) among subjects event-free until time ",x$time,"\n\n",sep="")
    print.listof(rlist[length(rlist)],quote=FALSE)
}
#----------------------------------------------------------------------
### print.reclassification.R ends here
