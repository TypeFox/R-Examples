
print.epoce <- function(x, digits = 3, ...)
{
        if(class(x)!="epoce"){
                stop("The object x must be a class epoce.")
        }else{
#               if (!is.null(cl <- x$call)){
#                       cat("Call:\n")
#                       dput(cl)
#                       cat("\n")
#               }
                if (x$new.data){
                        out <- matrix(x$mpol,nrow=1,ncol=length(x$pred.times),byrow=TRUE)
                        rownames(out) <- c("mpol")
                }else{
                        out <- matrix(c(x$mpol,x$cvpol),nrow=2,ncol=length(x$pred.times),byrow=TRUE)
                        rownames(out) <- c("mpol","cvpol")
                }
                colnames(out) <- round(x$pred.times,3)
                print(out)
        }
}

