
print.predLongi <- function(x, digits = 3, ...)
{
        if(class(x)!="predLongi"){
                stop("The object x must be a class predLongi.")
        }else{
                if (!is.null(cl <- x$call)){
                        cat("Call:\n")
                        dput(cl)
                        cat("\n")
                }

                        cat("\n")
                        cat("--------- Prediction given the exact history of longitudinal outcome ---------\n")
                        cat("\n")
                        print(x$pred,row.names=F,digits=digits)


        }
}
