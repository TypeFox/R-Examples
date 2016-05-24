print.MgSplitMethods <- function(x){
    if (x$name!="full data"){
        cat("Cross-Validation splitMethod:",x$name,"\n")
        if (x$internal.name=="crossval")
            cat("Repeat: ",x$B,"\n")
        else{
            if (x$M<x$N)
                cat("Type: subsampling\n",x$M," out of ",x$N,"\n\n")
            cat("No. bootstrap samples: ",x$B,"\n")
        }
    }
    cat("Sample size: ",x$N,"\n")
}
