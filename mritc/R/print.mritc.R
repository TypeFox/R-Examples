print.mritc <- function(x, ...){
    n <- nrow(x$prob)
    k <- ncol(x$prob)
    class <- max.col(x$prob)
    prop <- sapply(1:k, function(i) sum(class==i)/n)
    prop <- round(prop , 2) * 100
    prop[k] <- 100 - sum(prop[1:(k-1)])
    if(k==3){
        cat(paste("Using method ", x$method, ", among ", n,
                  " voxles, \n", prop[1], "% are classified to CSF,\n",
                  prop[2], "% are classified to GM,\n", 
                  prop[3], "% are classified to WM.\n", sep=""))
    }
    else{
        if(k==5){
            cat(paste("Using method ", x$method, ", among ", n,
                      " voxles, \n", prop[1], "% are classified to CSF,\n",
                      prop[2], "% are classified to CG,\n",
                      prop[3], "% are classified to GM,\n",
                      prop[4], "% are classified to GW,\n",
                      prop[5], "% are classified to WM.\n", sep=""))

        }
        else{
            cat(paste("Using method ", x$method, ", among ", n,
                      " voxles, \n",
                      "The proportions of different tissue types are: \n",
                      prop, sep=""))
        }
    }
}
