summary.mritc <- function(object, ...){
    n <- nrow(object$prob)
    k <- ncol(object$prob)
    class <- max.col(object$prob)
    prop <- sapply(1:k, function(i) sum(class==i)/n)
    prop <- round(prop , 2) * 100
    prop[k] <- 100 - sum(prop[1:(k-1)])
    if(k==5){
        prop.new <- prop
        prop.new[2] <- prop[3]
        prop.new[3] <- prop[5]
        prop.new[4] <- prop[2]
        prop.new[5] <- prop[4]
        prop=prop.new
    }

    result <- matrix(0, nrow=k, ncol=3)
    result[,1] <- prop
    result[,2][1:3] <- round(object$mu, 2)
    result[,3][1:3] <- round(object$sigma, 2)
    if(k==3){
        rownames(result) <- c("CSF", "GM", "WM")
    }
    else{
        result[,2][4:k] <- NA
        result[,3][4:k] <- NA
        if(k==5){
            rownames(result) <- c("CSF", "GM", "WM", "CG", "GW")
        }
        else{
            rownames(result) <- NULL
        }
    }
    colnames(result) <- c(paste("proportion", "%", sep=""), "mean", "standard deviation")
    cat(paste("Using method ", object$method, ", the results are: \n", sep=""))
    print(result)
}
