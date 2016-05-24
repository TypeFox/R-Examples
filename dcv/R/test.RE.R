test.RE <-
function(x, y){
    if(!length(x) == length(y)){
        stop("\"x\" and \"y\" have different length.")
    }
    sum <- 0
    sum1 <- 0
    sum2 <- 0
    n <- length(x)
    for(i in 1:n){
        sum <- sum + (x[i] - y[i])^2  
        sum1 <- sum1 + x[i]^2
        sum2 <- sum2 + (x[i] - mean(x))^2
    }
    #Defined by Wu Xiangding
    #RE <- 1-(sum/sum1)
    
    #Defined by foreign RE
    RE <- 1-(sum/sum2) 
    
    #Mean squared erro of validation
    MSE <- sum/n 
    
    # Root mean squared error of validation
    RMSE <- sqrt(sum/n) 
    res <- list(RE, MSE, RMSE)
    class(res) <- "RE"
    return(res)
}

