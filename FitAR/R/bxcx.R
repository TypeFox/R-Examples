`bxcx` <-
function(x, lambda, InverseQ=FALSE, type="BoxCox"){
if (type=="BoxCox"){    
    if (!InverseQ) {
        if (min(x) <= 0){
            cat("min data value <= 0, 0.25-min(x) added to data", fill=TRUE)
            x<- x + 0.25-min(x)
            }    
        if (abs(lambda)<1e-6)
            log(x)
        else
            (x^lambda-1)/lambda
        }
    else
        if (abs(lambda)<1e-6)
            exp(x)
        else {
            y <- lambda*x+1
            if (min(y) <= 0) 
                cat("Warning: modified inverse Box-Cox transformation used", fill=TRUE)
            abs(y)^(1/lambda)
            }
    }
else {    #simple power transformation
    if (!InverseQ) {
        if (min(x) <= 0){
            cat("min data value <= 0, 0.25-min(x) added to data", fill=TRUE)
            x<- x + 0.25-min(x)
            }    
        if (abs(lambda)<1e-6)
            log(x)
        else
            x^lambda
        }
    else
        if (abs(lambda)<1e-6)
            exp(x)
        else {
            y <- x
            if (min(y) <= 0) 
                cat("Warning: modified inverse Box-Cox transformation used", fill=TRUE)
            abs(y)^(1/lambda)
            }
    }
}

