## Model a nls class model
## This is an auxiliar function that returns a function 
## for estimating the derivatives in the loq_derivatives function
getPredict <- function(model) {
    stopifnot(inherits(model,"nls"))
    lhs <- as.character(model$m$formula()[[2]])
    parameters <- names(model$m$getPars())
    allobj <- ls(model$m$getEnv())
    rhs <- allobj[-match(c(parameters,lhs),allobj)]
    fx <- function(x) {
        newlist <- list(x)
        names(newlist) <- rhs
        predict(model,newdata=newlist)
    }
    assign("model", model, envir=environment(fx))
    assign("rhs", rhs, envir=environment(fx))
    return(fx)  
}
