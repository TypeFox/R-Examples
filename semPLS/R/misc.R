# for object of class: plsm
exogenous <- function(model){
    if(!inherits(model, "plsm")) stop("Model must be of class 'plsm'!")
    ret <- names(which(colSums(model$D)==0))
    return(ret)
}

endogenous <- function(model){
    if(!inherits(model, "plsm")) stop("Model must be of class 'plsm'!")
    ret <- names(which(colSums(model$D)!=0))
    return(ret)
}

formative <- function(model){
    if(!inherits(model, "plsm")) stop("Model must be of class 'plsm'!")
    ret <- names(which(lapply(model$blocks, function(x){attr(x, "mode")})=="B"))
    return(ret)
}

reflective <- function(model){
    if(!inherits(model, "plsm")) stop("Model must be of class 'plsm'!")
    ret <- names(which(lapply(model$blocks, function(x){attr(x, "mode")})=="A"))
    return(ret)
}

indicators <- function(model, LV){
    if(!inherits(model, "plsm")) stop("Model must be of class 'plsm'!")
    if(!LV %in% model$latent) stop("The LV must be contained in the model!")
    ret <- model$blocks[[LV]]
    return(ret)
}

# used in 'pathWeighting'
predecessors <- function(model){
    if(!inherits(model, "plsm")) stop("Model must inherit from class 'plsm'!")
    D <- model$D
    foo <- function(x) names(which(x==1))
    pred <- apply(D, 2, foo)
    return(pred)
}

successors <- function(model){
    if(!inherits(model, "plsm")) stop("Model must inherit from class 'plsm'!")
    D <- model$D
    foo <- function(x) names(which(x==1))
    succ <- apply(D, 1, foo)
    return(succ)
}
