or.relimp.default <- function(model, ui, ci = NULL, index = 2:ncol(model), 
    meq = 0, tol = sqrt(.Machine$double.eps), ...) 
{
    ## check input model
    if (!(is.matrix(model))) 
        stop("ERROR: model must be of class lm or a covariance matrix.")
    else
    if (!(nrow(model)==ncol(model)))
        stop("ERROR: If it is not a linear model, model must be a quadratic matrix.")
    else 
    if (!(all(eigen(model,TRUE,only.values=TRUE)$values>0)))
        stop("ERROR: matrix model must be positive definite.")
    
    namen <- colnames(model)
    if (is.null(namen)) namen <- c("y",paste("X",1:(ncol(model)-1),sep=""))
    
    ## work is done by functions all.R2 from this package 
    ## and function Shapley.value from package kappalab
    ## output is currently very limited
    
    ## prepare data for calculation of sub models
    aus <- Shapley.value(set.func(all.R2(model, ui, ci = ci, index = index, 
         meq = meq, tol = tol, ...)))
    names(aus) <- namen[-1]
    aus
}