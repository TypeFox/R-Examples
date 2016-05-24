or.relimp.lm <- function (model, ui, ci = NULL, index = 2:length(coef(model)), 
    meq = 0, tol = sqrt(.Machine$double.eps), ...) 
{
    ## check input model
    if (!("lm" %in% class(model))) 
        stop("ERROR: model must be of class lm.")
    if (length(model$xlevels)>0) stop("model must not contain any factors!")
    if (max(attr(model$terms,"order")) != 1) 
        stop ("model must not contain higher order terms")
    
    namen <- names(coef(model))
    
    ## work is done by functions all.R2 from this package 
    ## and function Shapley.value from package kappalab
    ## output is currently very limited
    
    ## prepare data for calculation of sub models
    resp <- attr(model$terms, "response")
    xcol <- which(rowSums(attr(model$terms, "factors")) > 0)
    DATA <- as.data.frame(model$model[, c(resp, xcol)])
    wt <- weights(model)
    if (is.null(wt)) wt <- rep(1/nrow(DATA),nrow(DATA))
    
    aus <- Shapley.value(set.func(all.R2(cov.wt(DATA,wt=wt)$cov, ui, ci = ci, index = index, 
         meq = meq, tol = tol, ...)))
    names(aus) <- namen[-1]
    aus
}