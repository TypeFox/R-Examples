boot.orlm <- function (model, B = 1000, fixed = FALSE, ui, ci, index, meq) 
{
    ## check for admissible model
    if (!("lm" %in% class(model))) 
        stop("ERROR: model must be of class lm.")
    ## prepare data for bootstrap sampling
    resp <- attr(model$terms, "response")
    xcol <- which(rowSums(attr(model$terms, "factors")) > 0)
    DATA <- as.data.frame(model$model[, c(resp, xcol)])
    wt <- weights(model)
    if (is.null(wt)) 
        wt <- rep(1/nrow(DATA), nrow(DATA))
    if (!fixed) 
        booterg <- boot(cbind(wt = wt, DATA), orlm.forboot, B, 
            ui = ui, ci = ci, index = index, meq = meq)
    else {
        e <- model$residuals
        fit <- model$fitted.values
        booterg <- boot(data.frame(DATA, fit = fit, e = e), orlm.forboot.fixed, 
            B, ui = ui, ci = ci, index = index, meq = meq)
    }
    booterg
}
