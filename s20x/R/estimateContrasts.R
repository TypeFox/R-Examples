estimateContrasts <- function (contrast.matrix, fit, row = TRUE, alpha = 0.05,L=NULL) 
{
    if (!inherits(fit, "lm")) 
        stop("Second input is not an \"lm\" object")
 	if (length(dimnames(fit$model)[[2]]) == 2) 
		estimateContrasts1(contrast.matrix, fit, alpha = alpha,L)
	else estimateContrasts2(contrast.matrix, fit, alpha = alpha, row,L)
}

