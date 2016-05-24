predict20x = function(object, newdata, cilevel = 0.95, digit = 3,
                     print.out = TRUE, ...){
## prediction which allows for factors
# and a data frame with data entered in the same order as the data frame that was used
# in fitting the above model (note: the variable names do not need to be specified)

    if (!inherits(object, "lm"))
        stop("First input is not an \"lm\" object")
    
    if (!is.data.frame(newdata))
        stop("Argument \"newdata\" is not a data frame!")


  name.row = paste("pred",1:nrow(newdata),sep=".")
	name.row = 1:nrow(newdata)
        #name.col = attr(object$terms,"term.labels")
	x = attr(object$terms,"term.labels")

	y = unlist(strsplit(x,"factor\\("))
	z = unlist(strsplit(y,"\\)"))
	name.col = z
	#print(name.col)

    if (ncol(newdata) != length(name.col))
        stop("Incorrectly input the new data!")


        dimnames(newdata) = list(name.row,name.col)

    pred = predict.lm(object, newdata, se.fit = TRUE, ...)
    Predicted = pred$fit
    percent = 1 - (1 - cilevel)/2
    Conf.lower = pred$fit - qt(percent, pred$df) * pred$se.fit
    Conf.upper = pred$fit + qt(percent, pred$df) * pred$se.fit
    pred.se = sqrt(pred$residual.scale^2 + pred$se.fit^2)
    Pred.lower = pred$fit - qt(percent, pred$df) * pred.se
    Pred.upper = pred$fit + qt(percent, pred$df) * pred.se
    mat = cbind(Predicted, Conf.lower, Conf.upper, Pred.lower,
        Pred.upper)
    mat = round(mat, digit)
    mat.df = as.data.frame(mat)
    dimnames(mat.df)[[1]] = dimnames(newdata)[[1]]
    dimnames(mat.df)[[2]] = c("Predicted", " Conf.lower", "Conf.upper",
        " Pred.lower", " Pred.upper")
  
    if (print.out)
        print(mat.df)
  
    invisible(list(frame = mat.df, fit = pred$fit, se.fit = pred$se.fit,
        residual.scale = pred$residual.scale, df = pred$df, cilevel = cilevel))
}

