farforecast <- function (object, h = 10, var_type = "const", level = 80, PI = FALSE) 
{
  	x = as.numeric(rownames(object$y$y))
  	order = ncol(object$basis) - 1
  	if(requireNamespace("vars", quietly = TRUE)) 
  	{
    	var_pred = predict(vars::VAR(object$coeff[, 2:(order + 1)], lag.max = 5, type = var_type), 
    						n.ahead = h, ci = level/100)
  	}
  	else 
  	{
    	stop("Please install vars")
  	}
  	qconf <- qnorm(0.5 + level/200)
  	meanfcast <- varfcast <- matrix(NA, nrow = h, ncol = order)
  	for(i in 1:order) 
  	{
    	var_fit_pred = var_pred$fcst[[i]]
    	meanfcast[, i] = var_fit_pred[, 1]
    	varfcast[, i] = ((var_fit_pred[, 3] - var_fit_pred[, 2])/(2 * qconf))^2
  	}
  	point_fore = object$basis[, 2:(order + 1)] %*% t(meanfcast) + object$basis[, 1]
  	colnames(point_fore) = 1:h
  	point_fore_fts = fts(x, point_fore, yname = "Forecasts", xname = object$y$xname)
  	if(PI == TRUE) 
  	{
    	n.curve = ncol(object$y$y)
    	L = max(round(n.curve/5), order)
    	insample_fore = matrix(, nrow(object$y$y), (ncol(object$y$y) - L))
    	for(i in 1:(ncol(object$y$y) - L)) 
    	{
	      	dum = ftsm(fts(object$y$x, object$y$y[, 1:(L + i - 1)]), order = order)
      		var_pred = predict(vars::VAR(dum$coeff[, 2:(order + 1)], lag.max = 5, type = var_type), n.ahead = 1)
      		meanfcast = matrix(NA, nrow = 1, ncol = order)
	      	for(j in 1:order) 
	      	{
        		var_fit_pred = var_pred$fcst[[j]]
        		meanfcast[, j] = var_fit_pred[, 1]
      		}
      		insample_fore[, i] = dum$basis[, 2:(order + 1)] %*% t(meanfcast) + dum$basis[, 1]
    	}
    	insample_test = object$y$y[, (L + 1):ncol(object$y$y)]
    	resi = insample_test - insample_fore
    	lb_resi = apply(resi, 1, quantile, (100 - level)/200, na.rm = TRUE)
    	ub_resi = apply(resi, 1, quantile, (100 + level)/200, na.rm = TRUE)
    	lb = point_fore + lb_resi
    	ub = point_fore + ub_resi
	    colnames(lb) = colnames(ub) = 1:h
    	PI_lb = fts(x, lb, yname = "Lower bound", xname = object$y$xname)
    	PI_ub = fts(x, ub, yname = "Upper bound", xname = object$y$xname)
	    return(list(point_fore = point_fore_fts, PI_lb = PI_lb, PI_ub = PI_ub))
  	}
  	else 
  	{
    	return(point_fore_fts)
  	}
}
