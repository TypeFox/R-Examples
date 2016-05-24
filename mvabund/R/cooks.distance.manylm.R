############################################################################
# return Cook's distance of a manylm object			        	
# h_ii*r_i^2/(1-hii)^2 divided by k*s^2     	
############################################################################

cooks.distance.manylm <- function (model, 
    res = as.matrix(residuals(model)),
    sd = sqrt(deviance(model)/df.residual(model)), 
    hat = diag(model$hat.X), ...) 
{
        if (!is.null(model$weighted.residuals))
            res <- as.matrix(model$weighted.residuals)
	p  <- model$rank
	sd <- matrix(rep.int(sd, times=nrow(res)), nrow=nrow(res),
          ncol=length(sd), byrow=TRUE)
	dighat <- diag(1- hat)
	res    <- (((res/sd) / (1 - hat))^2 * hat)/p
  res[is.infinite(res)] <- NaN
  res
    
}

