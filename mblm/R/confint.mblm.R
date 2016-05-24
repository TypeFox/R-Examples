"confint.mblm" <-
function (object, parm, level = 0.95, ...) 
{
	res = c(0,0,0,0); dim(res) = c(2,2);
	rownames(res) = names(object$coefficients)
	colnames(res) = as.character(c((1-level)/2,1-(1-level)/2))	
	res[2,] = wilcox.test(object$slopes,conf.int=TRUE,conf.level=level)$conf.int
	res[1,] = wilcox.test(object$intercepts,conf.int=TRUE,conf.level=level)$conf.int

res

}

