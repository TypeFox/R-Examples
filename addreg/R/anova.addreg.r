anova.addreg <- function(object, ..., test = NULL)
{
	dotargs <- list(...)
	named <- if(is.null(names(dotargs))) rep(FALSE, length(dotargs)) else(names(dotargs) != "")
	if(any(named))
		warning("the following arguments to 'anova.addreg' are invalid and dropped: ",
			paste(deparse(dotargs[named]), collapse = ", "))
	dotargs <- dotargs[!named]
	is.addreg <- unlist(lapply(dotargs, function(x) inherits(x, "addreg")))
	dotargs <- dotargs[is.addreg]
	if(length(dotargs))
		return(anova.addreglist(c(list(object), dotargs), test = test))
	else
		stop('anova.addreg does not support anova for a single model. Fit nested models manually and input to anova.addreg')
}