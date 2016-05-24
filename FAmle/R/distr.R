distr <-
function(x,dist,param,type='d',model=NULL,...)
{
	if(!is.null(model))
	{
		param <- model$par.hat
		dist <- model$dist
	}
	if(is.array(param) & length(x) > 1)
	{
		input <- as.list(as.data.frame(cbind(x,param)))
		names(input) <- NULL
	}
	else if(is.array(param) & length(x) == 1)
	{
		input <- as.list(as.data.frame(param))
		n.temp <- length(input)
		temp <- input
		input <- list(x)
		for(i in 1:n.temp) input[[i+1]] <- temp[[i]]
		names(input) <- NULL		
	}
	else
	{
#		input <- as.list(c(NaN,param))
		input <- as.list(c(NaN,as.numeric(param)))
		input[[1]] <- x
	}
	n.arg <- length(input)
	extra <- list(...)
	n.extra <- length(extra)
	if(n.extra != 0) for(i in 1:n.extra) input[[names(extra)[i]]] <- extra[[i]]
	do.call(paste(type,dist,sep=''),input)
}

