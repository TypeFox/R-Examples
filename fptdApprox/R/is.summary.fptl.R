is.summary.fptl <-
function (obj) 
{
    if (inherits(obj, "summary.fptl") & is.list(obj) & all(is.element(names(attributes(obj)), c("Call", "FPTLCall", "dp", "vars", "id", "class")))){
		if (all(sapply(obj, is.list)) & all(sapply(obj, length) == 2L) & all(sapply(lapply(obj, names), identical, c("instants", "FPTLValues"))) & 
		is.list(attr(obj, "Call")) & is.list(attr(obj, "FPTLCall")) & is.diffproc(attr(obj, "dp")) & (is.list(attr(obj, "vars")) | is.null(attr(obj, "vars"))) & 
		(((length(obj) == 1) & is.null(attr(obj, "id"))) | ((length(obj) > 1L) & is.list(attr(obj, "id")) & (length(attr(obj, "id")) == 4L)))){
			if (all(sapply(obj, function(x) is.numeric(x$instants) & is.matrix(x$instants) & (ncol(x$instants) == 5L))) & 
			all(sapply(obj, function(x) is.numeric(x$FPTLValues) & is.matrix(x$FPTLValues) & (ncol(x$FPTLValues) == 5L))) & 
			all(sapply(obj, function(x) nrow(x$instants) == nrow(x$FPTLValues))) & all(sapply(attr(obj, "Call"), is.call)) &
			all(sapply(attr(obj, "FPTLCall"), is.call))){
				if ((length(unique(lapply(attr(obj, "FPTLCall"), function(x) as.list(x)[-5]))) == 1L)) return(TRUE)		
			}
		}
	}
	return(FALSE)
}
