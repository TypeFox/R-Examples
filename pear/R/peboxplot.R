`peboxplot` <-
function(z, ...)
{
	ztitle <- attr(z, "title")
	z.names <- attr(z, "period.abb")
	if(is.null(z.names))
		z.names <- as.character(unique(cycle(z)))
	boxplot(split(z, cycle(z)), names = z.names, ...)
	if(!is.null(ztitle))
		title(main = ztitle)
}

