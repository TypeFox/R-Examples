# last modified 25 January 2009 by J. Fox

hetcor.default <-
	function (data, ..., ML = FALSE, std.err = TRUE, bins = 4, pd = TRUE) 
{
	dframe <- data.frame(data, ...)
	if (!missing(...)) names(dframe)[1] <- deparse(substitute(data))
	hetcor(dframe, ML = ML, std.err = std.err, bins = bins, pd = pd)
}

