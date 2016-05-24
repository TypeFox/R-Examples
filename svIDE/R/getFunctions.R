getFunctions <- function (pos)
{
	## Get a list of all R functions in a certain position
	return(as.character(lsf.str(pos = pos, all.names = TRUE)))
}
