env_returnSimIteration <-
function(xmlResultFile)
{
	# xmlResultData[[1]] gives all the simulation tag information
	# Get the attributes of this tag, and return the iteration (the first attribute)
	return(as.numeric(xmlAttrs(xmlResultFile[[1]])[[1]]))
}
