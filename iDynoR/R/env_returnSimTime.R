env_returnSimTime <-
function(xmlResultFile)
{
	# xmlResultData[[1]] gives all the simulation tag information
	# Get the attributes of this tag, and return the time (the second attribute)
	return(as.numeric(xmlAttrs(xmlResultFile[[1]])[[2]]))
}
