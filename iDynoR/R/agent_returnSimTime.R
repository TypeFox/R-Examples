agent_returnSimTime <-
function(xmlResultData)
{
	# xmlResultData[[1]] gives all the simulation tag information
	# Get the attributes of this tag, and return the time (the second attribute)
	return(as.numeric(xmlAttrs(xmlResultData[[1]])[[2]]))
}
