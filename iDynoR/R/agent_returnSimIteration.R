agent_returnSimIteration <-
function(xmlResultData)
{
	# xmlResultData[[1]] gives all the simulation tag information
	# Get the attributes of this tag, and return the iteration (the first attribute)
	return(as.numeric(xmlAttrs(xmlResultData[[1]])[[1]]))
}
