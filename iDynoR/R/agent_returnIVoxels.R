agent_returnIVoxels <-
function(xmlResultData)
{
	# xmlResultData[[1]] gives all the simulation tag information
	# The grid information is at the [[1]] element of this
	# Get the attributes of this tag, and return nI (the second attribute)
	return(as.numeric(xmlAttrs(xmlResultData[[1]][[1]])[[2]]))
}
