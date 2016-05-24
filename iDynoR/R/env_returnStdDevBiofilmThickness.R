env_returnStdDevBiofilmThickness <-
function(xmlResultFile)
{
	# xmlResultData[[1]] gives all the simulation tag information
	# xmlResultData[[1]][[1]] gives the simulation tag of the env_state file
	# xmlResultData[[1]][[1]][[1]] gives the thickness tag set
	# xmlResultData[[1]][[1]][[1]][[1]] gives the mean thickness value	
	return(xmlResultFile[[1]][[1]][[2]][[1]])
}
