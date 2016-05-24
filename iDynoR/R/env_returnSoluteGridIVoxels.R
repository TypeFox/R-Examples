env_returnSoluteGridIVoxels <-
function(xmlResultFile,soluteRequested)
{
	return(as.numeric(xmlAttrs(xmlResultFile[[1]][[soluteRequested+1]])[4][[1]]))
}
