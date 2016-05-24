env_returnSoluteGridJVoxels <-
function(xmlResultFile,soluteRequested)
{
	return(as.numeric(xmlAttrs(xmlResultFile[[1]][[soluteRequested+1]])[5][[1]]))
}
