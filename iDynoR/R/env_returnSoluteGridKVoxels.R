env_returnSoluteGridKVoxels <-
function(xmlResultFile,soluteRequested)
{
	return(as.numeric(xmlAttrs(xmlResultFile[[1]][[soluteRequested+1]])[6][[1]]))
}
