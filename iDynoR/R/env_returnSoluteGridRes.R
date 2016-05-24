env_returnSoluteGridRes <-
function(xmlResultFile,soluteRequested)
{
	return(as.numeric(xmlAttrs(xmlResultFile[[1]][[soluteRequested+1]])[3][[1]]))
}
