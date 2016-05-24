applyStyle <-
function(style.name,network.suid,port.number=1234)
{
	check<-checkCytoscapeRunning(port.number)
	base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")
	apply.style.url = paste(base.url, "apply/styles", style.name, toString(network.suid), sep="/")
	if (requireNamespace("httr",quietly=TRUE)) {res <- httr::GET(apply.style.url)} else {stop("httr package must be installed to use this function")}
}
