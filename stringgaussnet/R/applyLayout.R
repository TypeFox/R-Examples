applyLayout <-
function(network.suid,layout.name,port.number=1234)
{
	check<-checkCytoscapeRunning(port.number)
	base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")
	apply.layout.url = paste(base.url, "apply/layouts",layout.name, toString(network.suid), sep="/")
	if (requireNamespace("httr",quietly=TRUE)) {res <- httr::GET(apply.layout.url)} else {stop("httr package must be installed to use this function")}
}
