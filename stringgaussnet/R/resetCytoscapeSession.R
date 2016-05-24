resetCytoscapeSession <-
function(port.number=1234)
{
	check<-checkCytoscapeRunning(port.number)
	base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")
	reset.url <- paste(base.url,"session",sep="/")
	if (requireNamespace("httr",quietly=TRUE)) {res<-httr::DELETE(reset.url)} else {stop("httr package must be installed to use this function")}
}
