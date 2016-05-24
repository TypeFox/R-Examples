addGraphToCytoscape <-
function (Network,Collection=class(Network),Name=deparse(substitute(Network)),LayoutName="force-directed",StyleName=Collection,port.number=1234)
{
	igraphobj<-as.igraph.stringgaussnet(Network)
	igraphobj$name<-Name
	
	cygraph<-toCytoscape(igraphobj)
	
	check<-checkCytoscapeRunning(port.number)
	base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")
	
	network.url = paste(base.url, paste("networks?collection=",Collection,sep=""), sep="/")
	if (requireNamespace("httr",quietly=TRUE)) {res <- httr::POST(url=network.url, body=cygraph, encode="json")} else {stop("httr package must be installed to use this function")}
	if (requireNamespace("RJSONIO",quietly=TRUE)) {network.suid = unname(RJSONIO::fromJSON(rawToChar(res$content)))} else {stop("RJSONIO package must be installed to use this function")}
	applyLayout(network.suid,LayoutName,port.number)
	applyStyle(StyleName,network.suid,port.number)
	
	return(network.suid)
}
