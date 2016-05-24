checkCytoscapeRunning <-
function(port.number=1234)
{
	base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")
	version.url = paste(base.url, "version", sep="/")
	if (requireNamespace("RCurl",quietly=TRUE)) {if (!RCurl::url.exists(version.url)){stop("Cytoscape with cyREST plugin is not running or you entered a wrong port number.\nPlease see the cyREST help page here: http://rpubs.com/keiono/24510")}} else {stop("RCurl package must be installed to use this function")}
	if (requireNamespace("httr",quietly=TRUE)) {cytoscape.version = httr::GET(version.url)} else {stop("httr package must be installed to use this function")}
	if (requireNamespace("RJSONIO",quietly=TRUE)) {cy.version = RJSONIO::fromJSON(rawToChar(cytoscape.version$content))} else {stop("RJSONIO package must be installed to use this function")}
	return(TRUE)
}
