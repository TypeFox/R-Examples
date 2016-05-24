saveCytoscapeSession <-
function (filepath="stringgaussnet_networks",overwrite=F,absolute=F,port.number=1234)
{
	check<-checkCytoscapeRunning(port.number)
	filepath<-paste(gsub("\\.cys$","",filepath),"cys",sep=".")
	if (file.exists(filepath) & !overwrite){stop(paste(filepath,"already exists."))}
	if(!absolute){filepath<-paste(getwd(),filepath,sep="/")}
	filepath<-URLencode(filepath)
	base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")
	save.url<-paste(base.url,paste("session?file=",filepath,sep=""),sep="/")
	if (requireNamespace("httr",quietly=TRUE)) {res <- httr::POST(save.url)} else {stop("httr package must be installed to use this function")}
}
