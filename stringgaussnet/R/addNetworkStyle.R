addNetworkStyle <-
function (style.name,style.class,Annotations=F,points.size.map="PValue",min.points.value=0.05,max.points.value=0,points.fill.map="FC",min.points.fill=-2,max.points.fill=2,port.number=1234)
{
	check<-checkCytoscapeRunning(port.number)
	base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")
	
	# Defaults
	defaults <- list()
	defaults <- addSkeletonDefaults(defaults,Annotations)
	
	# Mappings
	mappings<-list()
	mappings <- eval(parse(text=paste("add",style.class,"Mappings(mappings)",sep="")))
	mappings<-addSkeletonMappings(mappings,Annotations,points.size.map,min.points.value,max.points.value,points.fill.map,min.points.fill,max.points.fill)
	
	style <- list(title=style.name, defaults = defaults, mappings = mappings)
	if (requireNamespace("RJSONIO",quietly=TRUE)) {style.JSON <- RJSONIO::toJSON(style)} else {stop("RJSONIO package must be installed to use this function")}
	
	style.url = paste(base.url, "styles", sep="/")
	if (requireNamespace("httr",quietly=TRUE)) {res<-httr::POST(url=style.url, body=style.JSON, encode = "json")} else {stop("httr package must be installed to use this function")}
}
