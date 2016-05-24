getNodeLabel<-function(graph,type="symbol",displayNumber=1){
     
     if(displayNumber<1) stop("displayNumber should >=1")
     graphics_name<-V(graph)$graphics_name
     node.names<-V(graph)$names
	 
     current.org<-unlist(strsplit(graph$org,";"))[1]
	 if(type=="symbol"){
	     if(current.org=="hsa"){
	         node.label<-lapply(node.names, function(x) getSymbolFromGene(unlist(strsplit(x,"[ ;]"))))
	     }else{stop("It is not ec, ko, or org graph.")}
	 }
	 node.label.new<-sapply(node.label,function(x) ifelse(length(x)>displayNumber,
	                  paste(paste(x[1:displayNumber],collapse=","),"...",sep=","),paste(x,collapse=",")))
	 for(i in seq(node.label.new)){
	    if(node.label.new[i]==""){
		   node.label.new[i]<-graphics_name[i]
		}
	 }
     return(node.label.new)
}