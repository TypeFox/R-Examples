plotGraphL<-function(graph,margin=0,vertex.label.cex=0.6,vertex.label.font=1,vertex.size=8,
  vertex.size2=6,edge.arrow.size=0.2,edge.arrow.width=3,vertex.label=V(graph)$graphics_name,
  vertex.shape=V(graph)$graphics_type,layout=getLayout(graph),vertex.label.color="black",
  vertex.color=V(graph)$graphics_bgcolor,vertex.frame.color="dimgray",edge.color="dimgray",
  edge.label=getEdgeLabel(graph),edge.label.cex=0.6,edge.label.color="dimgray",edge.lty=getEdgeLty(graph),
  axes=FALSE,xlab="",ylab="",sub=NULL,main=NULL,...){
    if(class(graph)!="igraph") stop("the graph should be a igraph graph.")
    if(vcount(graph)==0){
         print("the graph is an empty graph.")
    }else{	 
    vertex.shape<-replace(vertex.shape,which(vertex.shape %in% c("roundrectangle","line")),"crectangle")
	vertex.shape<-ifelse(vertex.shape=="sphere","circle","rectangle")
	tempList<-get.edgelist(graph)[,1]
	loca<-length(tempList[grep("[A-Z]",tempList)])
    vertex.color<-replace(vertex.color,which(vertex.color %in% c("unknow","none")),"white")
    if(length(vertex.shape)==0) vertex.shape<-NULL
    if(length(vertex.color)==0) vertex.color<-NULL  
    if(length(vertex.label)==0) vertex.label<-NULL 
    if(length(layout)==0) layout<-NULL 
    if(length(edge.label)==0) edge.label<-NULL
    if((axes==FALSE)&&xlab==""&&ylab==""&&is.null(sub)&&is.null(main)){
         old.mai<-par(mai=c(0.01,0.25,0.01,0.3))
         #old.mai<-par(mai=0.01+c(0,0,0,0))
         on.exit(par(mai=old.mai), add=TRUE)
    }
    plot(graph,margin=margin,vertex.label.cex=vertex.label.cex,vertex.label.font=vertex.label.font,
	      vertex.size=vertex.size,vertex.size2=vertex.size2,
         edge.arrow.size=edge.arrow.size,edge.arrow.width=edge.arrow.width,vertex.label=vertex.label,
         vertex.shape=vertex.shape,layout=layout,vertex.label.color=vertex.label.color,
         vertex.color=vertex.color,vertex.frame.color=vertex.frame.color,edge.color=edge.color,
		 edge.label=edge.label,edge.label.cex=edge.label.cex,edge.label.color=edge.label.color,
		 edge.lty=edge.lty,axes=axes,xlab=xlab,ylab=ylab,sub=sub,main=main,...)
    }
	   
}
