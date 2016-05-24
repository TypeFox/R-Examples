PlotPathwayGraph<-function(graph,margin=0,vertex.label.cex=0.6,vertex.label.font=1,vertex.size=8,
  vertex.size2=6,vertex.shape="rectangle",layout=layout.random,vertex.label.color="black",edge.color="dimgray",
  vertex.color="#C1FFC1",vertex.frame.color="dimgray",axes=FALSE,xlab="",ylab="",sub=NULL,main=NULL,...){
    Split<-unlist(strsplit(as.character(graph[,2]),split="|",fixed=T))
    graph<-graph.data.frame(as.data.frame(cbind(first=Split[(1:(length(Split)/2))*2-1],sencond=Split[(1:(length(Split)/2))*2],graph[,c(4,6)]),stringsAsFactors=FALSE),directed=F)
    if(class(graph)!="igraph") stop("the graph should be a igraph graph.")
    if(vcount(graph)==0){
         print("the graph is an empty graph.")
    }else{	 
	       E(graph)$color<-ifelse(E(graph)$CORE_ENRICHMENT=='YES', 'red', 'dimgray')
		   if(length(vertex.shape)==0) vertex.shape<-NULL
           if(length(vertex.color)==0) vertex.color<-NULL  
           if(length(layout)==0) layout<-NULL 
           if((axes==FALSE)&&xlab==""&&ylab==""&&is.null(sub)&&is.null(main)){
          old.mai<-par(mai=c(0.01,0.25,0.01,0.3))
          #old.mai<-par(mai=0.01+c(0,0,0,0))
          on.exit(par(mai=old.mai), add=TRUE)
    }
		   plot(graph,margin=margin,vertex.label.cex=vertex.label.cex,vertex.label.font=vertex.label.font,vertex.size=vertex.size, vertex.size2=vertex.size2,vertex.shape=vertex.shape,layout=layout,vertex.label.color=vertex.label.color,edge.color=E(graph)$color,vertex.color=vertex.color,vertex.frame.color=vertex.frame.color)

  }
  }