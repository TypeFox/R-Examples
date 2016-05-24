write.graph.to.sonia <- function(source_graph,fileN)
{
	export_graph <- source_graph

	listVertexAtts <-list.vertex.attributes(export_graph)
	for(i in 1:length(listVertexAtts))
	{
		export_graph<-remove.vertex.attribute(export_graph,listVertexAtts[i])
	}

	listEdgeAtts <-list.edge.attributes(export_graph)
	for(i in 1:length(listEdgeAtts))
	{
		export_graph<-remove.edge.attribute(export_graph,listEdgeAtts[i])
	}
	
	#rebuild acceptable .son attributes
	maxVertex <-vcount(source_graph)-1
	vertexId <-c(1:vcount(source_graph))
	
	#ensure numerical lables
	#V(export_graph)$name <-vertexId
	V(export_graph)$NodeId <-vertexId
	export_graph <-attachVertexAtt(export_graph,source_graph,
		"Label","label")
	export_graph <-attachVertexAtt(export_graph,source_graph,
		"ColorName","frame.color")
	# convert rectangle to square
	V(source_graph)$frame.shape[V(source_graph)$frame.shape=="rectangle"]="square"
	export_graph <-attachVertexAtt(export_graph,source_graph,
		"NodeShape","frame.shape")
	export_graph <-attachVertexAtt(export_graph,source_graph,
		"NodeSize","vertex.size")
	export_graph <-attachVertexAtt(export_graph,source_graph,
		"StartTime","start.time")
	export_graph <-attachVertexAtt(export_graph,source_graph,
		"EndTime","end.time")

	export_graph <-attachEdgeAtt(export_graph,source_graph,
		"ArcWeight","weight")
	export_graph <-attachEdgeAtt(export_graph,source_graph,
		"ArcWidth","arrow.size")
	export_graph <-attachEdgeAtt(export_graph,source_graph,
		"ColorName","color")
	export_graph <-attachEdgeAtt(export_graph,source_graph,
		"StartTime","start.time")
	export_graph <-attachEdgeAtt(export_graph,source_graph,
		"EndTime","end.time")

	#build the vertex string
	
	sonFile<-list()
	sonVertex<-list()
	sonVertexLine<-list()
	
	#set vertex data
	listVertexAtts <-list.vertex.attributes(export_graph)
	
	for(i in 0:maxVertex)
	{
		for(ii in 1:length(listVertexAtts))
		{	
			if(ii>1)
			{
				sonVertexLine <-paste(sonVertexLine,get.vertex.attribute(export_graph,listVertexAtts[ii],index=i),sep="\t")
			}
			else
			{
				sonVertexLine <-paste(get.vertex.attribute(export_graph,listVertexAtts[ii],index=i),sep="\t")
			}
		}
		if(i>0)
		{
			sonVertex <- paste(c(sonVertex,sonVertexLine),sep="\r\n")
		}
		else
		{
			sonVertex <- paste(c(sonVertexLine),sep="\r\n")
		}	
	}
	sonVertexLine<-list()
	
	# set vertex header, must loop because of igraph oddities
	sonVertexHeader<-list()
	for(i in 1:length(listVertexAtts))
	{	
		if(i>1)
		{
			sonVertexHeader <-paste(sonVertexHeader,listVertexAtts[i],sep="\t")
		}
		else
		{
			sonVertexHeader <-paste(listVertexAtts[i],sep="\t")
		}
	}
	
	sonVertex <- paste(c(sonVertexHeader,sonVertex),sep="\r\n")
	
	#ensure numerical lables
	V(export_graph)$name <-vertexId
	sonEdge<-list()
	sonEdgeLine<-list()
	maxEdge <-ecount(source_graph)-1
	listEdgeAtts <-list.edge.attributes(export_graph)
	sonEdgeList <- get.edgelist(export_graph)
	#set egde data
	for(i in 0:maxEdge)
	{
		for(ii in 1:length(listEdgeAtts))
		{	
			if(ii>1)
			{
				sonEdgeLine <-paste(sonEdgeLine,get.edge.attribute(export_graph,listEdgeAtts[ii],index=i),sep="\t")
			}
			else
			{
				sonEdgeLine <-paste(get.edge.attribute(export_graph,listEdgeAtts[ii],index=i),sep="\t")
			}
		}
		sonEdgeLine <-paste(sonEdgeList[i+1,1],sonEdgeList[i+1,2],sonEdgeLine ,sep="\t")
		
		
		if(i>0)
		{
			sonEdge <- paste(c(sonEdge,sonEdgeLine),sep="\r\n")
		}
		else
		{
			sonEdge <- paste(c(sonEdgeLine),sep="\r\n")
		}	
	}


	# set edgelist header, must loop because of igraph oddities
	sonEdgeHeader<-list()
	
	for(i in 1:length(listEdgeAtts))
	{	
		if(i>1)
		{
			sonEdgeHeader <-paste(sonEdgeHeader,listEdgeAtts[i],sep="\t")
		}
		else
		{
			sonEdgeHeader <-paste("FromId","ToId",listEdgeAtts[i],sep="\t")
		}
	}

	sonEdge <- paste(c(sonEdgeHeader,sonEdge),sep="\r\n")

	
	sonfile<-paste(c(sonVertex,sonEdge),sep="\r\n")

	outfile <- file(fileN, "w", encoding="UTF-8")
	writeLines(sonfile, con = outfile, sep = "\n", useBytes = FALSE)
	close(outfile)
	print("file exported")
	
}

attachVertexAtt <- function(export_graph,source_graph,
	exportAttName,sourceAttName)
{

	attValue <-get.vertex.attribute(source_graph,sourceAttName)
	if(length(attValue>0))
	{
		export_graph<-set.vertex.attribute(export_graph,exportAttName, 
			value=c(get.vertex.attribute(source_graph,sourceAttName)))
	}
	
	return(export_graph)	
}

attachEdgeAtt <- function(export_graph,source_graph,
	exportAttName,sourceAttName)
{

	attValue <-get.edge.attribute(source_graph,sourceAttName)
	if(length(attValue>0))
	{
		export_graph<-set.edge.attribute(export_graph,exportAttName, 
			value=get.edge.attribute(source_graph,sourceAttName))
	}
	
	return(export_graph)	
}


