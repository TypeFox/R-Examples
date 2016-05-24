##' @title Network visualization
##' 
##' @description This function takes an igraph formatted network as input and exports a picture of the visulization. Mode of the layout, and size, color, shape of the nodes, edges and labels can be set according to the experiment data.
##' 
##' @param graph An igraph object.
##' @param layout Mode of the layout, possible values are \code{fruchterman.reingold}, \code{reingold.tilford}, \code{random}, \code{circle}, \code{kamada.kawai}, \code{lgl} and \code{sphere}. See \pkg{igraph} for more information.
##' @param node.fill.color The fill color of the vertex. Default value is \code{SkyBlue2}.
##' @param node.border.color Border color the the vertex. Default value is \code{Black}.
##' @param node.shape Shape of the vertex. currently \code{circle}, \code{square}, \code{csquare}, \code{rectangle}, \code{crectangl}, \code{vrectangle} and \code{none} are supported.
##' @param node.size Size of the vertex. Default value is \code{15}.
##' @param node.label Label of the vertex. Specify \code{NA} to omit vertex labels. Default value is the vertex name.
##' @param node.label.color Color of labels of the vertex. Default value is \code{1}.
##' @param node.label.size Font size of the vertex label. Default value is \code{0.8}. 
##' @param node.label.position Positon of the vertex label. Default value is \code{0}.
##' @param edge.shape Shape of the edge. Default value is \code{1}.
##' @param edge.width Width of the edge. Default value is \code{1}.
##' @param edge.color Color of the edge. Default value is \code{gray1}.
##' @param ... Other arguments.
##' @return A plot of the network.
##' @export
##' @examples
##' local<-data.frame(1:5,2:6)
##' attribute<-data.frame(node=1:6,value=c(2.2,5.3,1.2,4.5,6.2,0.6))
##' net<-construction(input=local,local.net=TRUE,node.attribute=attribute)
##' visualization(net,layout="reingold.tilford")
##' visualization(net,layout="reingold.tilford",node.size=V(net)$value)

visualization<-function(graph,
                        layout=c("reingold.tilford","circle","random","fruchterman.reingold","sphere","kamada.kawai","lgl"),
                        node.fill.color="SkyBlue2",node.border.color="Black",
                        node.shape=c("circle","square","sphere","csquare","rectangle","crectangle","vrectangle","none","pie","raster"),
                        node.size=15,node.label="name",node.label.color=1,
                        node.label.size=0.8,node.label.position=0,edge.shape=1,
                        edge.width=1,edge.color="gray1",...)
{
	if(!is.igraph(graph)){
    stop("Not an igraph object")
	}
	if(missing(layout)){
    layout<-layout.fruchterman.reingold(graph)
	}else if(is.function(layout)){
    layout<-layout(graph)
	}else if(layout %in% c("reingold.tilford","circle","random","fruchterman.reingold","sphere","kamada.kawai","lgl")){
		layout<-eval(parse(text=paste("layout.",layout,"(graph)",sep="")))
	}
	layout<-layout.norm(layout,-1,1,-1,1)

	############################################
	## node.fill.color
	if(is.character(node.fill.color) & length(node.fill.color)==1){
		Expression<-get.vertex.attribute(graph,node.fill.color)
		if(!is.null(Expression)){
			node.fill.color<-Expression2color(Expression)
		}
	}
	############################################
	## node.size
	if(is.character(node.size) & length(node.size)==1){
		attrib<-get.vertex.attribute(graph,node.size)
		if(!is.null(attrib)){
			node.size<-attrib2Vs(attrib,method="unif",Vmax=15,Vmin=3)
		}
	}
	############################################
	## edge.width
	if(is.character(edge.width) & length(edge.width)==1){
		attrib<-get.edge.attribute(graph,edge.width)
		if(!is.null(attrib)){
			edge.width<-attrib2E(attrib,method="unif",Emax=6,Emin=1)
		}
	}	
	############################################
	## edge.color
	if(is.character(edge.color) & length(edge.color)==1){
		attrib<-get.edge.attribute(graph,edge.color)
		if(!is.null(attrib)){
			edge.color<-Expression2color(attrib)
		}
	}
	############################################
	## node.label
	if(is.character(node.label) & length(node.label)==1){
		attrib<-get.vertex.attribute(graph,node.label)
		if(!is.null(attrib)){
			node.label<-attrib
		}
	}
	############################################
	## node.label.position
	if(is.character(node.label.position) & length(node.label.position)==1){
		attrib<-get.vertex.attribute(node.label.position)
		if(!is.null(attrib)){
			node.label<-attrib2Vs(attrib,method="unif",Vmax=1,Vmin=-1)
		}
	}
	############################################
	## edge.shape	
	if(is.character(edge.shape) & length(edge.shape)==1){
		attrib<-get.vertex.attribute(edge.shape)
		if(!is.null(attrib)){
			edge.shape<-attrib2Vs(attrib,method="unif",Vmax=6,Vmin=0)
		}
	}
	###################################################
	##	node shape
	node.shape<-match.arg(node.shape)
	
	if(any(is.na(node.shape))){
		if(is.numeric(node.shape)){
			plot(graph,layout=layout,vertex.size=0,edge.lty=edge.shape,
           edge.width=edge.width,vertex.label=NULL,rescale=FALSE,...)

      par(new=T)
			points(layout,type="p",cex=node.size,col=node.border.color,pch=node.shape,bg=node.fill.color)

			par(new=T)
			plot(graph,layout=layout,vertex.size=0,vertex.label.color=node.label.color,
           vertex.label.cex=node.label.size,vertex.label.dist=node.label.position,
           edge.lty=0,vertex.label=node.label,rescale=FALSE,...)
      
			return(invisible(NULL))
		}
	}
  
	plot(graph,layout=layout,vertex.color=node.fill.color,vertex.frame.color=node.border.color,
       vertex.shape=node.shape,vertex.size=node.size,vertex.label.color=node.label.color,
       vertex.label.cex=node.label.size,vertex.label.dist=node.label.position,
       edge.lty=edge.shape,edge.width=edge.width,vertex.label=node.label,
       edge.color=edge.color,...)
  
	invisible(NULL)
}

## Set Color Based on Expression Value.
Expression2color<-function(expV,Colors=c("red", "green", "yellow"))
{
	color<-rep("white",length(expV))
	index<-which(!is.na(expV))
	expression.value<-as.character(expV[index])
	s.e.v<-unique(expression.value)
	seriation<-order(s.e.v)
	names(seriation)<-s.e.v
	colorP<-colorRampPalette(Colors)(length(seriation))
	color[index]<-colorP[seriation[expression.value]]
	return(color)
}

## Set vertex size Based on vertex attribute.
attrib2Vs<-function(attrib,method=c("unif","power"),Vmax=10,Vmin=3,Order=0.5)
{
	if(Vmax<=Vmin){
    stop("max(Node.size) less than min(Node.size)")
	}
	method=match.arg(method)
	VertexSize=rep(Vmin,length(attrib))
	indx<-which(!is.na(attrib))
	attrib<-attrib[indx]
	if(method=="unif"){
		VertexSize[indx]<-(attrib-min(attrib))/(max(attrib)-min(attrib))*(Vmax-Vmin)+Vmin
	}else if(method=="power"){
		VertexSize[indx]<-attrib^Order
	}
  
	return(VertexSize)
}

## Set edge size or edge width Based on edge attribute.
attrib2E<-function(attrib,method=c("unif","power"),Emax=1,Emin=0.1,Order=0.5)
{
	if(Emax<=Emin){
    stop("max(edge.width) less than min(edge.width)")
  }
	method=match.arg(method)
	EdgeSize<-rep(Emin,length(attrib))
	indx<-which(!is.na(attrib))
	attrib<-attrib[indx]
	if(method=="unif"){
		EdgeSize[indx]<-(attrib-min(attrib))/(max(attrib)-min(attrib))*(Emax-Emin)+Emin
	}else if(method=="power"){
		EdgeSize[indx]<-attrib^Order
	}
  
	return(EdgeSize)
}
