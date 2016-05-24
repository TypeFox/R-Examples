##' @title Hierarchical plot of biological network
##' 
##' @description Hierarchical plot of biological network according to the elements' subcellular location.
##' 
##' @param graph An igraph object.
##' @param species The species name, currently only \code{human} and \code{ath} are available.
##' @param vertex.color Color of the vertex. Default value is \code{1}.
##' @param vertex.size Size of the vertex. Default value is \code{10}.
##' @param vertex.shape Shape of the vertex. Default value is \code{square}.
##' @param edge.color Color of the edge. Default value is \code{1}.
##' @param label.add Boolean value, whether to add label to the plot (if \code{TRUE}) or not (if \code{FALSE}).
##' @param colorbar.add Boolean value, whether to add colorbar to the plot (if \code{TRUE}) or not (if \code{FALSE}).
##' @param xlim A vector indicating the range of x axis.
##' @param ylim A vector indicating the range of y axis.
##' @param ... Other arguments.
##' @return A hierarchical plot of biological network.
##' @references Barsky A, Gardy JL, Hancock REW, and Munzner T. (2007) Cerebral: a Cytoscape plugin for layout of and interaction with biological networks using subcellular localization annotation. Bioinformatics 23(8):1040-2.
##' @export
##' @examples
##' gene<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construction(input=gene,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)
##' location(net,species="human",vertex.color="vertex.hierarchy")

location<-function(graph,species=c("human","ath"),vertex.color=1,
                   vertex.size=10,vertex.shape="square",edge.color=1,
                   label.add=TRUE,colorbar.add=TRUE,xlim=c(-1,1),ylim=c(-1,1),...)
{
  locDB<-NULL
  data(locDB, envir = environment())
  
	species<-match.arg(species)
    
  loc<-locDB[[species]][,c(3,1)]
	id<-match(V(graph)$name,loc[,1],nomatch=0)
	graph<-set.vertex.attribute(graph,"category",index=as.character(loc[id,1]),
                            value=as.character(loc[id,2]))
  
	idna<-which(!is.na(V(graph)$category))
	if(length(idna)>0){
		graph<-induced.subgraph(graph,which(!is.na(V(graph)$category)))
		warning("There are some 'NA' values in graph attribute 'category'!")
	}
	##  set vertex color
	colorpalette<-vertex.color
	if(is.character(vertex.color)){
		attrib<-get.vertex.attribute(graph,vertex.color)
		if(!is.null(attrib)){
			colr<-attrib2color(attrib)
			colorpalette<-colr$colorpalette
			vertex.color<-colr$color
		}else{
      warning("Make sure right vertex attribute for vertex color!")
		}
	}
	op<-par(mai=c(0,0,0,0))
	treeplot(graph,species=species,vertex.color=vertex.color,vertex.size=vertex.size,
           vertex.shape=vertex.shape,colorbar.add=colorbar.add,colorpalette=colorpalette,
           xlim=xlim,ylim=ylim,label.add=label.add,...)
	par(op)
}

## Hierarchical visulization.
treeplot<-function(g,species=c("human","ath"),vertex.size=10,vertex.shape="circle",
                   vertex.color=3,colorbar.add=TRUE,colorpalette,xlim=c(-1,-1),ylim=c(-1,1),label.add,...)
{
	species<-match.arg(species)
	num_attrib<-num_category(V(g)$category,species)
	num_line<-unlist(lapply(num_attrib,function(x){sum(x,na.rm=T)}))
	num_line[num_line==0]<-2
	line_range<-(num_line)^(0.6)
	line_range<-line_range/sum(line_range)
	min_range<-min(line_range)
	##  line coordinate
	line_range<-ylim[2]-cumsum(line_range)*(ylim[2]-ylim[1]);
	##  check the min_range's job
	line_range<-cbind(line_range,c(ylim[2],line_range[-length(num_line)])-min_range*0.1)
	##  set frame 
	frames<-Frames2(line_range,num_attrib,vertex.size=vertex.size,xlim=xlim,ylim=ylim)
	L<-layout.Frame(V(g)$category,frames,g,vertex.size=vertex.size,num_attrib=num_attrib)
	vertex.label<-get.vertex.attribute(g,"name")
	if(is.null(vertex.label))vertex.label<-1:vcount(g)
	plot(g,layout=L,vertex.size=vertex.size,vertex.shape=vertex.shape,
       vertex.color=vertex.color,vertex.label=vertex.label,rescale=FALSE,
       xlim=c(xlim[1],xlim[2]+0.04*diff(xlim)),
       ylim=c(ylim[1]-0.04*diff(ylim),ylim[2]),...)
	##	plot frames
	plot.frames(frames,species,label.add,xlim,ylim,L,g)
	rm(L)
	##    hline
	hline_coordinate<-line_range[-nrow(line_range),1]
	for(i in hline_coordinate) lines(xlim,c(i,i),lwd=3)
	##  plot  color bar 
	if(colorbar.add){
		ucord<-par()$usr
		pin<-par()$pin
    dy<-0.02*(ucord[4]-ucord[3])
		dx<-0.2*pin[2]*(ucord[2]-ucord[1])/(pin[1])
		xs<-seq(0,dx,,length(colorpalette)+1)+xlim[1]
		ys<-seq(0,dy,,1+1)+ylim[1]-dy
		image(xs,ys,matrix(1:length(colorpalette),length(colorpalette),1),
          add=TRUE,col=colorpalette)
	}
}

##  Calculate the number of vertices belong to each category
num_category<-function(attrib,species)
{
	if(length(attrib)==0){
	  stop("Not a vector")
	}
	result<-list()
	attrib_table<-table(attrib)
	result<-lapply(seq_along(.category[[species]]),function(i){
		res<-attrib_table[.category[[species]][[i]]]
		names(res)<-.category[[species]][[i]]
		res})
	return(result)
}

## Calculate Frame 
Frame<-function(line_range,num_attrib,vertex.size,xlim,ylim)
{
	if(!is.matrix(line_range)){
		stop("Not a matrix")
	}
	nr<-nrow(as.data.frame(line_range))
	if(nr!=length(num_attrib)){
		stop("category is wrong ")
	}
	line_range[,1]<-line_range[,1]+0.020*diff(ylim)
	line_range[,2]<-line_range[,2]-0.010*diff(ylim)
	xgap <- 0.010*diff(xlim)
	result<-lapply(1:nr,function(i){
		if(i!=2)
			ylimit<-matrix(line_range[i,],length(num_attrib[[i]]),2,byrow=T)
		else 
			ylimit<-matrix(c(line_range[i,1]-0.01*diff(ylim),line_range[i,2]),
                     length(num_attrib[[i]]),2,byrow=T)
		n<-length(num_attrib[[i]])
		if(n>1){
			xlen<-diff(xlim)-(n-1)*xgap
			id<-which(is.na(num_attrib[[i]]))
			sum_attrib<-sum(num_attrib[[i]],na.rm=TRUE)+length(id)
			num_attrib[[i]][num_attrib[[i]]<=2]<-max(2,max(vertex.size)/200*sum_attrib/xlen)
			num_attrib[[i]][id]<-2
			sum_attrib<-sum(num_attrib[[i]])
			x<-c(xlim[1],xlim[1]+num_attrib[[i]][1]/sum_attrib*xlen)
			for(j in 2:n){
				x<-c(x,x[length(x)]+xgap)
				x<-c(x,x[length(x)]+num_attrib[[i]][j]/sum_attrib*xlen)
			}
			xmin<-x[2*(1:n)-1]
			xmax<-x[2*(1:n)]
		}else{
			xmin<-xlim[1]
			xmax<-xlim[2]
		}
		result<-cbind(xmin,xmax,ylimit)
		rownames(result)<-names(num_attrib[[i]])
		return(result)
	})
	return(result)
}

## Calculate Frames2
Frames2<-function(line_range,num_attrib,vertex.size,xlim,ylim)
{
	if(!is.matrix(line_range)){
		stop("Not a matrix")
	}
	nr<-nrow(as.data.frame(line_range))
	if(nr!=length(num_attrib)){
		stop("category is wrong ")
	}
	line_range[,1]<-line_range[,1]+0.020*diff(ylim)
	line_range[,2]<-line_range[,2]-0.010*diff(ylim)
	xgap<-0.010*diff(xlim)
	result<-lapply(1:nr,function(i){
		if(i!=2)
			ylimit<-matrix(line_range[i,],length(num_attrib[[i]]),2,byrow=T)
		else 
			ylimit<-matrix(c(line_range[i,1]-0.01*diff(ylim),line_range[i,2]),
                     length(num_attrib[[i]]),2,byrow=T)
		n<-length(num_attrib[[i]])
		if(n>1){
			xlen<-diff(xlim)-(n-1)*xgap
			id<-which(is.na(num_attrib[[i]]))
			sum_attrib<-sum(num_attrib[[i]],na.rm=TRUE)+length(id)
			num_attrib[[i]][num_attrib[[i]]<=2]<-max(2,max(vertex.size)/200*sum_attrib/xlen)
			num_attrib[[i]][id]<-2
			sum_attrib<-sum(num_attrib[[i]])
			if(n==4){
				xlen<-diff(xlim)-xgap
				ygap<-min(0.010*diff(ylim),xgap)
				ylen<-max(ylimit)-min(ylimit)-ygap
				s <- num_attrib[[i]]/sum_attrib
				ylimit[1:2,1]<-ylimit[1,2]-sum(s[1:2])/sum(s)*ylen
				ylimit[3:4,2]<-ylimit[1:2,1]-ygap
				x1<-xlim[1]+c(0,s[1]/sum(s[1:2])*xlen)
				x2<-c(x1[2]+xgap,xlim[2])
				x3<-xlim[1]+c(0,s[2]/sum(s[3:4])*xlen)
				x4<-c(x3[2]+xgap,xlim[2])
				xmin<-c(x1[1],x2[1],x3[1],x4[1])
				xmax<-c(x1[2],x2[2],x3[2],x4[2])
			}else{
				x<-c(xlim[1],xlim[1]+num_attrib[[i]][1]/sum_attrib*xlen)
				for(j in 2:n){
					x<-c(x,x[length(x)]+xgap)
					x<-c(x,x[length(x)]+num_attrib[[i]][j]/sum_attrib*xlen)
				}
				xmin<-x[2*(1:n)-1]
				xmax<-x[2*(1:n)]
			}
		}else{
			xmin<-xlim[1]
			xmax<-xlim[2]
		}

		result<-cbind(xmin,xmax,ylimit)
		rownames(result)<-names(num_attrib[[i]])
		return(result)
	})
  
	return(result)
}
## Calculate layout based on Frame
layout.Frame<-function(attrib,frames,g,vertex.size,num_attrib){
	if(!is.character(attrib)){
    stop("Not a character")
	}
	if(!is.list(frames)){
    stop("Not a list")
	}
	frames<-do.call(rbind,frames)
	id<-which(is.na(unlist(num_attrib)))
	if(length(id)){
    frames<-frames[-id,]
	}
	catgry_name<-rownames(frames)
	mvs<-c(max(vertex.size)/200+1/100)
	if(length(catgry_name)){
		frames[,3]<-frames[,3]+mvs
		frames[,4]<-frames[,4]-mvs
		frames[,1]<-frames[,1]+mvs
		frames[,2]<-frames[,2]-mvs
		L<-matrix(0,length(attrib),2)
		for(i in seq_along(catgry_name)){
			idx<-which(attrib==catgry_name[i])
			Ltemp<-layout.fruchterman.reingold(induced.subgraph(g,idx))
			Lframe<-c(min(Ltemp[,1]),max(Ltemp[,1]),min(Ltemp[,2]),max(Ltemp[,2]))
			Ltemp[,1]<-mean(frames[i,1:2])+(Ltemp[,1]-mean(Lframe[1:2]))/(Lframe[2]-Lframe[1]+0.00001*mvs)*(frames[i,2]-frames[i,1])
			Ltemp[,2]<-mean(frames[i,3:4])+(Ltemp[,2]-mean(Lframe[3:4]))/(Lframe[4]-Lframe[3]+0.00001*mvs)*(frames[i,4]-frames[i,3])
			L[idx,]<- Ltemp
		}
		return(L)
	}
}

##	Plot frames
plot.frames<-function(frames,species,label.add,xlim,ylim,L,g){
	if(!is.list(frames)){
    stop("Not a list")
	}
	frames<-data.frame(do.call(rbind,frames))
	if(species=="human"){
	  id<-match(rownames(frames),c(.category[["human"]][[2]],.category[["human"]][[5]][2:3]))
	}else if(species=="ath"){ 
	  id<-match(rownames(frames),c(.category[["ath"]][[2]],.category[["ath"]][[5]][2:4]))
	}
	for(i in which(!is.na(id))){
		rect(frames[i,1],frames[i,3],frames[i,2],frames[i,4],border="deepskyblue",
         col=rgb(0,1,1,0.1),lty=1,lwd=2)
	}
	frames["cytoplasm",3]<-frames["chloroplast",3]
	if(label.add){
		d1<-0.08*diff(xlim)
		d2<-0.020*diff(ylim)
		updown<-1
		for(i in which(!is.na(id))){
			if((frames[i,2]-frames[i,1]<d1) & updown==0){
					idx<-match(V(g)$category,rownames(frames)[(i-1):i])
					temp<-c(frames[i,3],sort(L[!is.na(idx),2]),frames[i,4])
					idx<-which.max(diff(temp))
					text(0.618*frames[i,1]+0.382*frames[i,2],mean(temp[idx:(idx+1)]),rownames(frames)[i],cex=0.7)
					updown<-1
					rm(temp,idx)
			}else{
				text(0.618*frames[i,1]+0.382*frames[i,2],frames[i,4]+d2*0.1,rownames(frames)[i],cex=0.7)
				if(frames[i,2]-frames[i,1]<d1){
					updown<-0
				}
			}
		}
		for(i in which(is.na(id))){
			text(xlim[2]-d1*0.25,frames[i,3]-d2*0.5,rownames(frames)[i],cex=0.7)
		}
	}
}

## Set color based on attribute Value
attrib2color<-function(expV,Colors=c("red", "green", "yellow"))
{
	color<-rep("white",length(expV))
	index<-which(!is.na(expV))
	expression.value<-as.character(expV[index])
	s.e.v<-unique(expression.value)
	seriation<-order(s.e.v)
	names(seriation)<-s.e.v
	colorP<-colorRampPalette(Colors)(length(seriation))
	color[index]<-colorP[seriation[expression.value]]
	return(list(color=color,colorpalette=colorP))
}

##	List of location categories in "human" and "ath".
.category<-list()
.category[["human"]]<-list("plasma membrane",c("lysosome","secretory vesicles","endosome"),
                           "Golgi apparatus","endoplasmic reticulum",c("cytoplasm","mitochondrion","peroxisome"),"nucleus")
.category[["ath"]]<-list("plasma membrane",c("lysosome","endosome"),"Golgi apparatus","endoplasmic reticulum",
                         c("cytoplasm","mitochondrion","peroxisome","chloroplast"),"nucleus")
