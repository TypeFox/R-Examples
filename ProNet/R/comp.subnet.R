##' @title Comparing a sub network to the whole one
##' 
##' @description Comparing the topological parameters of a sub network and the whole one.
##' 
##' @param x Vertex of the sub network.
##' @param graph An igraph object.
##' @param topology.parameters Logical value, indicating whether to do basic compariring (if \code{TRUE}) or not (if \code{FALSE}).
##' @param degree Logical value, indicating whether to do degree comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param cc Logical value, indicating whether to do clustering coefficient comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param betweenness Logical value, indicating whether to do betweenness comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param eccentricity Logical value, indicating whether to do eccentricity comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param ave.path.len Logical value, indicating whether to do average path comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param figure.type Type of the plot. See \code{\link{plot}} for more information.
##' @param method Test method, currently only \code{utest} is supported.
##' @param legendname Legend name for the plot.
##' @return A list of compared parameters and plot.
##' @seealso \code{\link{comp.rand.subnet}}, \code{\link{net.comparing}}
##' @export
##' @examples
##' g<-barabasi.game(100,power=0.8,directed = FALSE)
##' id<-sample(1:100, 20)
##' res<-comp.subnet(id,g)
##' res<-comp.subnet(id,g,topology.parameters=TRUE)

comp.subnet<-function(x,graph,topology.parameters=TRUE,degree=FALSE,
                      cc=FALSE,betweenness=FALSE,eccentricity=FALSE,
                      ave.path.len=FALSE,figure.type=1,method="utest",
                      legendname=c(substitute(x),substitute(graph)))
{
	len<-length(x)
	if(len==0){
    stop("x should have more than one element")
	}
	method<-match.arg(method)
	if(is.character(x)){
		id<-match(x,V(graph)$name,nomatch=0)
		id<-id[id!=0]
		len<-length(id)
	}else if(is.numeric(x)){	
		id<-x[which(x>0 & x<=vcount(graph))]
	}else if(is.igraph(x)){
		id<-match(V(x)$name, V(graph)$name,nomatch=0)
		id<-id[id!=0]
		len<-length(id)
	}
	
	if(figure.type==1 | figure.type==3){
    type="p"
	}else if(figure.type==2 | figure.type==4){
    type="l"
	}
  
	result<-list()
	
  if(topology.parameters){
		result$topology<-topology_parameters2(graph,id)
		colnames(result$topology)<-legendname
	}
  
	if(degree){
		dg<-degree(graph)
		result$dg.p.value<-wilcox.test(dg,dg[id])$p.value
		comparing_plot(dg,dg[id],xlab="Degree",ylab="Density",type=type,
                   legendname=legendname,figure.type=figure.type)
	}
	if(cc){
		cc<-transitivity(graph,type="local",vids=V(graph))
		result$cc.p.value<-wilcox.test(cc,cc[id])$p.value
		comparing_plot(cc,cc[id],xlab="Cluster Coefficient",ylab="Density",type=type,
                   legendname=legendname,figure.type=figure.type)
	}
	if(betweenness){
		bw<-betweenness(graph)
		result$btw.p.value<-wilcox.test(bw,bw[id])$p.value
		comparing_plot(bw,bw[id],xlab="Betweenness",ylab="Density",type=type,
                   legendname=legendname,figure.type=figure.type)
	}
	if(ave.path.len){
	  if(is.directed(graph)){
	    avep1<-as.vector(shortest.paths(induced.subgraph(graph,id)))
		  avep1<-avep1[avep1!=0]
		  avep2<-as.vector(shortest.paths(graph))
		  avep2<-avep2[avep2!=0]
  	}else{
			avep1<-shortest.paths(graph)[lower.tri(induced.subgraph(graph,id))]
			avep2<-shortest.paths(graph)[lower.tri(shortest.paths(graph))]
		}
		result$ave.path.len.p.value<-wilcox.test(avep1,avep2)$p.value
		comparing_plot(x1=avep1,x2=avep2,xlab="Ave.path.len",ylab="Density",type=type
		 ,legendname=legendname,figure.type=figure.type)
	}
	if(eccentricity){
		ec<-eccentricity(graph)
		result$ec.p.value<-wilcox.test(ec,ec[id])$p.value
		comparing_plot(ec,ec[id],xlab="Eccentricity",ylab="Density",type=type,
                   legendname=legendname,figure.type=figure.type)
	}
  
	return(result)
}

## Plot the comparison of two graphs' topological parameters.
comparing_plot<-function(x1,x2,pch=c(17,19),col=c(2,3),type="p",
                         legendname=c(substitute(graph1),substitute(graph)),
                         figure.type,...)
{
	x1<-as.data.frame(table(x1)/length(x1),stringsAsFactors=FALSE,responseName="Density")
	x2<-as.data.frame(table(x2)/length(x2),stringsAsFactors=FALSE,responseName="Density")
	x1[,1]<-as.numeric(x1[,1])
	x2[,1]<-as.numeric(x2[,1])
	if(figure.type<3){
		ylim<-as.numeric(c(min(min(x1[,2]),min(x2[,2])),max(max(x1[,2]),max(x2[,2]))))
		xlim<-as.numeric(c(min(min(x1[,1]),min(x2[,1])),max(max(x1[,1]),max(x2[,1]))))
		xlim[1]<-xlim[1]-0.01
		plot(x1,type=type,pch=pch[1],col=col[1],xlim=xlim,ylim=ylim,bty="l",...)
		par(new=TRUE)
		plot(x2,type=type,pch=pch[2],col=col[2],xlab="",ylab="",main="",
         axes=FALSE,xlim=xlim,ylim=ylim)
		legend("topright",legend=legendname,pch=pch,col=col)
	}else{
		op<-par(mfrow=c(1,2),mar=c(5,4,4,1))
		plot(x1,type=type,pch=pch[1],col=col[1],...)
		plot(x2,type=type,pch=pch[2],col=col[2],...)
		par(op)
	}
}

## Calculate simple topological parameters of the sub- and total graphs.
topology_parameters2<-function(graph,id){
	if(!is.igraph(graph)){
    stop("not a igraph object")
	}
	bt<-betweenness(graph)
	res1<-c(vcount(graph),ecount(graph),sum(degree(graph)==0),clusters(graph)$no,
          mean(neighborhood.size(graph,1))-1,mean(degree(graph)),
          transitivity(graph,type="average"),mean(bt))
	subg<-induced.subgraph(graph,id)
	res2<-c(length(id),ecount(subg),sum(degree(subg)==0),clusters(subg)$no,
          mean(neighborhood.size(graph,1)[id])-1,mean(degree(graph)[id]),
          mean(transitivity(graph,type="local",vids=V(graph))[id]),mean(bt[id]))
	names(res1)<-names(res2)<-c("Number of nodes","Number of edges","Isolated nodes",
                              "Connected components","Avg. number of neighbors",
                              "Ave. degree","Avg. clustering coefficient","Avg. betweenness")
	return(cbind(res2,res1))
}

