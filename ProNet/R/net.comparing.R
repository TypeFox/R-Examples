##' @title Comparing two networks
##' 
##' @description Comparing two networks' topological parameters.
##' 
##' @param graph1 An igraph object.
##' @param graph2 An igraph object.
##' @param topology.parameters Logical value, indicating whether to do basic compariring (if \code{TRUE}) or not (if \code{FALSE}).
##' @param degree Logical value, indicating whether to do degree comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param cc Logical value, indicating whether to do clustering coefficient comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param betweenness Logical value, indicating whether to do betweenness comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param eccentricity Logical value, indicating whether to do eccentricity comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param ave.path.len Logical value, indicating whether to do average path comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param figure.type Type of the plot. See \code{\link{plot}} for more information.
##' @param method Test method, only \code{utest} is supported.
##' @param legendname Legend name for graph.
##' @return A list of compared results and plot.
##' @references Dehmer,M. et al. A large scale analysis of information-theoretic network complexity measures using chemical structures. PLoS ONE, (2009a), 4, e8057.
##' @seealso \code{\link{comp.subnet}} and \code{\link{comp.rand.subnet}}.
##' @export
##' @examples
##' g1<-barabasi.game(100,power=0.8,directed = FALSE)
##' g2<-erdos.renyi.game(100,p=0.01,directed = FALSE)
##' res<-net.comparing(g1,g2,degree=TRUE,cc=TRUE,betweenness=FALSE)

net.comparing<-function(graph1,graph2,topology.parameters=FALSE,
                        degree=FALSE,cc=FALSE,betweenness=FALSE,
                        eccentricity=FALSE,ave.path.len=FALSE,
                        figure.type=1,method="utest",legendname)
{
	if((!is.igraph(graph1))||(!is.igraph(graph2))){
		stop("No igraph object!")
	}
	method<-match.arg(method)
	if(figure.type==1 | figure.type==3){
    type="p"
	}else if(figure.type==2 | figure.type==4){
    type="l"
	}
	
	legendname<-c(substitute(graph1),substitute(graph2))
	
	result<-list()
	if(topology.parameters){
		result$topology<-cbind(topology_parameters(graph1),
                           topology_parameters(graph2))
		colnames(result$topology)<-legendname
		result$topology<-round(result$topology,4)
	}
	if(degree){
		dg1<-degree(graph1)
		dg2<-degree(graph2)
		result$dg.p.value<-wilcox.test(dg1,dg2)$p.value
		comparing_plot(dg1,dg2,xlab="Degree",ylab="Density",type=type,
                   legendname=legendname,figure.type=figure.type)
	}
	if(cc){
		cc1<-transitivity(graph1,type="local",vids=V(graph1))
		cc2<-transitivity(graph2,type="local",vids=V(graph2))
		result$cc.p.value<-wilcox.test(cc1,cc2)$p.value
		comparing_plot(cc1,cc2,xlab="Cluster Coefficient",ylab="Density",type=type,
                   legendname=legendname,figure.type=figure.type)
	}
	if(betweenness){
		bw1<-betweenness(graph1)
		bw2<-betweenness(graph2)
		result$btw.p.value<-wilcox.test(bw1,bw2)$p.value
		comparing_plot(x1=bw1,x2=bw2,xlab="Betweenness",ylab="Density",type=type,
                   legendname=legendname,figure.type=figure.type)
	}
	if(ave.path.len){
		if(is.directed(graph1) & is.directed(graph2)){
			avep1<-as.vector(shortest.paths(graph1))
			avep1<-avep1[avep1!=0]
			avep2<-as.vector(shortest.paths(graph2))
			avep2<-avep2[avep2!=0]
		}else{
			avep1<-shortest.paths(graph1)[lower.tri(shortest.paths(graph1))]
			avep2<-shortest.paths(graph2)[lower.tri(shortest.paths(graph2))]
		}
		result$ave.path.len.p.value<-wilcox.test(avep1,avep2)$p.value
		comparing_plot(x1=avep1,x2=avep2,xlab="Ave.path.len",ylab="Density",
                   type=type,legendname=legendname,figure.type=figure.type)
	}
	if(eccentricity){
		ec1<-eccentricity(graph1)
		ec2<-eccentricity(graph2)
		result$ec.p.value<-wilcox.test(ec1,ec2)$p.value
		comparing_plot(x1=ec1,x2=ec2,xlab="Eccentricity",ylab="Density",
                   type=type,legendname=legendname,figure.type=figure.type)
	}
  
	return(result)
}

## Plot the comparison of the graphs' topological parameters.
comparing_plot<-function(x1,x2,pch=c(17,19),col=c(2,3),
                         type="p",legendname=c("graph1","graph2"),figure.type,...)
{
	x1<-as.data.frame(table(x1)/length(x1),stringsAsFactors=FALSE,responseName="Density")
	x2<-as.data.frame(table(x2)/length(x2),stringsAsFactors=FALSE,responseName="Density")
	x1[,1]<-as.numeric(x1[,1])
	x2[,1]<-as.numeric(x2[,1])
	if(figure.type<3){
		ylim<-as.numeric(c(min(min(x1[,2]),min(x2[,2])),max(max(x1[,2]),max(x2[,2]))))
		xlim<-as.numeric(c(min(min(x1[,1]),min(x2[,1])),max(max(x1[,1]),max(x2[,1]))))
		xlim[1]<-xlim[1]-0.01
		plot(x1,type=type,pch=pch[1],col=col[1],
         xlim=xlim,ylim=ylim,bty="l",...)
		par(new=TRUE)
		plot(x2,type=type,pch=pch[2],col=col[2],xlab="",ylab="",main="",
         axes=FALSE,xlim=xlim,ylim=ylim)
		legend("topright",legend=legendname,pch=pch,col=col)
	}else{
		op<-par(mfrow=c(1,2),mar=c(5,4,4,1))
		plot(x1,type=type,pch=pch[1],col=col[1],main=legendname[1],...)
		plot(x2,type=type,pch=pch[2],col=col[2],main=legendname[2],...)
		par(op)
	}
}

##  calculate basic topological parameters.
topology_parameters<-function(graph)
{
	if(!is.igraph(graph)){
    stop("not a igraph object")
	}
	result<-c(vcount(graph),ecount(graph),sum(degree(graph)==0),clusters(graph)$no,
            diameter(graph),average.path.length(graph),
            mean(neighborhood.size(graph,1))-1,mean(degree(graph)),
            transitivity(graph,type="average"),mean(betweenness(graph)))
	names(result)<-c("Number of nodes","Number of edges","Isolated nodes",
                   "Connected components","Network diameter","Average path length",
                   "Avg. number of neighbors","Ave. degree",
                   "Avg. clustering coefficient","Avg. betweenness")
	return(result)
}


