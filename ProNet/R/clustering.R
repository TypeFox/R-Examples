##' @title Network clustering
##' 
##' @description Network clustering based on different methods.
##'
##' @param graph An igraph object.
##' @param method Clustering method, possible values are \code{FN}, \code{MCL}, \code{LINKCOMM} and \code{MCODE}.
##' @param expansion Numeric value > 1 for the expansion parameter, if \code{method} is \code{MCL}. See \pkg{MCL} for more information.
##' @param inflation Numeric value > 0 for the inflation power coefficient, if \code{method} is \code{MCL}. See \pkg{MCL} for more information.
##' @param hcmethod A character string naming the hierarchical clustering method to use. Default value is \code{average}. See \pkg{linkcomm} for more information.
##' @param directed Logical value, indicating whether the network is directed (if \code{TRUE}) or not (if \code{FLASE}).
##' @param outfile File to save the clustering result.
##' @param plot Logical value, indicating whether to plot summary output (if \code{TRUE}) or not (if \code{FLASE}). 
##' @param layout Mode of the layout, possible values are \code{fruchterman.reingold}, \code{reingold.tilford}, \code{random}, \code{circle}, \code{kamada.kawai}, \code{lgl} and \code{sphere}. See \pkg{igraph} for more information.
##' @param ... Other arguments.
##' @return A summary and visulization of the clustering.
##' @seealso \code{\link{mcode}}
##' @references A Clauset, MEJ Newman, C (2004) Moore: Finding community structure in very large networks, \url{http://www.arxiv.org/abs/cond-mat/0408187}.
##' @references van Dongen, S.M. (2000) Graph Clustering by Flow Simulation. Ph.D. thesis, Universtiy of Utrecht.
##' @references Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. Bioinformatics 27 (14), 2011-2012.
##' @references Bader GD, Hogue CW. An automated method for finding molecular complexes in large protein interaction networks. BMC Bioinformatics. 2003 Jan 13;4(1):2.
##' @export
##' @examples
##' nlocal<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construction(input=nlocal,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)
##' cluster(net,method="MCODE",layout="fruchterman.reingold")
##' cluster(net,method="FN",layout="fruchterman.reingold")

cluster<-function(graph,method=c("FN","MCL","LINKCOMM","MCODE"),
                  expansion=2,inflation=2,hcmethod="average",directed= FALSE,outfile=NULL,plot=TRUE,
                  layout=c("reingold.tilford","circle","random","fruchterman.reingold","sphere","kamada.kawai","lgl"),...)
{
  method<-match.arg(method)
  if(method=="FN"){
	  graph<-simplify(graph)
	  fc=fastgreedy.community(graph,merges=TRUE,modularity=TRUE)
	  membership<-membership(fc)
	  if(!is.null(V(graph)$name)){
      names(membership)<-V(graph)$name
	  }
	  if(plot){
			color<-c("aquamarine","bisque","chartreuse2","chocolate1","cyan","deeppink","goldenrod1","plum1","magenta","slateblue1","seashell","red1")
			if(length(unique(membership))>12){
        color<-rainbow(length(unique(membership)))
			}
			color<-color[membership]
			visualization(graph=graph,layout=layout,node.fill.color=color,...)
	  }
	  if(!is.null(outfile)){
		  cluster.save(cbind(names(membership),membership),outfile=outfile)		
	  }else{
		  return(membership)
	  }
  }else if(method=="LINKCOMM"){
	  edgelist<-get.edgelist(graph)
	  if(!is.null(E(graph)$weight)){
      edgelist<-cbind(edgelist,E(graph)$weight)
	  }
	  lc<-getLinkCommunities(edgelist,plot=FALSE,directed=directed,hcmethod=hcmethod)
	  if(!is.null(outfile)){
		  cluster.save(lc$nodeclusters,outfile=outfile)
		  if(plot){
			  layout<-match.arg(layout)
			  layout<-paste("layout.",layout,sep="")
			  layout<-eval(parse(text=layout))
			  op<-par(mai=c(rep.int(0.1,4)))
			  plotLinkCommGraph(lc,layout=layout,...)
			  par(op)
		  }
	  }else{
		  if(plot){
			  layout<-match.arg(layout)
			  layout<-paste("layout.",layout,"(graph)",sep="")
			  layout<-eval(parse(text=layout))
			  plotLinkCommGraph(lc,layout=layout,...)
		  }
		  return(lc$nodeclusters)
	  }
  }else if(method=="MCL"){
    adj<-matrix(rep(0,length(V(graph))^2),nrow=length(V(graph)),ncol=length(V(graph)))
    for(i in 1:length(V(graph))){
      neighbors<-neighbors(graph,v=V(graph)$name[i],mode="all")
      j<-match(neighbors$name,V(graph)$name,nomatch=0)
      adj[i,j]=1
    }
    lc<-mcl(adj,addLoops=TRUE,expansion=expansion,inflation=inflation,allow1=TRUE,max.iter=100,ESM=FALSE)
    lc$name<-V(graph)$name
    lc$Cluster<-lc$Cluster
    if(plot){	
      color<-c("white","aquamarine","bisque","chartreuse2","chocolate1","cyan","deeppink","goldenrod1","plum1","magenta","slateblue1","seashell","red1")
      if(length(unique(lc$Cluster))>12){
        color<-rainbow(length(lc$Cluster))
      }
      color<-color[lc$Cluster]
      color[lc$Cluster==1]<-"white"
      visualization(graph=graph,node.fill.color=color,outfile=outfile,layout=layout,...)
    }
    if(!is.null(outfile)){
      cluster.save(cbind(lc$name,lc$Cluster),outfile=outfile)		
    }else{
      result<-lc$Cluster
      names(result)<-V(graph)$name
      return(result)
    }
  }else if(method=="MCODE"){
	  compx<-mcode(graph,vwp=0.9,haircut=T,fluff=T,fdt=0.1)
	  index<-which(!is.na(compx$score))
	  membership<-rep(0,vcount(graph))
	  for(i in 1:length(index)){
		  membership[compx$COMPLEX[[index[i]]]]<-i
	  }
	  if(!is.null(V(graph)$name)) names(membership)<-V(graph)$name
	  if(plot){	
		  color<-c("white","aquamarine","bisque","chartreuse2","chocolate1","cyan","deeppink","goldenrod1","plum1","magenta","slateblue1","seashell","red1")
		  if(length(unique(membership))>12){
        color<-rainbow(length(unique(membership)))
		  }
		  color<-color[membership+1]
		  color[membership==1]<-"white"
		  visualization(graph=graph,node.fill.color=color,layout=layout,...)
	  }
	  if(!is.null(outfile)){
		  cluster.save(cbind(names(membership),membership),outfile=outfile)
		  invisible(NULL)
	  }else{
		  return(membership)
	  }
  }
}

## Save the clustering result.
cluster.save<-function(membership,outfile){
	wd<-dirname(outfile)
	wd<-ifelse(wd==".",paste(wd,"/",sep=""),wd)
	filename<-basename(outfile)
	if((filename=="")||(grepl(":",filename))){
		filename<-"membership.txt"
	}else if(grepl("\\.",filename)){
		filename<-sub("\\.(?:.*)",".txt", filename)
	}
	write.table(membership,file=paste(wd,filename,sep="/"),
              row.names=FALSE,col.names=c("node","cluster"),quote =FALSE)
}
	
	
	
	
