##' @title Comparing a sub network to the randomly simulated ones
##' 
##' @description Comparing a sub network with randomly simulated ones from the whole network.
##'  
##' @param subgraph An igraph object.
##' @param graph An igraph object. The whole one for random simulation.
##' @param nsim Times for simulation. Default value is \code{1000}.
##' @param degree Logical value, indicating whether to do vertex degree comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param betweenness Logical value, indicating whether to do betweenness comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param ave.path.len Logical value, indicating whether to do average path comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param eccentricity Logical value, indicating whether to do eccentricity comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param cc Logical value, indicating whether to do clustering coefficient comparing (if \code{TRUE}) or not (if \code{FALSE}).
##' @param method Test method, currently only \code{utest} is supported.
##' @param FDR False discovery rate. Default value is \code{0.05}.
##' @return A matrix of compared parameters and plots.
##' @references Y Benjamini, Y Hochberg. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 57, No. 1. (1995), pp. 289-300.
##' @seealso \code{\link{net.comparing}}, \code{\link{comp.subnet}}
##' @export
##' @examples
##'	g<-barabasi.game(100,power=0.8,directed=FALSE)
##'	subg<-induced.subgraph(g,sample(1:100,30))
##'	comp.rand.subnet(subg,g)
##' comp.rand.subnet(subg,g,degree=TRUE)

comp.rand.subnet<-function(subgraph,graph,nsim=1000,degree=FALSE,betweenness=FALSE,
                           ave.path.len=FALSE,eccentricity=FALSE,cc=FALSE,method="utest",
                           FDR=0.05)
{
	if((!is.igraph(subgraph)) | (!is.igraph(graph))){
		stop("subgraph and graph must be igraph object")
	}
	result<-NULL
	##  sample from graph with the same node number as subgraph
	n<-vcount(subgraph)
	graph_sample<-lapply(1:nsim,function(i){induced.subgraph(graph,sample(1:vcount(graph),n))})
  
	if(degree){
		dg_sub<-degree(subgraph)
		dg_sam<-lapply(graph_sample,function(i){degree(i)})
		pvalues<-multi_test(dg_sub,dg_sam)
		dg_sam_mean<-unlist(lapply(dg_sam,mean))
		dg_sub_mean<-mean(dg_sub)
		result<-rbind(result,Degree=c(dg_sub_mean,mean(dg_sam_mean),pvalue(pvalues,FDR=FDR)))
		hist(dg_sam_mean,col="blue",xlab="Mean degree",main="")
		abline(v=dg_sub_mean,col="red")
	}
	if(betweenness){
		bt_sub<-betweenness(subgraph)
		bt_sam<-lapply(graph_sample,function(i){betweenness(i)})
		pvalues<-multi_test(bt_sub,bt_sam)
		bt_sam_mean<-unlist(lapply(bt_sam,mean))
		bt_sub_mean<-mean(bt_sub)
		result<-rbind(result,Betweenness=c(bt_sub_mean,mean(bt_sam_mean),pvalue(pvalues,FDR=FDR)))
		hist(bt_sam_mean,col="blue",xlab="Mean betweenness",main="")
		abline(v=bt_sub_mean,col="red")
	}
	if(ave.path.len){
		apl_sub<-average.path.length(subgraph)
		apl_sam<-lapply(graph_sample,function(i){average.path.length(i)})
		pvalues<-multi_test(apl_sub,apl_sam)
		apl_sam_mean<-unlist(lapply(apl_sam,mean))
		apl_sub_mean<-mean(apl_sub)
		result<-rbind(result,Ave.path.len=c(apl_sub_mean,mean(apl_sam_mean),pvalue(pvalues,FDR=FDR)))
		hist(apl_sam_mean,col="blue",xlab="Mean ave.path.len",main="")
		abline(v=apl_sub_mean,col="red")
	}
	if(eccentricity){
		ec_sub<-eccentricity(subgraph)
		ec_sam<-lapply(graph_sample,function(i){eccentricity(i)})
		pvalues<-multi_test(ec_sub,ec_sam)
		ec_sam_mean<-unlist(lapply(ec_sam,mean))
		ec_sub_mean<-mean(ec_sub)
		result<-rbind(result,Eccentricity=c(ec_sub_mean,mean(ec_sam_mean),pvalue(pvalues,FDR=FDR)))
		hist(ec_sam_mean,col="blue",xlab="Mean eccentricity",main="")
		abline(v=ec_sub_mean,col="red")
	}
	if(cc){
	  cc_sub<-transitivity(subgraph)
	  cc_sam<-lapply(graph_sample,function(i){transitivity(i)})
	  pvalues<-multi_test(cc_sub,cc_sam)
	  cc_sam_mean<-unlist(lapply(cc_sam,mean))
	  cc_sub_mean<-mean(cc_sub)
	  result<-rbind(result,Clust.coeffi=c(cc_sub_mean,mean(cc_sam_mean),pvalue(pvalues,FDR=FDR)))
	  hist(cc_sam_mean,col="blue",xlab="Mean clustering coefficient",main="")
	  abline(v=cc_sub_mean,col="red")
	}
	if(!is.null(result)){
    colnames(result)<-c(substitute(subgraph),substitute(graph),"pvalue")
	}
  
	return(result)
}

##  Multiple test with pvalues output.
multi_test<-function(sample1,sample_list,method="utest"){
	unlist(lapply(sample_list,function(i){wilcox.test(sample1,i)$p.value}))
}

##  Calculate P-value with given FDR.
pvalue<-function(pvalues,FDR=0.1){
	m<-length(pvalues)
	pvalues<-sort(pvalues)
	k<-which((pvalues-(1:m)/m*FDR)<=0)
	if(length(k)){
		k<-k[length(k)]
		return(pvalues[k])
	}else{
		return(pvalues[length(pvalues)])
	}
}
