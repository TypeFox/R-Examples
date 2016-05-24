plotBicluster <- function(x,dataset,col=gray(seq(from=1,to=0,length=100))){
	heatmap(dataset[x$features,c(x$samples,setdiff(1:ncol(dataset),x$samples))],col=col,Colv=NA,Rowv=NA,ColSideColors=c(rep("black",length(x$samples)),rep("white",ncol(dataset)-length(x$samples))))
}
