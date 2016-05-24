as.mcmc.blca.gibbs <-
function(x, ...){

#	if(class(fit)[1]!= "blca.gibbs") stop("Invalid Class, must be blca.gibbs object")
	if(inherits(x, "blca.gibbs", TRUE) == 0) stop("Invalid Class, must be blca.gibbs object")
	else{
  G<-length(x$classprob)
  M<-dim(x$itemprob)[2]
  y<-x$samples$classprob
  colnames(y)<-paste(rep("ClassProb",G),1:G)
  for(g in 1:G) y<-cbind(y,x$samples$itemprob[,g,])
  colnames(y)[(G+1):(G+G*M)]<-paste(rep("ItemProb",G*M), rep(1:G, each=M), rep(1:M,G) )
  as.mcmc(y, ...)
  }
}
