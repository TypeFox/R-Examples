ltsk <-
function(query,obs,th,xcoord='x',ycoord='y',tcoord='t',zcoord='z',
		vth=NULL,vlen=NULL,llim=c(3,3),verbose=T,Large=2000,future=T,cl=NULL)
{
 seed <- round(runif(1) * 1000000)
 l.query <- check_input(query,xcoord,ycoord,tcoord,zcoord)
 l.query <- check_na(l.query[,c(xcoord,ycoord,tcoord)],'query')
 l.obs <- check_input(obs,xcoord,ycoord,tcoord,zcoord)
 l.obs <- check_na(l.obs,'observed')
 
 if(is.null(cl)){
	set.seed(seed=seed)
	## Use single core mode
 	out <- apply(l.query,1,working.ltsk,obs=l.obs,th=th,vth=vth,vlen=vlen,
		llim=llim,verbose=verbose,Large=Large,future=future)
 	out <- t(out)
 }
 else if ("cluster" %in% class(cl)){
	## Use multiple core mode
 	clusterSetRNGStream(cl,seed)
	pwd <- getwd()
 	clusterCall(cl,setwd,dir=pwd)
	clusterEvalQ(cl,library(ltsk))

	res <- partUtil(l.obs,l.query,length(cl),th,xcoord=xcoord,ycoord=ycoord,tcoord=tcoord)
	ll.query <- vector('list',length(res$query))
	ll.obs <- vector('list',length(res$obs))
	ll.order <- vector('list',length(res$query))
	toMat <- function(m,i)
	{
	  if(length(i)>1){
	    mo <- m[i,]
	  }
	  else{
	    mo <- matrix(m[i,],ncol=ncol(m))
	  }
	  mo
	}
	for(i in 1:length(ll.query)){
	  ll.query[[i]] <- toMat(l.query,res$query[[i]])
	  ll.obs[[i]] <- toMat(l.obs,res$obs[[i]])
	  ll.order[[i]] <- res$query[[i]]
	}
	
	ll.args <- list(th=th,vth=vth,vlen=vlen,llim=llim,verbose=verbose,Large=Large,future=future)
	out1 <- clusterMap(cl=cl,fun=working.ltsk.par,ll.query,ll.obs,MoreArgs=ll.args)
	r.order <- do.call(c,ll.order)
	out <- matrix(unlist(out1),ncol=3,byrow=T)
	out <- out[order(r.order),]
	
 }
 else{
   stop(cl," is not a cluster\n")
 }
 colnames(out) <- c('fit','se','flag')
 cbind(query,out)
}
