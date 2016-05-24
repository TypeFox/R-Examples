cltsk <-
function(query,obs,th,nbins,xcoord='x',ycoord='y',tcoord='t',zcoord='z',
		vth=NULL,vlen=NULL,llim=c(3,3),verbose=T,Large=2000,future=T,cl=NULL)
{
 seed <- round(runif(1) * 1000000)
 l.query <- check_input(query,xcoord,ycoord,tcoord,zcoord)
 l.query <- check_na(l.query[,c(xcoord,ycoord,tcoord)],'query')
 l.obs <- check_input(obs,xcoord,ycoord,tcoord,zcoord)
 l.obs <- check_na(l.obs,'observed')
  
 ## prepare Bins
 dbins <- seq(0,th[1],len=nbins[1]+1)
 tbins <- seq(0,th[2],len=nbins[2]+1)
 bins <- expand.grid(dth = dbins[-1], tth = tbins[-1])
 bins <- as.matrix(bins)
 
 if(is.null(cl)){
	set.seed(seed=seed)
	## Use single core mode
 	out <- apply(l.query,1,working.cltsk,obs=l.obs,th=th,bins=bins,vth=vth,vlen=vlen,
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
  
	ll.args <- list(th=th,bins=bins,vth=vth,vlen=vlen,llim=llim,verbose=verbose,Large=Large,future=future)
	out1 <- clusterMap(cl=cl,fun=working.cltsk.par,ll.query,ll.obs,MoreArgs=ll.args)
	r.order <- do.call(c,ll.order)
	out <- t(do.call(cbind,out1))
	out <- out[order(r.order),]
	
 }
 else{
   stop(cl," is not a cluster\n")
 }
 txt <- paste('D',rep(1:nbins[1],nbins[2]),'T',rep(1:nbins[2],each=nbins[1]),sep='_')
 colnames(out) <- txt
 legend <- data.frame(varName=txt,bins)
 list(krig=out,legend=legend)
}
