lk <- function(query,obs,th,xcoord='x',ycoord='y',zcoord='z',vlen=15,cl=NULL)
{
  l.query <- check_input(query,xcoord,ycoord,'tt__tt',zcoord)
  l.obs <- check_input(obs,xcoord,ycoord,'tt__tt',zcoord)
  seed <- round(runif(1) * 1000000)
  
  if(is.null(cl)){
    r <- working.lk.par(l.query,l.obs,th,xcoord=xcoord,ycoord=ycoord,zcoord=zcoord,vlen=vlen)
  }
  else if ("cluster" %in% class(cl)){
    clusterSetRNGStream(cl,seed)
    pwd <- getwd()
    clusterCall(cl,setwd,dir=pwd)
    clusterEvalQ(cl,library(ltsk))
    
    res <- partSpUtil(l.obs,l.query,length(cl),th,xcoord=xcoord,ycoord=ycoord)
    ll.query <- vector('list',length(res$query))
    ll.obs <- vector('list',length(res$obs))
    ll.order <- vector('list',length(res$query))
    for(i in 1:length(ll.query)){
      ll.query[[i]] <- l.query[res$query[[i]],]
      ll.obs[[i]] <- l.obs[res$obs[[i]],]
      ll.order[[i]] <- res$query[[i]]
    }
    ll.args <- list(th=th,xcoord=xcoord,ycoord=ycoord,zcoord=zcoord,vlen=vlen)
    out1 <- clusterMap(cl=cl,fun=working.lk.par,ll.query,ll.obs,MoreArgs=ll.args)
    r <- do.call(rbind,out1)
    r.order <- do.call(c,ll.order)
    r <- r[order(r.order),]
  
  }
  else{
    stop(cl," is not a cluster\n")
  }
  r
}
