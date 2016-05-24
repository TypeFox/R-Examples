ltsk.cv <- function(nfold,obs,th,nbins,part=NULL,zcoord='z',...)
{
  nfold <- min(nfold,nrow(obs))
  
  ## prepare Bins
  dbins <- seq(0,th[1],len=nbins[1]+1)
  tbins <- seq(0,th[2],len=nbins[2]+1)
  bins <- expand.grid(dth = dbins[-1], tth = tbins[-1])
  bins <- as.matrix(bins)
  
  residual <- matrix(NA,nrow(obs),nrow(bins))
  if(is.null(part)){
    part <- sample(1:nfold, nrow(obs), replace = TRUE)
  }
  
  for(i in 1:nfold){
    sel <- (part != i)
    m.model <- obs[sel, ]
    m.valid <- obs[!sel, ]
    tmp<- try(cltsk(query=m.valid,obs=m.model,th=th,nbins=nbins,zcoord=zcoord,...),silent=T)
    if(class(tmp)=="try-error"){
      browser()
      save(m.valid,m.model,sel,file='dump.rda')
      stop(attr(tmp,"condition"))
    }
    else{
      residual[!sel,] <- m.valid[,zcoord]-tmp$krig 
    }
  }
 
  colnames(residual) <- colnames(tmp$krig)
  
  stat <- tmp$legend
  stat$n  <- apply(residual,2,function(v) sum(!is.na(v)))
  stat$SSE  <- apply(residual,2,function(v) sum(v[!is.na(v)]^2))
  stat$MSE  <- with(stat,SSE/n)
  
  list(residual=residual,stat=stat)
}