polarisation.EGR <-
function(alpha,beta,rho,y,f=NULL,dist=NULL,weights=NULL,pm0=NA,lower=NULL,upper=NULL,...){
  if(is.null(weights)) weights<-rep(1/length(y),length(y))
  ER       <- polarisation.ER(alpha,rho)
  if(is.null(f)&is.null(dist)){
    dist.out <- "empirical distribution used"
    P        <- ER$P-beta*(weighted.gini(x=y,w=weights)$Gini-weighted.gini(ER$means,ER$shares)$Gini)
  }else if(!is.null(f)){
    dist.out <- "user specified density used"
    P        <- ER$P-beta*(gini.den(incs=y,dens=f,pm0=pm0,lower=lower,upper=upper)$Gini-weighted.gini(ER$means,ER$shares)$Gini)
  }else if(!is.null(dist)){
    dist.out <- "parametric distribution used"
    ddist    <- get(paste("d", dist, sep = ""), mode = "function", 
                    envir = parent.frame())
    fm<-matrix(NA,length(y),dim(rho)[1])
    for(i in 1:length(y)){
      fm[i,]       <- ddist(y[i],...)
    } 
    fm <- fm%*%diag(rho$shares) 
    f  <- apply(fm,1,sum)
    P        <- ER$P-beta*(gini.den(incs=y,dens=f,pm0=pm0,lower=lower,upper=upper)$Gini-weighted.gini(ER$means,ER$shares)$Gini)
  }
  P  <- as.vector(P)
  PG <- P+beta
  list(P=P,PG=PG,alpha=alpha,beta=beta,dist=dist.out)
}
