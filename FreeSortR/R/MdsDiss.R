

################################################################
#__________________    Function  MdsDiss     __________________
#               Mds of a dissimilarity matrix
################################################################

MdsDiss<-function(MatDissimil,ndim=2,metric=TRUE,ties="primary",itmax=5000,eps=1e-06){
  
  #library(smacof)
  # Metric Mds of MatDissimil data
  res<-smacofSym(MatDissimil,ndim=ndim,type="ratio",verbose=FALSE,itmax=itmax,eps=eps)
  Config<-scale(res$conf,center=TRUE,scale=FALSE)  
  # Rotation of the solution
  W<- Config%*%t(Config)                    
  bid<-svd(W)                                               
  Config<-bid$u[,1:ndim]%*%sqrt(diag(bid$d[1:ndim]))                                  
  
  if (metric==FALSE){
    # Non netric Mds of MatDissimil data using metric initialisation
    res<-smacofSym(MatDissimil,ndim=ndim,type="ordinal",ties=ties,init=Config,verbose=FALSE,itmax=itmax,eps=eps)
    Config<-scale(res$conf,center=TRUE,scale=FALSE)  
    # Rotation of the solution
    W<- Config%*%t(Config)                    
    bid<-svd(W)                               
    bid$d[1:ndim]/sum(bid$d)                  
    Config<-bid$u[,1:ndim]%*%sqrt(diag(bid$d[1:ndim])) 
  }
  # Percentage of inertia
  Percent<-bid$d/sum(bid$d)
  #Kruskal Stress
  Stress<-sqrt(sum((res$obsdiss-res$confdiss)^2)/sum(res$confdiss^2))
  
  res<-list(Config=Config,Percent=Percent,Stress=Stress)
  return(res)
}  
