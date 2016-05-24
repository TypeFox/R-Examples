descript<-function(part,rec,cl,matsim,iter=FALSE) 
{  
  k<-length(levels(as.factor(part)))
  
  Z<-rec$Z
  X<-rec$X	
  #G<-rec$G
  gc <- rec$g #gravity center
  sdev <- rec$s #standard deviations
  p1<-rec$p1	
  indexj<-rec$indexj 
  n<-rec$n
  p<-rec$p
  indexk<-NULL
  for (i in 1:length(indexj))
    indexk[i]<-part[indexj[i]]
  latent.var <- matrix(,n,k)
  sv<-rep(NA,k)#singular values
  
  coef <- structure(vector(mode = "list", length = k), names = paste("cluster", 1:k, sep = ""))
  var <- structure(vector(mode = "list", length = k), names = paste("cluster", 1:k, sep = ""))
  sim <- structure(vector(mode = "list", length = k), names = paste("cluster", 1:k, sep = ""))
  
  for (g in 1:k)  {
    Zclass<-as.matrix(Z[,which(indexk==g)])
    gclass <- gc[which(indexk==g)]
    sclass <- sdev[which(indexk==g)]
    latent<-clusterscore(Zclass)
    latent.var[,g]<- latent$f
    sv[g]<-latent$sv
    v <- 	latent$v
    clus<-which(part==g)
    #correlation ratio
    C <- matrix(NA, length(clus), 1)
    colnames(C)[1]<-"squared loading"
    rownames(C)<-names(clus)
    for (i in 1:length(clus))
      C[i,1] <- mixedVarSim(latent.var[,g],X[,clus[i]])
    var[[g]] <- C
    #coefficients of the score function
    beta <- matrix(NA,length(v)+1,1)
    beta[1,1] <- -sum(v*gclass/sclass)
    beta[2:(length(v)+1),1] <- v/sclass
    rownames(beta) <- c("const",colnames(Z)[which(indexk==g)])
    coef[[g]] <- beta
    #similarity matrix
    tabcos2 <- matrix(NA, length(clus), length(clus))
    colnames(tabcos2)<-names(clus)
    rownames(tabcos2)<-names(clus)
    diag(tabcos2)<-1
    if ((length(clus)>1) && (matsim==TRUE)) {
      for (i in 1:(length(clus)-1))
        for(j in (i+1):length(clus))  {
          tabcos2[i,j] <- mixedVarSim(X[,clus[i]],X[,clus[j]])
          tabcos2[j,i]<- tabcos2[i,j]
        }
    }		
    if (matsim==TRUE) sim[[g]] <- tabcos2 else sim[[g]] <- NULL
  }
  
  wss<-sv^2 
  size<-rep(NA,k)
  for (g in 1:k) size[g]<-length(which(part==g))
  colnames(latent.var) <- paste("cluster", 1:k, sep = "")
  rownames(latent.var) <-rownames(X)
  Ht<-clusterscore(Z)$sv^2
  E<-(sum(wss)-Ht)/(p-Ht)*100
  return(list(call = cl,var=var,coef=coef,sim=sim,cluster=part,wss=wss,E=E,size=size,scores=latent.var,rec=rec,k=k,iter=iter))
}
