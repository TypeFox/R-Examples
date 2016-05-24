"scalebypanelist" <- function(matrice,center=TRUE,scale=FALSE,col.p,col.j,firstvar,lastvar=ncol(matrice),method="coeff"){

if ((method!="coeff")&(method!="none")&(method!="average")) stop(paste("The method",method,"is unknown. Use coeff or average"))
for (j in 1 :(firstvar-1))  matrice[,j] <- as.factor(matrice[,j])
  labjuge=levels(matrice[,col.j])
  levels(matrice[,col.j])=1:length(labjuge)
  oo <- order(matrice[,col.p])
  matrice <- matrice[oo,]
  oo <- order(matrice[,col.j])
  matrice <- matrice[oo,]
  matrice <- matrice[,c(col.j,col.p,firstvar:lastvar)]

  nbjuge <- nlevels(matrice[,1])
  if (0 %in% summary(matrice[,1])) nbjuge <- nbjuge - 1
  nbprod <- length(levels(matrice[,2]))
  nbdesc <- ncol(matrice)-2

#  #### Calculate the average of the scores by product and by panelist ####
  moy.aux <- as.data.frame(matrix(0,((nbjuge+1)*nbprod),(nbdesc+2)))
  moy.aux[,2]<-as.factor(rep(levels(matrice[,2])[1:nbprod],(nbjuge+1)))
  moy.aux[,1]<-as.factor(rep(0:nbjuge,rep(nbprod,(nbjuge+1))))
  for (juge in 1:nbjuge) moy.aux[(nbprod*juge+1):(nbprod*(juge+1)),3:ncol(moy.aux)]<-sapply(as.data.frame(matrice[matrice[,1]==juge,-(1:2)]),function(vec,fac) tapply(vec,fac,mean,na.rm=TRUE),fac=matrice[matrice[,1]==juge,2])

  if (any(is.na(moy.aux))&(method!="none")){
  #### Missing value are replaced by the average of the average for the descriptor
    if (method=="average") {
      for (i in 3:(nbdesc+2)) moy.aux[,i]<-replace(moy.aux[,i],is.na(moy.aux[,i]),mean(moy.aux[1:nbprod,i]))
    }  
    #### Calculate the average by product and the average by product and by panelist, 
    #### Missing values are replaced by the adjusted coefficient alpha_prod+beta_panelist 
    old.contr = options()$contrasts
    options(contrasts=c("contr.sum", "contr.sum"))
    if (method=="coeff")  {
     coeff.prod <- matrix(0,nbprod,nbdesc)
     coeff.juge <- matrix(0,nbjuge,nbdesc)
     for (i in 1:nbdesc){
       # If all data are missing for a panelist and a descriptor, put a score, why not 0
       for (juge in 1:nbjuge) if (is.na(mean(moy.aux[moy.aux[,1]==juge,i+2],na.rm=TRUE))) moy.aux[which(moy.aux[,1]==juge)[1],i+2]<-0
       res <- summary.lm(aov(moy.aux[(nbprod+1):nrow(moy.aux),i+2]~ as.factor(as.character(moy.aux[(nbprod+1):nrow(moy.aux),2])) + as.factor(as.character(moy.aux[(nbprod+1):nrow(moy.aux),1])),na.action=na.exclude))$coef
       coeff.prod[1:(nbprod-1),i] <- res[1]+res[2:nbprod]
       coeff.prod[nbprod,i] <- res[1]-sum(res[2:nbprod])
       coeff.juge[1:(nbjuge-1),i] <- res[(nbprod+1):(nbjuge+nbprod-1)]
       coeff.juge[nbjuge,i] <- -sum(res[(nbprod+1):(nbjuge+nbprod-1)])
     }
     dimnames(coeff.prod) <- list(moy.aux[1:nbprod,2],colnames(matrice)[3:ncol(matrice)])
     dimnames(coeff.juge) <- list(1:nbjuge,colnames(matrice)[3:ncol(matrice)])
     moy.aux[1:nbprod,3:dim(moy.aux)[2]] <- coeff.prod
     for (i in 1:nbdesc){
      for (k in 1:nbjuge){
       for (kk in 1:nbprod) moy.aux[k*nbprod+kk,i+2]<-replace(moy.aux[k*nbprod+kk,i+2],is.na(moy.aux[k*nbprod+kk,i+2]),coeff.prod[kk,i]+coeff.juge[k,i])
      }
    }
  }}
  moy.aux <- as.data.frame(moy.aux)

  if (center|scale){
  moy.aux.center<-matrix(0,((nbjuge+1)*nbprod),nbdesc+1)
  moy.auxi<-as.matrix(moy.aux[,-2])
  moy.auxi<-matrix(as.double(moy.auxi),nrow=(nbjuge+1)*nbprod,ncol=(nbdesc+1))
  for (i in 1:nbdesc){
    vec<-tapply(moy.auxi[,1+i],moy.auxi[,1],mean,na.rm=TRUE)
    vec<-rep(vec,rep(nbprod,(nbjuge+1)))
    moy.aux.center[,1+i]<-moy.auxi[,1+i]-vec
  }
  moy.aux<-cbind.data.frame(moy.aux[,1:2],moy.aux.center[,-1])
  }

  if (scale){
  #### Calculate the scaled scores by panelist ####
  moy.aux.scale<-matrix(0,((nbjuge+1)*nbprod),(nbdesc+2))
  moy.aux.scale[1:nbprod,3:(nbdesc+2)]<-apply(moy.aux[1:nbprod,3:(nbdesc+2)],2,scale)*sqrt(nbprod/(nbprod-1))
  for (juge in 1:nbjuge){
    for (k in 3:(nbdesc+2)){
      if (var(moy.aux[(moy.aux[,1]==juge),k],na.rm=TRUE)==0) moy.aux.scale[(nbprod+(juge-1)*nbprod)+1:nbprod,k]<-0
      if (var(moy.aux[(moy.aux[,1]==juge),k],na.rm=TRUE)!=0) moy.aux.scale[(nbprod+(juge-1)*nbprod)+1:nbprod,k]<-(moy.aux[(moy.aux[,1]==juge),k]-mean(moy.aux[(moy.aux[,1]==juge),k],na.rm=TRUE))/sqrt((nbprod-1)/nbprod*var(moy.aux[(moy.aux[,1]==juge),k],na.rm=TRUE))
  }}
  moy.aux<-cbind.data.frame(moy.aux[,1:2],moy.aux.scale[,-(1:2)])
  }
  
  for (i in 1:nbprod) moy.aux[i,-(1:2)] <- apply(moy.aux[(moy.aux[,1]!=0)&(moy.aux[,2]==moy.aux[i,2]),-(1:2)],2,mean,na.rm=TRUE)
  dimnames(moy.aux)[2] <- dimnames(matrice)[2]
  levels(moy.aux[,1])=c("0",labjuge)
  return(moy.aux)
  options(contrasts=old.contr)
}
