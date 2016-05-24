theta.it <-
function(bug, sims, max.phi=8, max.psi=1){
  if(is.null(colnames(sims)))
    stop("columns of sims can not be NULL. names should correspond to parameters")
  if(class(bug)!="tsbugs")
    stop("bug must be a object of class tsbugs")
  
  k<-nrow(sims)
  param<-colnames(sims)
  phi<-matrix(0,k,max.phi+1)
  colnames(phi)<-paste0("phi",0:max.phi)
  id<-param[grep("phi", param)]
  if(!is.null(id))
    phi[,which(colnames(phi) %in% param)]<-data.matrix(sims[,sort(id)])
  if(bug$info$variance=="CV"){
    sigma<-NULL
    if("sigma" %in% param)  sigma<-sims[,"sigma"]
    sigma2<-NULL
    if("sigma2" %in% param)  sigma2<-sims[,"sigma2"]
    isigma2<-NULL
    if("isigma2" %in% param)  isigma2<-sims[,"isigma2"]
    temp<-list(phi=phi, sigma=sigma, sigma2=sigma2, isigma2=isigma2)
  }
  if(bug$info$variance=="SV"){
    psi<-matrix(0,k,max.psi+1)
    colnames(psi)<-paste0("psi",0:max.psi)
    id<-param[grep("psi", param)]
    if(length(grep(".star", id))>0) id<-id[-grep(".star", id)]
    psi[,which(colnames(psi) %in% param)]<-data.matrix(sims[,sort(id)])
    psi.star<-NULL
    if(length(grep(".star", param))>0){
      psi.star<-matrix(0,k,max.psi+1)
      colnames(psi.star)<-paste0("psi",0:max.psi,".star")
      id<-param[grep(".star", param)]
      psi.star[,which(colnames(psi.star) %in% param)]<-data.matrix(sims[,sort(id)])
    }
    tau<-NULL
    if("tau" %in% param)  tau<-sims[,"tau"]
    if("itau2" %in% param)  itau2<-sims[,"itau2"]
    id<-param[grep("h.", param,fixed=TRUE)]
    h<-data.matrix(sims[,id])
    temp<-list(phi=phi, psi=psi, psi.star=psi.star, tau=tau, itau2=itau2, h=h)
  }
  if(bug$info$variance=="RV"){
    sig0<-NULL
    if("sig0" %in% param)  sig0<-sims[,"sig0"]
    isig02<-NULL
    if("isig02" %in% param)  isig02<-sims[,"isig02"]
    id<-param[grep("beta.", param,fixed=TRUE)]
    beta<-data.matrix(sims[,id])
    id<-param[grep("delta.", param,fixed=TRUE)]
    delta<-data.matrix(sims[,id])
    epsilon<-NULL
    lambda<-NULL
    ilambda2<-NULL
    if("epsilon" %in% param)  epsilon<-sims[,"epsilon"]
    if("lambda" %in% param)  lambda<-sims[,"lambda"]
    if("ilambda2" %in% param)  ilambda2<-sims[,"ilambda2"]
    temp<-list(phi=phi, sig0=sig0, isig02=isig02, beta=beta, delta=delta, epsilon=epsilon, lambda=lambda, ilambda2=ilambda2)
  }
  temp
}
