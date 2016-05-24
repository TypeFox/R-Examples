#  Growth transition matrix via Chen et al. 2003

growtrans<-function(Lmin=NULL, Lmax=NULL, Linc=NULL,
             Linf=NULL,SELinf=NULL, K=NULL, SEK=NULL, rhoLinfK=NULL){
   if(is.null(Lmin)) stop ("Lmin is missing") 
   if(is.null(Lmax)) stop ("Lmax is missing") 
   if(is.null(Linc)) stop ("Linc is missing") 
   if(is.null(Linf)) stop ("Linf is missing") 
   if(is.null(SELinf)) stop ("SELinf is missing") 
   if(is.null(K)) stop ("K is missing") 
   if(is.null(SEK)) stop ("SEK is missing") 
   if(is.null(rhoLinfK)) stop ("rhoLinfK is missing") 
   if(Lmin>Lmax) stop ("Lmin is larger than Lmax") 
   Ln<-NULL;COV<-NULL;DL<-NULL;VL<-NULL;growmat<-NULL
   rhoLinfK<-abs(rhoLinfK)
   Ln<-seq(Lmin,Lmax+1,Linc)
   COV<-rhoLinfK*SEK*SELinf
   DL<-(Linf-Ln)*(1-exp(-K))
   DL<-DL[1:(length(which(DL>=0))+1)]
   Ln<-Ln[1:(length(which(DL>=0))+1)]
   DL<-ifelse(DL<0,0,DL)
   VL<-SELinf^2*(1-exp(-K))^2+(Linf-Ln)^2*SEK^2*exp(-2*K)+2*COV*(1-exp(-K))*(Linf-Ln)*exp(-K)
   growmat<-matrix(0,nrow=length(Ln)-1,ncol=length(Ln)-1)
   for(L in 1:as.numeric(length(Ln)-1)){
     for(m in L:as.numeric(length(Ln)-1)){
      growmat[L,m]<-pnorm(Ln[m+1],mean=Ln[L]+DL[L],sqrt(VL[m]))-pnorm(Ln[m],mean=Ln[L]+DL[L],sqrt(VL[m]))
    }
   }
  growmat<-growmat/rowSums(growmat,na.rm=T)
  dimnames(growmat)[[1]]<-Ln[1:as.numeric(length(Ln)-1)]
  dimnames(growmat)[[2]]<-Ln[1:as.numeric(length(Ln)-1)]
  return(growmat)
}


