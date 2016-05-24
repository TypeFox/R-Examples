genet.dist<-function(dat,diploid=TRUE,method="Dch"){
  cl<-match.call()
  if (is.genind(dat)) dat<-genind2hierfstat(dat)
  if(!is.na(pmatch(method,"Dch")))
    method<-"Dch"
  METHODS<-c("Dch","Da","Ds","Fst","Dm","Dr","Cp","X2","Nei87","WC84")
  method<-pmatch(method,METHODS)
  if (is.na(method)) stop("Invalid distance method")
  if (method==-1) stop("Ambiguous distance method")

  if (method==9) {gdist<-pairwise.neifst(dat,diploid); return(as.dist(gdist))}
  if (method==10) {gdist<-pairwise.WCfst(dat,diploid);return(as.dist(gdist))}
  
    
  pft<-pop.freq(dat,diploid)
  nl<-length(pft)
  npop<-dim(pft[[1]])[2]
  dist.loc<-array(numeric(npop^2*nl),c(npop,npop,nl))
  if (method==1)    temp<-lapply(pft,function(x) 2*(1-t(x)^0.5%*%x^0.5))
  if (method==2) temp<-lapply(pft,function(x) t(x)^0.5%*%x^0.5)
  if (method %in% c(3,4,5)) temp<-lapply(pft,function(x) t(x)%*%x)
  if (method==7) temp<-lapply(pft,function(x) 
    matrix(apply
           (apply(x,1,
                  function(b) outer(b,b,
                                    FUN=function(y,z) abs(y-z)))
            ,1,sum,na.rm=TRUE),nrow=npop))
  if (method==6) temp<-lapply(pft,function(x) 
    matrix(apply
           (apply(x,1,
                  function(b) outer(b,b,
                                    FUN=function(y,z) (y-z)^2))
            ,1,sum,na.rm=TRUE)/2,nrow=npop)^.5)
  if (method==8) temp<-lapply(pft,function(x) 
    matrix(apply
           (apply(x,1,
                  function(b) outer(b,b,
                                    FUN=function(y,z) (y-z)^2/(y+z)))
            ,1,sum,na.rm=TRUE),nrow=npop))

  for (il in 1:nl){dist.loc[,,il]<-temp[[il]]}
  rm(temp)
  nlpp<-apply(dist.loc,c(1,2),function(x) sum(!is.na(x)))

    if (method==1) gdist<-2/pi/nlpp*apply(dist.loc,c(1,2),sum,na.rm=TRUE)
  if (method==2) gdist<-1-1/nlpp*apply(dist.loc,c(1,2),sum,na.rm=TRUE)
  if (method %in% c(3,4,5)) {
    Jxy<-1/nlpp*apply(dist.loc,c(1,2),sum,na.rm=TRUE)
  }
  if (method %in% c(4,5)){
    dum<-apply(dist.loc,3,function(x) outer(diag(x),diag(x),FUN=function(y,z) (y+z)/2))
    dum<-matrix(1/nlpp*apply(dum,1,sum,na.rm=TRUE),nrow=npop) 
  }
  if (method==3) {
    dum<-outer(diag(Jxy),diag(Jxy),FUN="*")
    gdist<- -log(Jxy/dum^.5);
  }
  if (method==4) gdist<-(dum-Jxy)/(1-Jxy)
  if (method==5) gdist<-(dum-Jxy)
  if (method==6) gdist<-1/nlpp*apply(dist.loc,c(1,2),sum,na.rm=TRUE)
  if (method==7) gdist<-1/2/nlpp*apply(dist.loc,c(1,2),sum,na.rm=TRUE)
  if (method==8) gdist<-2/nlpp*apply(dist.loc,c(1,2),sum,na.rm=TRUE)
  return(as.dist(gdist))      
}
