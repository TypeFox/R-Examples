#rm(list=ls())
#load("~/lavori/Rdevel/data/Dx.rda")
#load("~/lavori/Rdevel/data/PM.rda")
#Q <- NULL; Pparm <- NULL; fake.model=NULL;p=NULL
#Pparm <- list(p=c(0,0),gam=c(0,4),del=c(0,1.5))

partition.replacement <- function(Dx,PM,Q=NULL,Pparm=NULL,
                                  fake.model=NULL,p=NULL) {
  
  if ((!is.matrix(Dx))&(!is.data.frame(Dx))) stop("Dx must be a matrix or data-frame") 
  
  (LPM <- unique(as.vector(PM))) # numero di partizioni
  if (min(LPM[which(LPM>0)])!=1) stop("Partition numbers must start from one, check partition matrix.")
  
  if (is.null(Q)) Q <- max(Dx,na.rm=TRUE)
  if ((is.null(Pparm))&(is.null(fake.model))) {
    warning("zero replacements")
    rg <- length(LPM)-1
    Pparm <- list(p=matrix(0,rg,ncol(Dx)),
                  gam=matrix(1,rg,ncol(Dx)),del=matrix(1,rg,ncol(Dx)))
  } else {
    if (!is.null(fake.model)) {
      if (length(LPM[which(LPM!=0)])!=length(fake.model)) 
        stop("Length fake.model must be equal to number of partitions: check fake.model and PM")
      if (is.null(p)) {
        warning("zero replacements")
        p <- matrix(0,(length(LPM)-1),ncol(Dx))
      } else {
        if (length(LPM[which(LPM!=0)])!=nrow(p))
          stop("p rows must be equal to number of partitions: check p and PM")
      }
      
      gam <- NULL; del <- NULL
      for (j in 1:length(fake.model)) {
        FMpar <- model.fake.par(fake.model[j])
        gam <- cbind(gam,matrix(FMpar$gam,ncol=1))
        del <- cbind(gam,matrix(FMpar$del,ncol=1))
      }
      Pparm <- list(p=p,gam=gam,del=del)
    }
  }
  
  if (!is.matrix(Pparm$p)) Pparm$p <- matrix(Pparm$p,ncol=2)
  if (!is.matrix(Pparm$gam)) Pparm$gam <- matrix(Pparm$gam,ncol=2)
  if (!is.matrix(Pparm$del)) Pparm$del <- matrix(Pparm$del,ncol=2)
  
  if (length(LPM[which(LPM!=0)])!=nrow(Pparm$p)) 
    stop("Pparm rows must be equal to number of partitions: check Pparm and PM")
  
  Fx <- Dx
  for (h in sort(unique(as.vector(PM)))[-1]) {
    righe <- NULL; colonne <- NULL
    for (i in 1:nrow(PM)) {
      for (j in 1:ncol(PM)) {
        if (PM[i,j]==h) {
          righe <- c(righe,i)
          colonne <- c(colonne,j)
        }
      }
    }
    righe <- unique(righe)
    colonne <- unique(colonne)
    (K <- as.matrix(Dx[righe,colonne],length(righe),length(colonne)))
    
    (R <- replacement.matrix(Q,p=Pparm$p[h,],
                             gam=Pparm$gam[h,],del=Pparm$del[h,]))
    Df <- rdatarepl(K,R,FALSE)$Fx
    
    for (i in 1:length(righe)) {
      for (j in 1:length(colonne)) {
        Fx[righe[i],colonne[j]] <- Df[i,j]
      }
    }
  }
  
  Delta <- Dx-Fx
  Delta[Delta!=0] <- 1
  Fperc <- sum(Delta,na.rm=TRUE)/(prod(dim(Delta)))*100
  cat(paste(round(Fperc,2),"% of data replaced.",sep=""),"\n")
  return(Fx)
}

#R <- matrix(c(1,.3,.3,1),2,2)
#Dx <- rdatagen(n=20,R=R,Q=5)$data

## partition matrix
#PM <- matrix(0,nrow(Dx),ncol(Dx))
#PM[3:10,2] <- 1
#PM[3:10,1] <- 1
#partition.replacement(Dx,PM,Pparm=Pparm) # warning! zero replacements

#partition.replacement(Dx,PM,fake.model="uninformative",p=matrix(c(.3,.2),ncol=2))