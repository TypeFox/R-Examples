# reparticion de la inercia del ACS
# funcion para dividir inercia del ACS
# entra tab, rbl, cbl
inerepa <- function(tab,rbl,cbl)
{
  # control de entrada
  #------------------- pendiente
  # factores fila y columna
  rbl.fac <- factor(rep(names(rbl),rbl,each=TRUE))
  cbl.fac <- factor(rep(names(cbl),cbl,each=TRUE))
  # tabla T suma de bloques
  T <- rowsum(tab,rbl.fac) 
  T <- t(rowsum(t(T),cbl.fac)) 
  #- inter-intra 
  TJ <- rowsum(tab,rbl.fac) 
  #- intra-inter 
  TL <- data.frame(t(rowsum(t(tab),cbl.fac)))
  acs <- dudi.coa(data.frame(tab),scannf=FALSE)
  iner <- NULL
  iner<-c(iner,inT=sum(dudi.coa(data.frame(T),scannf=FALSE)$eig))
  iner<-c(iner,winTL=sum(wca(dudi.coa(data.frame(TL),scannf=FALSE),
                 rbl.fac,scannf=FALSE)$eig))
  iner<-c(iner,winTJ=sum(wca(dudi.coa(data.frame(t(TJ)),scannf=FALSE),
                 cbl.fac,scannf=FALSE)$eig))
  iner<-c(iner,inACI=sum(witwit.model(acs,rbl,cbl,scannf=FALSE)$eig))
  dimen <- dim(T)
  dimen <- rbind(dimen,dim(TL))
  dimen <- rbind(dimen,dim(TJ))
  dimen <- rbind(dimen,dim(tab))
  colnames(dimen)<-c("rows","columns")
  rownames(dimen)<- c("T","TL","TJ","tab")
  inACS <- sum(acs$eig)
  return(list(iner=iner,dimen=dimen,inACS=inACS))
}
# fin funcion  
