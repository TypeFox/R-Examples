pintsc <-
function(traits,control=NA){

  if(is.na(control[1]))  stop("Undefined control trait!")
  X<-cbind(traits,control)
  nas<-length(unique(which(is.na(X),arr.ind=T)[,1]))
  if(nas>0)
	{
  warning(paste("Rows containing missing data (",nas, if(nas==1) " row", if(nas>1) " rows",") has been removed to perform the analysis",sep=""))
  X<-na.exclude(X)
	}

  N<-ncol(X)

  traits<-X[,(1:ncol(X))[-N]] # redefino los traits
  c.trait<-X[,N]
  cor_X<-cor.par(traits, c.trait, silent=TRUE)
  eig_X<-eigen(cor_X, only.values=TRUE)$values
  d <- eig_X
  p <- length (d)
  n <- nrow(X)
  INT<-sum((d-1)^2)/(p)
  INT.c<-(INT-((p-1)/n))
  pref="PINTsc = "
  pref2="RelPINTsc = "
  perc<-(INT/(p-1))*100
  pref3="PINTsc.c = "
  pref5="N = "

names<-matrix(c(pref,pref2,pref3,pref5))
outs<-matrix(c(
round(INT, 3),
round(perc, 3),
round(INT.c, 3),
n
))
row.names(outs)<-names 
colnames(outs)<-""

OUT<-as.list(outs)
namesOUT<-c()
for(i in 1:nrow(outs))
namesOUT<-c(namesOUT,strsplit(row.names(outs),split=" =")[[i]][1])
names(OUT)<-namesOUT
OUT<-OUT
OUT
}
