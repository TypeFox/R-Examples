pint <-
function(traits){
 
  X<-traits
  nas<-length(unique(which(is.na(X),arr.ind=T)[,1]))
  if(nas>0)
	{
  warning(paste("Rows containing missing data (",nas, if(nas==1) " row", if(nas>1) " rows",") has been removed to perform the analysis",sep=""))
  X<-na.exclude(traits)
	}

  cor_X <-cor(X)
  eig_X<-eigen(cor_X, only.values=TRUE)$values
  d <- eig_X
  p <- length (d)
  n <- nrow(X)
  INT<-sum((d-1)^2)/(p)
  INT.c<-(INT-((p-1)/n))
  pref="PINT = "
  pref2="RelPINT = "
  pref3="PINT.c = "
  pref5="N = "

names<-matrix(c(pref,pref2,pref3,pref5))
outs<-matrix(c(round(INT, 3),round((INT/(p-1))*100, 3),round(INT.c, 3),n))
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

