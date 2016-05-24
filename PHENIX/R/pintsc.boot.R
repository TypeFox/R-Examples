pintsc.boot <-
function (traits,control=NA,replicates=1000){
  INT = list()
  INTC = list()

  X<-traits

  if(is.na(control[1]))  stop("Undefined control trait!")
  X<-cbind(traits,control)
  nas<-length(unique(which(is.na(X),arr.ind=T)[,1]))
 if(nas>0)
	{
  warning(paste("Rows containing missing data (",nas, if(nas==1) " row", if(nas>1) " rows",") has been removed to perform the analysis",sep=""))
  X<-na.exclude(X)
	}

  Y<-replicates
  length (INT) = Y

  Z<-ncol(X)

  for (i in 1:Y){
    t.sample<-X[sample(nrow(X), replace=TRUE),]
    traits<-t.sample[,(1:ncol(X))[-Z]]
    c.trait<-t.sample[,Z]
    cor_X<-cor.par(traits, c.trait,silent=TRUE)
    d <-eigen(cor_X, only.values=TRUE)$values
    p <- length (d)

    n <- nrow(t.sample) 
    Int<-sum((d-1)^2)/(p)
    INT[i]<-Int
    Int.c<-(Int-((p-1)/n))
    INTC[i]<-Int.c

	if(i==1) cat("\nStarting bootstrap...........\n")
	if(i==round(Y/4)) cat("\nPerforming bootstrap......25%\n")
	if(i==round(Y/2)) cat("\nPerforming bootstrap......50%\n")
	if(i==round(3*Y/4)) cat("\nPerforming bootstrap......75%\n")
	if(i==Y) cat("\nBootstrap finished.......100%\n")

  }
  Intphen1 <-as.numeric(INT)
  Intphen2 <-as.numeric(INTC)
  pref0="Mean = "
  pref1="Median ="
  pref2="SD = "
  pref3="SE = "
  se1<-(sd(Intphen1)/sqrt(nrow(X)))
  se2<-(sd(Intphen2)/sqrt(nrow(X)))
  pref4="Lower IC 99% = "
  pref5="Higher IC 99% = "
  pref6="Lower IC 95% = "
  pref7="Higher IC 95% = "
 pref8="Number of replicates = "


#Igual que antes:
names<-matrix(c(pref0,pref1,pref2,pref3,pref4,pref5,pref6,pref7,pref8))
outs<-cbind(
c(
round(mean(Intphen1), 3),
round(median(Intphen1), 3),
round(sd(Intphen1), 3),
round(se1, 3),
round(quantile(Intphen1, probs=0.5/100), 3),
round(quantile(Intphen1, probs=99.5/100), 3),
round(quantile(Intphen1, probs=2.5/100), 3),
round(quantile(Intphen1, probs=97.5/100), 3),
length(INT)
)
,
c(
round(mean(Intphen2), 3),
round(median(Intphen2), 3),
round(sd(Intphen2), 3),
round(se2, 3),
round(quantile(Intphen2, probs=0.5/100), 3),
round(quantile(Intphen2, probs=99.5/100), 3),
round(quantile(Intphen2, probs=2.5/100), 3),
round(quantile(Intphen2, probs=97.5/100), 3),
length(INT)
)
)
row.names(outs)<-names
colnames(outs)<-c("PINTSC","PINTSC.C")
outs

}
