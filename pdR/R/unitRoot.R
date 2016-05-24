## This code computes panel unit root based on 
## Chang(2002)'s simple IGF


pIGF <- function(datamat,maxp,ic,spec){
  N=ncol(datamat)
  IV_stat=NULL
  for (j in 1:N){
    igf=IGF(datamat[,j],maxp,ic,spec)$tstat.IGF
    IV_stat=rbind(IV_stat,igf)
  }
  Ptstat=sum(IV_stat)/sqrt(N)
  Pval=pnorm(Ptstat)
  results=list(panel.tstat=Ptstat,pvalue=Pval)
  return(results)  
  
}


IGF <- function(y,maxp,ic,spec){
  p=lagSelect(y,maxp,ic)
  tmp=y
  DY=as.matrix(embed(diff(tmp),p+1)[,-1])
  dyNAMES=paste("dy",1:ncol(DY),sep="")
  colnames(DY)=dyNAMES
  DEP=embed(tmp,p+2)[,1:2]
  tb=tbar(DEP[,2])
  y=DEP[,1]-tb
  y1=DEP[,2]-tb
  datz=cbind(y,y1,DY)
  Y=datz[,1]
  
  yy=diff(Y)
  CC=3*(1/sqrt(length(yy)))*(1/sqrt(var(yy)));
  fy=y1*exp(-CC*abs(y1))
  trend=1:nrow(fy)
  if (spec==0){
  z=cbind(y1,DY)
  w = cbind(fy,DY)}
  if (spec==1){
  z=cbind(y1,DY,1)
  w = cbind(fy,DY,1)}
  if (spec==2){
  z=cbind(y1,DY,1,trend)
  w = cbind(fy,DY,1,trend)}
  
  iwz = solve(t(w)%*%z)
  beta0 = iwz%*%t(w)%*%Y
  e=Y-z%*%beta0
  sd=sqrt(t(e)%*%e/(length(Y)-p-1))
  va = iwz%*%(t(w)%*%w)%*%t(iwz)
  sa = sd%*%sqrt(va[1,1])
  
  tstat=(beta0[1]-1)/sa
  results=list(tstat.IGF=tstat,beta=beta0,sdev=sa,cV=CC,p=p)
  return(results)
}


lagSelect <-function(y,maxp,ic){
  tmp=y;
  IC=NULL
  for (j in 1:maxp){
    DY=as.matrix(embed(diff(tmp),j+1)[,-1])
    if (ncol(DY)==1) {dyNAMES="dy"}
    dyNAMES=paste("dy",1:ncol(DY),sep="")
    colnames(DY)=dyNAMES
    DEP=embed(tmp,j+2)
    y=DEP[,1];y1=DEP[,2]
    datz=cbind(y,y1,DY)
    eq=as.formula(paste("y~y1+",paste(dyNAMES,collapse="+"),spe=""))
    
    if(ic=="AIC") {ic.tmp=AIC(lm(eq,data=as.data.frame(datz)))}
    if(ic=="BIC") {ic.tmp=BIC(lm(eq,data=as.data.frame(datz)))}
    
    IC=rbind(IC,ic.tmp)
  }
  rownames(IC)=NULL
  return(which(IC==min(IC)))
}


tbar <- function(x){
  T=length(x)
  pmean=NULL
  for (i in 1:T){
    pmean=rbind(pmean,mean(x[1:i]))
  }
  return(pmean)
}
