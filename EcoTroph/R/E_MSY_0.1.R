
E0.1=function(TL,Y,Fish_mort){
  y=Y[as.character(TL),];ff=Fish_mort[as.character(TL),]
  if(length(y[!y==0])==0){f.01=c(NA,NA)}else{
    #pente.origine=(y[2]-y[1])/(ff[2]-ff[1])
    pente=(y[2:length(y)]-y[1:(length(y)-1)])/(ff[2:length(y)]-ff[1:(length(ff)-1)])
    names(pente)=names(y[1:(length(y)-1)])
    dif=abs(pente-0.1*pente[1])# pente[1]: pente à l'origine 
    n=names(dif[dif==min(dif)]) ; if(length(n)>1){n=n[1]}
 f.01=c(ff[names(ff)==as.numeric(n)],n)}
 return(as.numeric(f.01))
}

Emsy=function(TL,Y,Fish_mort){
  y=Y[as.character(TL),];ff=Fish_mort[as.character(TL),]
  x=cbind(y,ff)
  if(length(y[!y==0])==0){emsy=c(NA,NA)}else{
  emsy=c(x[x[,1]==max(y)[1],'ff'][1],names(y[y==max(y)[1]])[1])}
  return(as.numeric(emsy))
}

E_MSY_0.1=function(data,Mul_eff=NULL,B.Input=NULL,Beta=NULL,TopD=NULL,FormD=NULL,TLpred=NULL,maxTL=NULL){
  #data=create.ETmain(ecopath_guinee)
  #Mul_eff=NULL;B.Input=NULL;Beta=NULL;TopD=NULL;FormD=NULL;TLpred=NULL;maxTL=NULL
  n.TL=nrow(data$ET_Main)
  if (is.null(Mul_eff)){Mul_eff=seq(0,10,.1)}
  if (is.null(B.Input)){B.Input <- FALSE}
  if (is.null(Beta)){Beta <- 0.2}
  if (is.null(TopD)){TopD <- rep(.4,n.TL)}else{if(length(TopD)==1){TopD=rep(TopD,n.TL)}}
  if (is.null(FormD)){FormD <- rep(.5,n.TL)}else{if(length(FormD)==1){FormD=rep(FormD,n.TL)}}
  if (is.null(TLpred)){TLpred <- 3.5}
  if (is.null(maxTL)){maxTL=5.5}
  
  diagn.list=create.ETdiagnosis(data=data,Mul_eff=Mul_eff,B.Input=B.Input,Beta=Beta,TopD=TopD,FormD=FormD,TLpred=TLpred,same.mE=T)
  fleet.of.interest=diagn.list[['fleet.of.interest']]
  if(!is.null(fleet.of.interest)){diagn.list=diagn.list[-length(diagn.list)]}
  names(diagn.list)=Mul_eff

  for(i in 1:length(Mul_eff)){
   fm=diagn.list[[as.character(Mul_eff[i])]][['Fish_mort']]
   y=diagn.list[[as.character(Mul_eff[i])]][['Catches.tot']]
  if(i==1){Fish_mort=fm;Y=y}else{Fish_mort=cbind(Fish_mort,fm);Y=cbind(Y,y)}  
  }
  colnames(Y)=Mul_eff;colnames(Fish_mort)=Mul_eff
  
  E_MSY=sapply(row.names(Y),Emsy,Y,Fish_mort)
  E_0.1=sapply(row.names(Y),E0.1,Y,Fish_mort)
  MSY_0.1=cbind(t(E_MSY),t(E_0.1))
  colnames(MSY_0.1)=c('F_MSY','E_MSY','F_0.1','E_0.1')
  MSY_0.1=as.data.frame(MSY_0.1)
  
  # Control, if E_MSY=max(Mul_eff)=> NA
  msy.na=rownames(MSY_0.1[MSY_0.1$E_MSY==max(Mul_eff) & !is.na(MSY_0.1$E_MSY),])
  if(!length(msy.na)==0){MSY_0.1[row.names(MSY_0.1)%in%msy.na,c(1,2)]=matrix(data=NA,ncol=2,nrow=length(msy.na))}
  # Control E_0.1
  m01.na=rownames(MSY_0.1[MSY_0.1$E_0.1==max(Mul_eff) & !is.na(MSY_0.1$E_0.1),])
  if(!length(m01.na)==0){MSY_0.1[row.names(MSY_0.1)%in%m01.na,c(1,2)]=matrix(data=NA,ncol=2,nrow=length(m01.na))}

  MSY_0.1=MSY_0.1[as.numeric(row.names(MSY_0.1))<=maxTL,] 
return(MSY_0.1)}