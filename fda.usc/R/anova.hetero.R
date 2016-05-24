anova.hetero=function(object=NULL,formula,pr=FALSE,contrast=NULL,...){
if (pr) print("INI ANOVA.HETERO")
  if (is.data.frame(object)) data=data.frame(object)
  else if (is.fdata(object)) data=data.frame(object[["data"]])
  nrow=nrow(data);ncol=ncol(data)
  fdata2=intercambio(data) #para evitrar transponer
  Terms <- if (missing(fdata2)) terms(formula, "Error")
  else terms(formula, "Error", fdata2 = fdata2)
  indError <- attr(Terms, "specials")$Error
   if (length(indError) > 1) stop(sprintf(ngettext(length(indError),
   "there are %d Error terms: only 1 is allowed",
  "there are %d Error terms: only 1 is allowed"), length(indError)), domain = NA)
  mf=model.frame(formula,fdata2)
  mfo=model.frame(formula,data) ####
  mt=model.matrix(formula,mf)
  ff=attr(Terms,"factors")
  ffcol=colnames(ff)
  ffrow=rownames(ff)
  ff=as.matrix(ff[-1,],nrow=nrow(ff[-1,]))
  rownames(ff)=ffrow[-1]
  colnames(ff)=ffcol
  nombres=attr(Terms,"term.labels")
  uniq=apply(mf,2,function(x){length(unique(x))})
  lis=sapply(mf,is.factor)
  fact=as.vector(lis)
#  fact=uniq<10
  nfactors= length(uniq[fact])
  nterms=length(Terms)
  aaa=dim(mf)[2]-1
  dddd=(dim(mf)[2]-1)
#  if (nfactors!=(dim(mf)[2]-1)) stop("Not Yet Implemented the use of
# covariables or number of levels great than 10")
 if (nfactors!=(dim(mf)[2]-1)) stop("Not Yet Implemented the use of
 covariables")
  if (is.null(contrast)) {
     ans=matrix(NA,ncol=4,nrow=length(nombres))
     rownames(ans)=nombres}
  else {
            contrast2=intercambio.l(contrast)
            b=length(contrast)
            cnombres=ncontrast=rep(0,len=b)
            tnombres=rep(FALSE,len=length(ncontrast))
            tgroups=rep(FALSE,len=(length(ffcol)+1))
            bb=length(uniq)
            for (i in 1:b)    {
                       a=which(colnames(mfo)==names(contrast)[i])
                       tgroups[a]=tnombres[a]=TRUE
                       if (is.vector(contrast[[i]])) {
                          ncontrast[a]=1
                          contrast[[i]]=matrix(contrast[[i]],ncol=1)
                          }
                       else ncontrast[a]=ncol(contrast[[i]])
                       names(ncontrast)[a]=names(contrast[i])
                       }
            name.contr=rep(NA,len=sum(ncontrast))
            j=1;ji=1;jk=1
            for (i in 1:length(ncontrast))    {
            if (tgroups[ji])  {
                name.contr[j:(j+ncontrast[i]-1)]=paste("C",j:(j+ncontrast[i]-1),".",
                names(contrast[jk]),sep="")
                colnames(contrast[[jk]])=name.contr[j:(j+ncontrast[i]-1)]
              j=j+ncontrast[i];jk=jk+1
              }
              ji=ji+1
              }
    nombres.esp=c(nombres,name.contr)
    ans=matrix(NA,nrow=length(nombres.esp),ncol=4)
    rownames(ans)=nombres.esp          }
  colnames(ans)=c("Est.","df1","df2","p.value")
  nlev=uniq[fact]
  if (pr) print(paste("Levels:",nlev))
  mf2=intercambio(mf)   #para evitar transponer
  n=table(mf2[,fact])
  m=tapply(mf2[,1],mf2[,fact],mean)
  sig=tapply(mf2[,1],mf2[,fact],var)
  sig=sig*(n-1)/n^2
  max.vecm=length(m)
  vecm=matrix(as.vector((m)),nrow=1)
  ssig=as.vector((sig))
  if (pr) {print("Means");print(m)}
  if (pr) {print("Variances");print(sig)}
  vecDelta=as.vector(1/(n-1))
  Delta=diag(vecDelta)
  SN=sum(n)*diag(ssig)
  if (pr) {print("Delta");print(Delta);print("SN");print(SN)}
  for (i in 1:length(nombres)){
      MM=1
      for (k in 1:nrow(ff)){
          J=matrix(1,ncol=nlev[k],nrow=nlev[k])/nlev[k]
          if (ff[k,i]==1) P=diag(1,nlev[k])-J else P=J
          MM=kronecker(MM,P)      #
          }
      DM=diag(diag(MM))
      if (pr) {print(paste("MM",nombres[i]));print(MM)}
      FNM=sum(n)*vecm%*%MM%*%t(vecm)/traza(DM%*%SN)
      f1=traza(DM%*%SN)^2/traza(MM%*%SN%*%MM%*%SN)
      f0=traza(DM%*%SN)^2/traza(DM%*%DM%*%SN%*%SN%*%Delta)
      pF=1-pf(FNM,f1,f0)
      ans[i,]=c(FNM,f1,f0,pF)
  }
  if (!is.null(contrast)){
jj2=ind.f=1
for (i in 1:b)    {
#ind.g=which(colnames(mf)==names(ncontrast[(i+1)]))#ok
ind.g=NA
ind.g=which(colnames(mfo)==names(contrast[i])) #ko
 if (is.vector(contrast[[i]])) vvf=matrix(contrast[[i]],ncol=1)
 else { if (is.list(contrast))    vvf=contrast[[i]]
        else  vvf=contrast}
for (jj in 1:ncol(vvf)) {
                   aa=mfo[,ind.g]
                   n.f=table(aa)
                   m.f=tapply(mfo[,1],mfo[,ind.g],mean)
                   sig.f=tapply(mfo[,1],mfo[,ind.g],var)
                   sig.f=sig.f*(n.f-1)/n.f^2
  max.vecm.f=length(m.f)
  vecm.f=matrix(as.vector(m.f),nrow=1)
  ssig.f=as.vector(sig.f)
#  MM.f=vvf%*%ginv(t(vvf)%*%vvf)%*%t(vvf)
  MM.f=vvf[,jj]%*%ginv(t(vvf[,jj])%*%vvf[,jj])%*%t(vvf[,jj])
  DM.f=diag(diag(MM.f))
  vecDelta.f=as.vector(1/(n.f-1))
  Delta.f=diag(vecDelta.f)
  SN.f=sum(n.f)*diag(ssig.f)
      FNM.f=sum(n.f)*vecm.f%*%MM.f%*%t(vecm.f)/traza(DM.f%*%SN.f)
      f1.f=traza(DM.f%*%SN.f)^2/traza(MM.f%*%SN.f%*%MM.f%*%SN.f)
      f0.f=traza(DM.f%*%SN.f)^2/traza(DM.f%*%DM.f%*%SN.f%*%SN.f%*%Delta.f)
      pF.f=1-pf(FNM.f,f1.f,f0.f)
      ans[length(nombres)+jj2,]=c(FNM.f,f1.f,f0.f,pF.f)
      jj2=jj2+1
}}}
res=list("ans"=ans)
if (!is.null(contrast)) res$contrast=contrast
if (pr) print(res)
class(res)="anova.hetero"
return(res)
}
intercambio.l=function(contrast){
   n=length(contrast)
   mdata=contrast
   avance=seq(1,n)
   retro=seq(n,1)
   mdata[avance]=contrast[retro]
   names(mdata)[avance]= names(contrast)[retro]
return((mdata))
}
intercambio=function(fdata2){
   n=ncol(fdata2)
   mdata=fdata2
   avance=seq(2,n)
   retro=seq(n,2)
   mdata[,avance]=fdata2[,retro]
  colnames(mdata)[avance]= colnames(fdata2)[retro]
return(data.frame(mdata))
}


