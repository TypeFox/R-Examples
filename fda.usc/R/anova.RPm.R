anova.RPm=function(object,formula,data.fac,RP=min(30,ncol(object)),alpha=0.95,
hetero=TRUE,pr=FALSE,w=rep(1,ncol(object)),nboot=0,contrast=NULL,...){
  if (is.data.frame(object)) data=as.matrix(object)
  else if (is.fdata(object)) data=object[["data"]]
  if (class(fdata)=="data.fac") data.fac=as.matrix(data.fac) #new
  min.data.fac<-min(table(data.fac))
  if (min.data.fac==0)  stop("Contingency table of factor levels (data.fac argument) contains 0 counts  values")
  nrow=nrow(data);ncol=ncol(data)
  bonf=(1-alpha)/RP
  nprRP=max(RP)
  terms.fd=attr(terms(formula),"term.labels")
  fml=as.formula(paste("value ~ ", paste(terms.fd, collapse= "+")))
  if (is.null(nrow) || is.null(ncol)) stop("fdata must be a matrix")
  nterms=length(terms.fd)+1
  modulo=function(z){sqrt(sum(z^2))}
  z=rnorm(ncol*nprRP)
  z=matrix(z,nrow=nprRP,ncol=ncol)
  z=t(t(z)*w)
  modu=apply(z,1,modulo)
  z=z/modu
  ff=attr(terms.fd,"factors")
  ffcol=colnames(ff)
  if (is.null(contrast)){
    mat=matrix(NA,ncol=(nterms-1),nrow=nprRP)
    colnames(mat)=c(terms.fd)
    ncontrast=0
  }
  else {
     b=length(contrast)
     mf2=model.frame(formula,data.fac)
     uniq=apply(mf2,2,function(x){length(unique(x))})
     ncontrast=rep(0,len=b)
     name.contr=rep(NA,len=sum(ncontrast))
     bb=length(uniq)
     cnombres=ncontrast=rep(0,len=b)
     tnombres=rep(FALSE,len=length(ncontrast))
     tgroups=rep(FALSE,len=(length(ffcol)+1))
if (pr) {print(contrast);print("contrast     contrast  contrast")}
     for (i in 1:b)    {
       a=which(names(contrast)[i]==colnames(mf2))
        tgroups[a]=tnombres[a]=TRUE
         if (is.vector(contrast[[i]])) {
              if (pr) print("Is vector")
              print(a)
                          ncontrast[a]=1
                          contrast[[i]]=matrix(contrast[[i]],ncol=1)
                          }
       else ncontrast[a]=ncol(contrast[[i]])
       names(ncontrast)[a]=names(contrast[i])
                       }
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
     mat=matrix(NA,ncol=(nterms+sum(ncontrast)-1),nrow=nprRP)
     colnames(mat)=c(terms.fd,name.contr)
if (pr) {
   print(ncontrast);print("ncontrast")
   print(contrast);print("contrast")
   print(name.contr);print("name.contr")
   }
  }
  for (j in 1:nprRP){
    value=data%*%z[j,]
    mdata=as.data.frame(cbind(value,data.fac))
    colnames(mdata)=c("value",colnames(data.fac))
    result=aov(fml,data=mdata)
    out=summary(result)
    if (!hetero){
       result3=lm(fml,data=mdata,contrasts=contrast,...)
       if (pr) {print(result);print(result3)}
       if (!is.null(contrast)) {
          out4=summary(result3)
          if (pr) {print(out4);print("summary lm contrast")}
          if (length(out)==1) {
            mat[j,1:(nterms-1)]=out[[1]][1:(nterms-1),5]
            ind=nterms;ind2=2
            for (i in 1:bb)    {
                mat[j,ind:(ind+ncontrast[i]-1)]=out4$coefficients[ind2:(ind2+ncontrast[i]-1),4]
                ind=ind+ncontrast[i];ind2=ind2+ncontrast[i]-1
         }        }
        if (pr) {print(mat);print("mat")}}
    else {  out=summary(result)
           if (length(out)==1) mat[j,1:(nterms-1)]=out[[1]][1:(nterms-1),5] }
    }
  else {              ###############################
      out=anova.hetero(object=mdata,fml,pr=pr,contrast=contrast)[[1]]
      mat[j,]=out[,4]
      if (pr) print(out)
      if (!is.null(contrast)) mat[j,nterms]=out[nterms,4]
    }
  }
  if (pr) {print(summary(mat));print(paste("Bonferroni:",round(bonf,6)))}
  tmat=matrix(NA,nrow=length(RP),ncol=ncol(mat));colnames(tmat)=colnames(mat);rownames(tmat)=paste("RP",RP,sep="")
  mins=matrix(NA,nrow=length(RP),ncol=ncol(mat));colnames(mins)=colnames(mat);rownames(mins)=rownames(tmat)
  tFDR=matrix(NA,nrow=length(RP),ncol=ncol(mat));colnames(tFDR)=colnames(mat);rownames(tFDR)=rownames(tmat)
  pFDR=matrix(NA,nrow=length(RP),ncol=ncol(mat));colnames(pFDR)=colnames(mat);rownames(pFDR)=rownames(tmat)
  for (l in 1:length(RP)){
      if (RP[l]==1)      tmat[l,]=mat[1,]
      else {
      tmat[l,1:(nterms-1+sum(ncontrast))]=
      apply(matrix(mat[1:RP[l],1:(nterms-1+sum(ncontrast))],ncol=(nterms-1+sum(ncontrast))),2,min)}
      if (RP[l]==1) mins[l,]=rep(1,len=ncol(mat))
      else {
      mins[l,]=apply(matrix(mat[1:RP[l],1:(nterms-1+sum(ncontrast))],ncol=(nterms-1+sum(ncontrast))),2,which.min)
      }
      if (RP[l]==1) tFDR[l,]=(mat[1,]<(1-alpha))
      else {
      tFDR[l,1:(nterms-1+sum(ncontrast))]=
      apply(matrix(mat[1:RP[l],1:(nterms-1+sum(ncontrast))],ncol=(nterms-1+sum(ncontrast))),2,FDR,alpha=alpha)}
      if (RP[l]==1) pFDR[l,]=mat[1,]
      else {
      pFDR[l,1:(nterms-1+sum(ncontrast))]=
      apply(matrix(mat[1:RP[l],1:(nterms-1+sum(ncontrast))],,ncol=(nterms-1+sum(ncontrast))),2,pvalue.FDR)}
  }
  tbonf = (tmat < matrix(bonf, nrow = length(RP), ncol = ncol(mat), byrow = TRUE))
  colnames(tbonf) = colnames(mat);rownames(tbonf) = rownames(tmat)
  pbonf=tmat*RP
  colnames(pbonf) = colnames(mat);rownames(pbonf) = rownames(tmat)
  if (is.null(contrast)) {
   resb=matrix(NA,nrow=length(RP),ncol=nterms-1)    }
  else { resb=matrix(NA,nrow=length(RP),ncol=(nterms-1+sum(ncontrast)))}
  if (pr) {print(tbonf);print(pbonf);print(tFDR);print(pFDR)}
  if (nboot>0){
     if (pr) print("bootstrap procedure") ############
     resboot=anova.RPm.boot(data,formula,data.fac,RP=RP,alpha=alpha,nboot=nboot,
     z=z,hetero=hetero,contrast=contrast,pr=pr,...)
     if (pr) print(resboot)
     if (is.null(contrast)){nbuc=nterms-1}
     else {nbuc=nterms-1+sum(ncontrast)}
     for (l in 1:length(RP)){
         for (k in 1:nbuc){
          if ((nterms-1)==1)    Fn=ecdf(resboot[[l]][k])
          else   Fn=ecdf(resboot[[l]][,k])
          resb[l,k]=Fn(tmat[l,k])
          }}
     rownames(resb)=paste("RP",RP,sep="")
  }
pFDR[pFDR>1]=1
pbonf[pbonf>1]=1
res=list("proj"=z,"mins"=mins,"result"=mat,
"test.Bonf"=tbonf,"p.Bonf"=pbonf,"test.FDR"=tFDR,"p.FDR"=pFDR)
if (nboot>0) {
   resb[resb>1]=1
   res$test.Boot=(resb<(1-alpha));res$p.Boot=resb
  colnames(res$p.Boot)=colnames(mat);colnames(res$test.Boot)=colnames(mat)}
if (!is.null(contrast)) res$contrast=contrast
 class(res)="anova.RPm"
return(res)
}

