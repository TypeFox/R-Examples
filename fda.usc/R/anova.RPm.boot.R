anova.RPm.boot=function(object,formula,data.fac,RP=min(30,ncol(object)),
alpha=0.95,z=NULL,nboot=500,hetero=FALSE,contrast=NULL,pr=FALSE,...){
  if (is.data.frame(object)) object=as.matrix(object)
  else if (is.fdata(object)) object=object[["data"]]
  if (class(object)=="data.fac") data.fac=as.matrix(data.fac) #new
  min.data.fac<-min(table(data.fac))
  if (min.data.fac==0)  stop("Contingency table of factor levels (data.fac argument) contains 0 counts  values")
  nrow=nrow(object);ncol=ncol(object)
 nprRP=max(RP)
 if (!is.null(z) & is.matrix(z)) nprRP=nrow(z)
 if (is.null(z)) {
        modulo=function(z){sqrt(sum(z^2))}
        z=rnorm(ncol*nprRP)
        z=matrix(z,nrow=nprRP,ncol=ncol)
        modu=apply(z,1,modulo)
        z=z/modu}
 terms.fd=attr(terms(formula),"term.labels")
 nterms=length(terms.fd)+1
 lterms=length(terms.fd)
 ff=attr(terms.fd,"factors")
 ffcol=colnames(ff)
 fml=as.formula(paste("object ~ ", paste(terms.fd, collapse= "+")))
 fmlb=as.formula(paste("value ~ ", paste(terms.fd, collapse= "+")))
 if (is.null(contrast))  {
      bb=array(NA,dim=c(nboot,nprRP,nterms-1))
      if (lterms==1)    dimnames(bb)[[3]]=list(terms.fd)
      else     dimnames(bb)[[3]]=(terms.fd) #
    ncontrast=0
 }
 else {
     b=length(contrast)
     mf2=model.frame(formula,data.fac)
     uniq=apply(mf2,2,function(x){length(unique(x))})
     ncontrast=rep(0,len=b)
     name.contr=rep(NA,len=sum(ncontrast))
     bb2=length(uniq)
     cnombres=ncontrast=rep(0,len=b)
     tnombres=rep(FALSE,len=length(ncontrast))
     tgroups=rep(FALSE,len=(length(ffcol)+1))
     for (i in 1:b)    {
       a=which(names(contrast)[i]==colnames(mf2))
        tgroups[a]=tnombres[a]=TRUE
         if (is.vector(contrast[[i]])) {
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
  bb=array(NA,dim=c(nboot,nprRP,(nterms+sum(ncontrast)-1)))
  }
 for (k in 1:ncol(data.fac)) assign(names(data.fac)[k],as.factor(data.fac[,k]))
 fit=manova(fml)
 pb=txtProgressBar(min=0,max=nboot,style=3)
 for (i in 1:nboot){
    setTxtProgressBar(pb,i-0.5)
    if (hetero){
        if (pr) print("hetero")
        term=intersect(attr(terms(formula),"term.labels"),names(data.fac))
        datafac2=as.matrix(data.fac[,term])
        uu=unique(datafac2)
        l=NULL
        for (j in 1:nrow(uu)){
        l2=which(apply(datafac2==rep(uu[j,],rep(nrow(datafac2),ncol(uu))),1,all))
        l=c(l,sample(l2,replace=TRUE))}
        }
    else {
        if (pr) print(" no hetero")
        l=sample(1:nrow(fit$residuals),replace=TRUE)        }
    funcboot=fit$residuals[l,]
    for (j in 1:nprRP){
       value=funcboot%*%z[j,]
       if (hetero){
           mdata=as.data.frame(cbind(value,datafac2)) #####################
           resb=anova.hetero(object=mdata,fmlb,pr=FALSE,contrast=contrast)[[1]]
           if (is.null(contrast)) {bb[i,j,1:(nterms-1)]=resb[,4]}
           else {bb[i,j,]=resb[,4]}
           }
       else {
        mdata=as.data.frame(cbind(value,data.fac))
        colnames(mdata)=c("value",colnames(data.fac))
        result=aov(fmlb,data=mdata)
        out=summary(result)
        if (pr) print(out)
        if (!is.null(contrast)) {
              result3=lm(fmlb,data=mdata,contrasts=contrast,...)
              out4=summary(result3)
              if (pr)     print(out4)
              if (length(out)==1) {
               bb[i,j,1:(nterms-1)]=out[[1]][1:(nterms-1),5] #p-valor
               if (pr) print(out4$coefficients)
               ind=nterms
               ind2=2
              for (ib in 1:bb2)    {
               if (ncontrast[ib]!=0) {
               bb[i,j,ind:(ind+ncontrast[ib]-1)]=out4$coefficients[ind2:(ind2+ncontrast[ib]-1),4]
                }
                ind=ind+ncontrast[ib];ind2=ind2+ncontrast[ib]-1
           }       }   }
           else   bb[i,j,1:(nterms-1)]=out[[1]][1:(nterms-1),5] #new
           }     }
        setTxtProgressBar(pb,i)
 }
 close(pb)
 resboot=vector("list",length(RP))
 for (k in 1:length(RP)){
      if (RP[k]==1)              resboot[[k]]=bb[,1,]
      else {
      if (lterms==1)       resboot[[k]]=apply(bb[,1:RP[k],],c(1),min)
      else   resboot[[k]]=apply(bb[,1:RP[k],],c(1,3),min)
      }
  }
 names(resboot)=paste("RP",RP,sep="")
 return(resboot)
}


