#to organize data
data.org<-function(x,y,pred,contmed=NULL,binmed=NULL,binref=rep(1,length(binmed)),catmed=NULL,
                   catref=rep(1,length(catmed)),jointm=NULL,biny=T, family1=binomial(link = "logit"),
                   binpred=T,predref=1,alpha=0.1,alpha2=0.1)
{cattobin<-function(x,cat1,cat2=rep(1,length(cat1)))
{ ad1<-function(vec)
{vec1<-vec[-1]
 vec1[vec[1]]<-1
 vec1
}
 dim1<-dim(x)
 catm<-list(n=length(cat1))
 g<-dim1[2]-length(cat1)
 ntemp<-names(x)[cat1]
 j<-1
 for (i in cat1)
 {a<-x[,i]
  d<-rep(0,dim1[1])
  b<-sort(unique(a[a!=cat2[j]]))
  l<-1
  for (k in b)
  {d[a==k]<-l
   l<-l+1}
  d[a==cat2[j]]<-l
  f<-matrix(0,dim1[1],l-1) 
  colnames(f)<-paste(ntemp[j],b,sep=".") #changed for error info
  hi<-d[d!=l & !is.na(d)]
  f[d!=l & !is.na(d),]<-t(apply(cbind(hi,f[d!=l & !is.na(d),]),1,ad1))
  f[is.na(d),]<-NA
  x<-cbind(x,f)
  catm<-append(catm,list((g+1):(g+l-1)))
  g<-g+length(b)
  j<-j+1
 }
 x<-x[,-cat1]
 list(x=x,catm=catm)
}
colnum<-function(vec,cutx) 
{z<-vec
 for (i in 1:length(vec))
   z[i]<-vec[i]-sum(vec[i]>cutx)
 z}

if(!is.null(jointm))  #identify mediators that should be jointly considered
{joint<-NULL 
 for (i in 2:(jointm[[1]]+1))
   joint=c(joint,jointm[[i]])
 joint1<-unique(joint)
 
 if(!is.null(contmed))
 {cont1<-rep(F,length(contmed))
  for (i in 1:length(contmed))
    cont1[i]<-ifelse(sum(contmed[i]==joint1)>0,T,F)
 }
 if(!is.null(binmed))
 {bin1<-rep(F,length(binmed))
  for (i in 1:length(binmed))
    bin1[i]<-ifelse(sum(binmed[i]==joint1)>0,T,F)
 }
 if(!is.null(catmed))
 {cat1<-rep(F,length(catmed))
  for (i in 1:length(catmed))
    cat1[i]<-ifelse(sum(catmed[i]==joint1)>0,T,F)
 }
}
else
{if(!is.null(contmed))
  cont1<-rep(F,length(contmed))
 
 if(!is.null(binmed))
   bin1<-rep(F,length(binmed))
 
 if(!is.null(catmed))
   cat1<-rep(F,length(catmed))
}

if(binpred)
  x[,pred]<-ifelse(x[,pred]==predref,0,1)

if(!is.null(binmed))
{j<-1
 for (i in binmed)
   x[,i]<-ifelse(x[,i]==binref[j],0,1)
 j<-j+1}

if(!is.null(catmed))
{tempx<-cattobin(x,catmed,catref)
 newx<-tempx$x
 catnewx<-tempx$catm}
else newx<-x

#delete variables that are not significant
fullmodel<-summary(glm(y~.,data=data.frame(newx),family=family1))
xname<-names(x)

covr.cont<-rep(F,length(contmed))
if(!is.null(contmed))
{for (i in 1:length(contmed))
  covr.cont[i]<-ifelse(fullmodel$coef[dimnames(fullmodel$coef)[[1]]==xname[contmed[i]],4]<alpha,T,F)
 covr.cont<-ifelse(covr.cont|cont1,T,F)
 cont2<-cont1[covr.cont]
 contmed1<-contmed[covr.cont]
}

covr.bin<-rep(F,length(binmed))
if(!is.null(binmed))
{for (i in 1:length(binmed))
  covr.bin[i]<-ifelse(fullmodel$coef[dimnames(fullmodel$coef)[[1]]==xname[binmed[i]],4]<alpha,T,F)
 covr.bin<-ifelse(covr.bin+bin1>0,T,F) 
 bin2<-bin1[covr.bin]
 binmed1<-binmed[covr.bin]}

covr.cat<-rep(F,length(catmed))
if(!is.null(catmed))
{for (i in 1:length(catmed))
  covr.cat[i]<-ifelse(min(fullmodel$coef[grep(xname[catmed[i]],rownames(fullmodel$coef)),4])<alpha,T,F) #version 6: change to fit names instead of line numbers
 covr.cat<-ifelse(covr.cat+cat1>0,T,F)
 cat2<-cat1[covr.cat]
 catmed1<-catmed[covr.cat]
 catref1<-catref[covr.cat]}

a1<-c(contmed,binmed,catmed)
a2<-c(covr.cont,covr.bin,covr.cat)

cutx<-a1[!a2]

if (sum(a2)==0)
  return ("no mediators found")
else if(length(cutx)==0)
{newx1<-x
 contm1<-contmed
 binm1<-binmed
 catm1<-catmed
 catref1<-catref
 pred1<-pred}
else {newx1<-x[,-cutx]
      if(sum(covr.cont)==0)
      {contm1<-NULL}  
      else 
      {contm1<-colnum(contmed[covr.cont],cutx)}
      if(sum(covr.bin)==0)
        binm1<-NULL
      else 
        binm1<-colnum(binmed[covr.bin],cutx)
      if(sum(covr.cat)==0)
      {catm1<-NULL
       catref1<-NULL}
      else 
      {catm1<-colnum(catmed[covr.cat],cutx)
       catref1<-catref[covr.cat]}
      pred1<-colnum(pred,cutx)
      if(!is.null(jointm))
        for(i in 2:(jointm[[1]]+1))
          jointm[[i]]<-colnum(jointm[[i]],cutx)
}

#delete nonmediators
rela_var<-NULL               #to store the relationships between mediators and predictor   
rela_p<-NULL
name_newx<-names(newx1)
if (binpred)             #for binary (x)
{contm2<-contm1
 if(length(contm1)>0)
 {med.cont<-rep(F,length(contm1))
  for (i in 1:length(contm1))
  {tempmodel<-summary(glm(newx1[,contm1[i]]~newx1[,pred1]))
   med.cont[i]<-ifelse(tempmodel$coef[2,4]<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[contm1[i]])
   rela_p<-c(rela_p,tempmodel$coef[2,4])
  }
  med.cont<-ifelse(med.cont+cont2>0,T,F)
  contm2<-contm1[med.cont]}
 
 binm2<-binm1
 if(length(binm1)>0) 
 {med.bin<-rep(F,length(binm1))
  for (i in 1:length(binm1))   
    {temp.p<-suppressWarnings(chisq.test(newx1[,pred1],newx1[,binm1[i]]))$p.value
     med.bin[i]<-ifelse(temp.p<alpha2,T,F)
     rela_var<-c(rela_var,name_newx[binm1[i]])
     rela_p<-c(rela_p,temp.p)
    }
  med.bin<-ifelse(med.bin+bin2>0,T,F)
  binm2<-binm1[med.bin]}
 
 catm2<-catm1
 if(length(catm1)>0) 
 {med.cat<-rep(F,length(catm1))
  for (i in 1:length(catm1))  
   {temp.p<-suppressWarnings(chisq.test(newx1[,pred1],newx1[,catm1[i]]))$p.value
    med.cat[i]<-ifelse(temp.p<alpha2,T,F)
    rela_var<-c(rela_var,name_newx[catm1[i]])
    rela_p<-c(rela_p,temp.p)
   }
  med.cat<-ifelse(med.cat+cat2>0,T,F)
  catm2<-catm1[med.cat]
  catref2<-catref1[med.cat]}
}
else
{contm2<-contm1
 if(length(contm1)>0)
 {med.cont<-rep(F,length(contm1))
  for (i in 1:length(contm1))
  {temp.p<-cor.test(newx1[,contm1[i]],newx1[,pred1])$p.value
   med.cont[i]<-ifelse(temp.p<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[contm1[i]])
   rela_p<-c(rela_p,temp.p)
  }
 med.cont<-ifelse(med.cont|cont2,T,F)
  contm2<-contm1[med.cont]}
 
 binm2<-binm1
 if(length(binm1)>0) 
 {med.bin<-rep(F,length(binm1))
  for (i in 1:length(binm1))   
  {tempmodel<-summary(glm(newx1[,pred1]~newx1[,binm1[i]]))
   med.bin[i]<-ifelse(tempmodel$coef[2,4]<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[binm1[i]])
   rela_p<-c(rela_p,tempmodel$coef[2,4])
  }    
  med.bin<-ifelse(med.bin|bin2,T,F)
  binm2<-binm1[med.bin]}
 
 catm2<-catm1
 if(length(catm1)>0) 
 {med.cat<-rep(F,length(catm1))
  for (i in 1:length(catm1))   
  {tempmodel<-summary(glm(newx1[,pred1]~newx1[,catm1[i]]))
   med.cat[i]<-ifelse(tempmodel$coef[2,4]<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[catm1[i]])
   rela_p<-c(rela_p,tempmodel$coef[2,4])
  }
  med.cat<-(med.cat|cat2)
  cat3<-cat2[med.cat]
  catm2<-catm1[med.cat]
  catref2<-catref1[med.cat]}
}

if(length(catm2)==0)
  catm2<-NULL
if(length(binm2)==0)
  binm2<-NULL
if(length(contm2)==0)
  contm2<-NULL

newx2<-newx1

if (binpred) 
{if(!is.null(catm2))
  for (i in 1:length(catm2))
    newx2[,catm2[i]]<-as.factor(newx2[,catm2[i]])
 if(!is.null(binm2))
   for (i in 1:length(binm2))
     newx2[,binm2[i]]<-as.factor(newx2[,binm2[i]])
 results<-list(x=newx2, dirx=pred1, contm=contm2, catm=c(binm2,catm2),
              jointm=jointm,fullmodel=fullmodel,rela=data.frame(name=rela_var,p=rela_p))
 class(results)<-"med_iden"
 return(results)
}
else
{catm<-NULL
 if(!is.null(catm2))
 {tempx<-cattobin(newx1,catm2,catref2)
  newx2<-tempx$x
  catm<-tempx$catm
  if(!is.null(binm2))
    binm2<-colnum(binm2,catm2)
  if(!is.null(contm2))
    contm2<-colnum(contm2,catm2)
  pred1<-colnum(pred1,catm2)
  if (!is.null(jointm))
  {if (sum(cat3)==0)
    for(i in 2:(jointm[[1]]+1))
      jointm[[i]]<-colnum(jointm[[i]],catm2)
   else
     for(i in 2:(jointm[[1]]+1))
     {a<-jointm[[i]]
      b<-NULL
      for (j in a)
        if(sum(j==catm2)==0)
          b<-c(b,colnum(j,catm2))
      else for (k in 1:length(catm2))
        if (j==catm2[k])
          b<-c(b,catm[[k+1]])
      jointm[[i]]<-b}
  }
 }
 results<-list(x=newx2,dirx=pred1,contm=contm2,binm=binm2,catm=catm, jointm=jointm, 
               fullmodel=fullmodel,rela=data.frame(name=rela_var,p=rela_p))
 class(results)<-"med_iden"
 return(results)
}
}

summary.med_iden<-function(object,...)
{var.name<-names(object$x)
 mediator<-var.name[c(object$contm,object$binm,object$catm)]
 covariates<-var.name[-c(object$dirx,object$contm,object$binm,object$catm)]
 tests<-object$fullmodel$coefficients
 tests<-cbind(tests[,4],rep(NA,nrow(tests)))
 temp<-dimnames(tests)[[1]]
 for (i in 1:nrow(object$rela))
   tests[grep(object$rela[i,"name"],temp)[1],2]<-object$rela[i,"p"]
 dimnames(tests)[[2]]<-c("P-Value 1", "P-Value 2")
 result<-list(mediator=mediator, covariates=covariates,tests=tests, results=object)
 class(result)<-"summary.med_iden"
 result
}

print.summary.med_iden<-function(x,...)  #version 6: changed typos in the function--tests->x$tests
{cat("Identified as mediators: \n")
 print(x$mediator)
 cat("Selected as covariates: \n")
 print(x$covariates)
 cat("Tests: \n")
 temp<-rownames(x$tests)
 for (z in x$mediator)
   dimnames(x$tests)[[1]][grep(z,temp)]<-paste(temp[grep(z,temp)],"*")
 if (!is.null(x$results$jointm))
  {tt<-NULL
   for (i in 1:(length(x$results$jointm)-1)) 
     tt<-c(tt,x$results$jointm[[i+1]])
   tt<-unique(tt)
   for (z in names(x$results$x)[tt])
     rownames(x$tests)[grep(z,temp)]<-paste(temp[grep(z,temp)],"-")}
 dimnames(x$tests)[[2]]<-c("P-Value in the Full Linear Model", "P-Value Relate to Predictor")
 print(round(x$tests,3))
 cat("----\n *:mediator,-:joint mediator")
}
#for binary predictor
med.binx<-function(x,y,dirx,contm=NULL,catm=NULL,jointm=NULL,allm=c(contm,catm),
                   n=20,seed=sample(1:1000,1),mart=F,nu=0.001,D=3,distn="bernoulli",family1=binomial("logit"))
{if (is.null(allm))
  stop("Error: no potential mediator is specified")
  
  te.binx<-function(full.model,x,y,dirx,best.iter1=NULL)       
{x1<-x[x[,dirx]==1,]
 x0<-x[x[,dirx]==0,]
 n3<-dim(x)[1]
 new1<-x1[sample(1:dim(x1)[1],floor(n3/2),replace=T),]
 new0<-x0[sample(1:dim(x0)[1],floor(n3/2),replace=T),]
 te<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
 te
}

med.binx.contm<-function(full.model,x,y,med,dirx,best.iter1=NULL)  
{x1<-x[x[,dirx]==1,]                               
 x0<-x[x[,dirx]==0,]
 n3<-dim(x)[1]
 marg.m<-c(x1[sample(1:dim(x1)[1],floor(n3/2),replace=T),med],x0[sample(1:dim(x0)[1],floor(n3/2),replace=T),med])
 marg.m<-marg.m[sample(2*floor(n3/2))]
 nom1<-x1[sample(1:dim(x1)[1],floor(n3/2),replace=T),]
 nom0<-x0[sample(1:dim(x0)[1],floor(n3/2),replace=T),]
 new1<-nom1
 new1[,med]<-marg.m[1:floor(n3/2)]
 new0<-nom0
 new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]
 dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
 dir.nom
}

med.binx.jointm<-function(full.model,x,y,med,dirx,best.iter1=NULL)  
{x1<-x[x[,dirx]==1,]                               
 x0<-x[x[,dirx]==0,]
 n3<-dim(x)[1]
 if (length(med)==1)                       #added for the new program, when there is only one mediator
   {sample.temp<-sample(c(sample(which(x[,dirx]==1),floor(n3/2),replace=T),
                   sample(which(x[,dirx]==0),floor(n3/2),replace=T)))
    marg.m<-x[sample.temp,med]}        #added for the new program
 else                                         #added for the new program
   {marg.m<-rbind(x1[sample(1:dim(x1)[1],floor(n3/2),replace=T),med],x0[sample(1:dim(x0)[1],floor(n3/2),replace=T),med])
    marg.m<-marg.m[sample(2*floor(n3/2)),]  }     
 nom1<-x1[sample(1:dim(x1)[1],floor(n3/2),replace=T),]
 nom0<-x0[sample(1:dim(x0)[1],floor(n3/2),replace=T),]
 new1<-nom1
 new0<-nom0
 if(length(med)==1)                                       #added for the new program, when there is only one mediator
   {new1[,med]<-marg.m[1:floor(n3/2)]                     #added for the new program 
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]}  #added for the new program
 else                                                     #added for the new program
   {new1[,med]<-marg.m[1:floor(n3/2),]
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2)),]}
 dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
 dir.nom
}

med.binx.catm<-function(full.model,x,y,med,dirx,best.iter1=NULL)  
{x1<-x[x[,dirx]==1,]                               
 x0<-x[x[,dirx]==0,]
 n3<-dim(x)[1]
 marg.m1<-x1[sample(dim(x1)[1],floor(n3/2),replace=T),med]
 marg.m2<-x0[sample(dim(x0)[1],floor(n3/2),replace=T),med]
 nom1<-x1[sample(1:dim(x1)[1],floor(n3/2),replace=T),]
 nom0<-x0[sample(1:dim(x0)[1],floor(n3/2),replace=T),]
 dir.nom<-0
 for (i in levels(x[,med]))
 {new1<-nom1
  new1[1:dim(new1)[1],med]<-i
  new0<-nom0
  new0[1:dim(new0)[1],med]<-i
  p<-0.5*mean(marg.m1==i,na.rm=T)+ 0.5*mean(marg.m2==i,na.rm=T)
  dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T))
 }
 dir.nom
}

  #1.fit the model
  if (mart)
  {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                       distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
   best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))
   while(full.model$n.trees-best.iter1<30){
     full.model<-gbm.more(full.model, 50)           # do another 50 iterations
     best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}}
  else
  {full.model<-glm(y~., data=x, family=family1)
   best.iter1=NULL}
  #2. prepare for the store of results
  set.seed(seed)
  te<-rep(0,n)
  if(!is.null(jointm))
  {denm<-matrix(0,n,1+length(c(contm,catm))+jointm[[1]])
   dimnames(denm)[[2]]<-c("de",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))}
  else
  {denm<-matrix(0,n,1+length(c(contm,catm)))
   dimnames(denm)[[2]]<-c("de",names(x)[c(contm,catm)])}
  ie<-denm
  #3. repeat to get the mediation effect
  for (k in 1:n)
  {#3.1 get the te         full.model,x,y,dirx,best.iter1=NULL
    te[k]<-te.binx(full.model,x,y,dirx,best.iter1)
    denm[k,1]<-med.binx.jointm(full.model,x,y,allm,dirx,best.iter1)
    j<-2
    #3.2 mediation effect from the continuous mediator
    if (!is.null(contm))
      for (i in contm)          #full.model,x,y,med,dirx,best.iter1=NULL
      {denm[k,j]<-med.binx.contm(full.model,x,y,i,dirx,best.iter1)
       j<-j+1}
    #3.3.mediation effect from the categorical mediator
    if (!is.null(catm))
      for (i in catm)           #full.model,x,y,med,dirx,best.iter1=NULL
      {denm[k,j]<-med.binx.catm(full.model,x,y,i,dirx,best.iter1)
       j<-j+1}
    #3.4 mediation effect from the joint mediators
    if (!is.null(jointm))
      for (i in 1:jointm[[1]])          #full.model,x,y,med,dirx,best.iter1=NULL
      {denm[k,j]<-med.binx.jointm(full.model,x,y,jointm[[i+1]],dirx,best.iter1)
       j<-j+1}
    
    #3.5 get the indirect effects
    ie[k,]<-te[k]-denm[k,]
  }
  a<-list(denm=denm,ie=ie,te=te,model=full.model,best.iter=best.iter1)
  return(a)
}

boot.med.binx<-function(x,y,dirx,contm=NULL,catm=NULL,jointm=NULL,
                        n=20,seed=sample(1:1000,1),n2=50,mart=F,nu=0.001,
                        D=3,distn="bernoulli",family1=binomial("logit"),
                        weight=rep(1,length(y)))
  #n2 is the time of bootstrap
{if (is.null(c(contm,catm)))
  stop("Error: no potential mediator is specified")
 allm=c(contm,catm)
 te<-rep(0,n2+1)
 de<-rep(0,n2+1)
 if(is.null(jointm))
 {ie<-matrix(0,n2+1,1+length(c(contm,catm)))
  dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)])}
 else
 {ie<-matrix(0,n2+1,1+length(c(contm,catm))+jointm[[1]])
  dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))}
 temp<-med.binx(x,y,dirx,contm,catm,jointm,allm,n,seed,mart,nu,D,distn,family1)
 te[1]<-mean(temp$te)
 de[1]<-mean(temp$denm[,1])
 ie[1,]<-apply(temp$ie,2,mean)  #first row is the estimated value
 model<-temp$model
 best.iter<-temp$best.iter
 for (i in 1:n2)
 {boots<-sample(1:length(y),replace=T,prob=weight)
  x1<-x[boots,]
  y1<-y[boots]
  temp<-med.binx(x1,y1,dirx,contm,catm,jointm,allm,n,seed+i,mart,nu,D,distn,family1)
  te[1+i]<-mean(temp$te)
  de[1+i]<-mean(temp$denm[,1])
  ie[1+i,]<-apply(temp$ie,2,mean)  #first row is the estimated value
 }
 a<-list(estimation=list(ie=ie[1,],te=te[1],de=de[1]),bootsresults=list(ie=ie[-1,],te=te[-1],de=de[-1]),model=list(MART=mart,model=model,best.iter=best.iter), 
         data=list(x=x,y=y,dirx=dirx,contm=contm,catm=catm,jointm=jointm,binpred=T))
 class(a)<-"mma"
 return(a)
}

#for continous predictor
med.contx<-function(x,y,dirx,binm=NULL,contm=NULL,catm=NULL, jointm=NULL, margin=1, n=20,seed=sample(1:1000,1),
                    mart=F,nu=0.001,D=3,distn="gaussian",family1=gaussian(link="identity"))
{if (is.null(c(binm,contm,catm)))
  stop("Error: no potential mediator is specified")
  
  anymissing<-function(vec)
{if(sum(is.na(vec))>0)
  return(F)
 else return(T)
}

dist.m.given.x<-function(x,dirx,binm=NULL,contm=NULL,catm=NULL)
{models<-NULL
 res<-NULL
 if(!is.null(catm))
 {for (i in 2:(catm$n+1))
   binm<-c(binm,catm[[i]])}
 
 z<-x[,dirx]
 j<-1
 if(!is.null(binm))
 {for(i in binm)
 {models[[j]]<-glm(x[,i]~z,family=binomial(link = "logit"))
  res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))
  j<-j+1}
 }
 for (i in contm)
 {models[[j]]<-glm(x[,i]~z,family=gaussian(link="identity"))
  res<-cbind(res,models[[j]]$res)
  j<-j+1
 }
 list(models=models,varmat=var(res))
}


sim.xm<-function(distmgivenx,x1,dirx,binm,contm,catm)
{mult.norm<-function(mu,vari,n) 
{  if (nrow(vari)!=ncol(vari)) 
{result<-c("Error: Variance matrix is not square")}  
else if (length(mu)!=nrow(vari)) 
{result<-c("Error: length mu is not right!")}  
else {   p<-length(mu)
         tmp1<-eigen(vari)$values
         tmp2<-eigen(vari)$vectors   
         result<-matrix(0,n,p)   
         for (i in 1:p)
         {result[,i]<-rnorm(n,mean=0,sd=sqrt(tmp1[i]))}   
         for (i in 1:n)
         {result[i,]<-tmp2%*%result[i,]+mu}
}  
result
}

match.margin<-function(vec)   
{range1<-vec[1:2]
 vec1<-vec[-(1:2)]
 range2<-range(vec1)
 vec1<-range1[1]+diff(range1)/diff(range2)*(vec1-range2[1])
 vec1
}

gen.mult<-function(vec)
{l<-1-sum(vec)
 l<-ifelse(l<0,0,l)
 rmultinom(1,size=1,prob=c(l,vec))[-1]
}

means<-NULL
z<-x1[,dirx]
binm1<-binm
if(!is.null(catm))
{for (i in 2:(catm$n+1))
  binm1<-c(binm1,catm[[i]])}

if(!is.null(binm1))
  for (i in 1:length(binm1))
    means<-cbind(means,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(z=z)))

if(!is.null(contm))
  for (i in (length(binm1)+1):length(c(binm1,contm)))
    means<-cbind(means,predict(distmgivenx$models[[i]],newdata=data.frame(z=z)))

if(dim(means)[2]==1)                                                   #added in the new program, in case there is only one mediator
  {sim.m<-rnorm(length(means),mean=means,sd=sqrt(distmgivenx$varmat))     #added in the new program
   sim.m2<-match.margin(c(range(means),sim.m))}                          #added in the new program   
else{
sim.m<-t(apply(means,1,mult.norm,vari=distmgivenx$varmat,n=1))

range.means<-apply(means,2,range)

sim.m2<-apply(rbind(range.means,sim.m),2,match.margin)    #to make the simulate fit the means' ranges
}
n<-dim(sim.m2)[1]
if(!is.null(binm))
  for (i in 1:length(binm))
    sim.m2[,i]<-rbinom(n,size=1,prob=sim.m2[,i])

if(!is.null(catm))
{j<-length(binm)+1
 for (i in 2:(catm$n+1))
 {a<-sim.m2[,j:(j+length(catm[[i]])-1)]
  sim.m2[,j:(j+length(catm[[i]])-1)]<-t(apply(a,1,gen.mult))
  j<-j+length(catm[[i]])}
}

x1[,c(binm1,contm)]<-sim.m2

x1
}

  
  if(is.null(catm))
  multi=jointm
  else if(is.null(jointm))
    multi=catm
  else {temp1<-catm
        temp2<-jointm
        temp1[[1]]=catm[[1]]+jointm[[1]]
        temp2[[1]]<-NULL
        multi=append(temp1,temp2)} 
 listm=list(single=c(contm,binm),multi=multi)
 
 nonmissing<-apply(cbind(y,x[,c(dirx,listm$single)]),1,anymissing)
  x<-x[nonmissing,]
  y<-y[nonmissing]
 
  #1.fit the model
  if(mart)
  {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                       distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
   best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))         
   while(full.model$n.trees-best.iter1<30){
     full.model<-gbm.more(full.model, 50)           # do another 50 iterations
     best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}
  }
  else
    {full.model<-glm(y~., data=x, family=family1)
     best.iter1=NULL}
  
  #2. prepare for the store of results
  set.seed(seed)
  te<-rep(0,n)
  mul<-ifelse(is.null(listm$multi),0,listm$multi[[1]])        #added in the new program, in case multi is null
  denm<-matrix(0,n,1+length(listm$single)+mul)                #added in the new program, in case multi is null
  if(!is.null(listm$multi))
    colnames(denm)<-c("de",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep=""))
  else 
    colnames(denm)<-c("de",names(x)[listm$single])
  ie<-denm
  
  #3. get the joint distribution of m given x
  distmgivenx<-dist.m.given.x(x,dirx,binm,contm,catm) 
  x1<-x
  x1[,dirx]<-x[,dirx]+margin
  ybar0<-mean(predict(full.model,x,best.iter1),na.rm=T)
  
  n1<-dim(x)[1]
  
  denm[1:n,1]<-mean(predict(full.model,x1,best.iter1),na.rm=T)-ybar0   #assume can change x without changing others
  
  #4. repeat to get the mediation effect
  for (k in 1:n)
  {new1<-sim.xm(distmgivenx,x1,dirx,binm,contm,catm)  #draw from the conditional distribution of m given x
   #4.1 get the te
   te[k]<-mean(predict(full.model,new1,best.iter1),na.rm=T)-ybar0
   
   j<-2
   #4.2 mediation effect from the single mediator
   if (!is.null(listm$single))
     for (i in listm$single)
     {new1.nm<-new1
      new1.nm[,listm$single[j-1]]<-x[sample(n1,replace=T),listm$single[j-1]]    #draw m from its original distribution
      denm[k,j]<-mean(predict(full.model,new1.nm,best.iter1),na.rm=T)-ybar0
      j<-j+1}
   
   #4.3.mediation effect from the joint mediator
   if (!is.null(listm$multi))
     for (i in 2:(listm$multi[[1]]+1))
     {new1.nm<-new1
      new1.nm[,listm$multi[[i]]]<-x[sample(n1,replace=T),listm$multi[[i]]]    #draw joint m from its original distribution
      denm[k,j]<-mean(predict(full.model,new1.nm,best.iter1),na.rm=T)-ybar0
      j<-j+1}
   
   #4.4 get the indirect effects
   ie[k,]<-te[k]-denm[k,]
  } 
  a<-list(denm=denm,ie=ie,te=te,model=list(MART=mart,model=full.model,best.iter=best.iter1))
  return(a)
}


boot.med.contx<-function(x,y,dirx,binm=NULL,contm=NULL,catm=NULL, 
                         jointm=NULL, margin=1, n=20,seed=sample(1:1000,1),
                         mart=F,nu=0.001,D=3,distn="gaussian",
                         family1=gaussian(link="identity"),n2=50,
                         weight=rep(1,length(y)))
{if (is.null(c(binm,contm,catm)))
  stop("Error: no potential mediator is specified")
  
 if(is.null(catm))
  {multi=jointm
   name1<-NULL                       #added in the new program
   if (!is.null(multi))              #added in the new program, in case that multi is NULL
   name1<-paste("j",1:multi[[1]],sep="")}
 else if(is.null(jointm))
  {multi=catm
   name1<-NULL
   for (i in 2:(catm[[1]]+1))
   name1<-c(name1,names(x)[multi[[i]][1]])}
 else {temp1<-catm
       temp2<-jointm
       temp1[[1]]=catm[[1]]+jointm[[1]]
       temp2[[1]]<-NULL
       multi=append(temp1,temp2)
       name1<-NULL
       for (i in 2:(catm[[1]]+1))
            name1<-c(name1,names(x)[multi[[i]][1]])
       name1<-c(name1,paste("j",1:jointm[[1]],sep=""))} 
 listm=list(single=c(contm,binm),multi=multi)
 
 te<-rep(0,n2+1)
 de<-rep(0,n2+1)
 mul<-ifelse(is.null(multi),0,multi[[1]])        #added in the new program, in case multi is null
 ie<-matrix(0,n2+1,1+length(listm$single)+mul)   #added in the new program
 if(!is.null(listm$multi))
   dimnames(ie)[[2]]<-c("all",names(x)[listm$single],name1)
 else 
   dimnames(ie)[[2]]<-c("all",names(x)[listm$single])
 
 temp<-med.contx(x,y,dirx,binm,contm,catm,jointm, margin,n,seed,mart,nu,D,distn,family1)
 
 te[1]<-mean(temp$te)
 de[1]<-mean(temp$denm[,1]) 
 ie[1,]<-apply(temp$ie,2,mean)  #first row is the estimated value
 model<-temp$model
 for (l in 1:n2)
 {boots<-sample(1:length(y),replace=T, prob=weight)
  x1<-x[boots,]
  y1<-y[boots]
  temp<-med.contx(x1,y1,dirx,binm,contm,catm,jointm, margin,n,seed+l,mart,nu,D,distn,family1) #added to the new codel, change the seed to make different results
  te[1+l]<-mean(temp$te)
  de[1+l]<-mean(temp$denm[,1])
  ie[1+l,]<-apply(temp$ie,2,mean)  #first row is the estimated value
 }
 a<-list(estimation=list(ie=ie[1,],te=te[1],de=de[1]),bootsresults=list(ie=ie[-1,],te=te[-1],de=de[-1]),model=model,
         data=list(x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm, jointm=jointm, binpred=F))
 class(a)<-"mma"
 return(a)
}


mma<-function(x,y,pred,contmed=NULL,binmed=NULL,binref=rep(1,length(binmed)),catmed=NULL,catref=rep(1,length(catmed)),jointm=NULL,biny=T,
              binpred=T,predref=1,alpha=0.1,alpha2=0.1, margin=1, n=20,seed=sample(1:1000,1),
              mart=F,nu=0.001,D=3,distn="bernoulli",family1=binomial(link = "logit"),n2=50,weight=rep(1,length(y)))
{if(binpred) 
  {data.bin<-data.org(x,y,pred,contmed,binmed,binref,catmed,catref,jointm,biny,family1,binpred=T,predref,alpha,alpha2)
   result<-boot.med.binx(x=data.bin$x,y,dirx=data.bin$dirx,contm=data.bin$contm,catm=data.bin$catm,jointm=data.bin$jointm,n,seed,n2,mart,nu,D,distn,family1,weight)
  }
 else
  {data.cont<-data.org(x,y,pred,contmed,binmed,binref,catmed,catref,jointm,biny,family1,binpred=F,predref=NULL,alpha,alpha2)
   result<-boot.med.contx(x=data.cont$x,y,dirx=data.cont$dirx,binm=data.cont$binm,contm=data.cont$contm,catm=data.cont$catm,jointm=data.cont$jointm,margin, n,seed,
                        mart,nu,D,distn,family1,n2,weight)
  }
 result
}


#classes and methods for mma
print.mma<-function(x,...)
{cat("MMA Analysis: Estimated Mediation Effects Using ")
 if (x$model$MART)
   cat ("MART\n")
 else cat("GLM\n")
 print(x$e)
}


summary.mma<-function(object,...,alpha=0.05,plot=TRUE)
{x<-object
 temp1<-x$boots
 temp2<-x$est
 a1<-alpha/2
 a2<-1-a1
 b1<-qnorm(a1)
 b2<-qnorm(a2)
 temp1.result<-list(indirect.effect=rbind(est=temp2$ie,mean=apply(temp1$ie,2,mean),sd=apply(temp1$ie,2,sd),
                                          upbd=apply(temp1$ie,2,mean)+b2*apply(temp1$ie,2,sd),
                                          lwbd=apply(temp1$ie,2,mean)+b1*apply(temp1$ie,2,sd),
                                          upbd_q=apply(temp1$ie,2,quantile,a2), lwbd_q=apply(temp1$ie,2,quantile,a1)),
                    total.effect=c(est=temp2$te,mean=mean(temp1$te),sd=sd(temp1$te),
                                   upbd=mean(temp1$te)+b2*sd(temp1$te),lwbd=mean(temp1$te)+b1*sd(temp1$te),
                                   upbd=quantile(temp1$te,a2),lwbd=quantile(temp1$te,a1)),
                    direct.effect=c(est=temp2$de,mean=mean(temp1$de),sd=sd(temp1$de),
                                    upbd=mean(temp1$de)+b2*sd(temp1$de),lwbd=mean(temp1$de)+b1*sd(temp1$de),
                                    upbd=quantile(temp1$de,a2),lwbd=quantile(temp1$de,a1)))
 result<-list(results=temp1.result,alpha=alpha,plot=plot,obj=x)
 class(result)<-"summary.mma"
 result
 }

print.summary.mma<-function(x,...)
{cat("MMA Analysis: Estimated Mediation Effects Using ")
 if (x$obj$model$MART)
  cat ("MART\n")
 else cat("GLM\n")
print(lapply(x$results,round,3))

if(x$plot)
{re<-c(x$results$indirect.effec[1,],x$results$dir[1])/x$results$tot[1]
 d<-order(re)
 name1<-c(names(x$results$indirect.effec[1,]),"de")[d]
 barplot(re[d],horiz=T,xlab="Relative Effects",names=name1,
        cex.names=0.9,beside=F,cex.axis=0.9,las=1,xlim=range(re),
        col = rainbow(length(d), start = 3/6, end = 4/6))}
}


plot.mma<-function(x,...,vari,xlim=range(data$x[,vari],na.rm=T))
{marg.den<-function(x,y)
{y<-y[!is.na(x)]
 x<-x[!is.na(x)]
 z1<-unique(x)
 z2<-rep(0,length(z1))
 for (i in 1:length(z1))
   z2[i]<-mean(y[x==z1[i]],na.rm=T)
 z3<-order(z1)
 cbind(z1[z3],z2[z3])
}
op <- par(no.readonly = TRUE) # the whole list of settable par's.
 if (x$model[1]==T) 
 {full.model=x$model$model
  best.iter=x$model[3]
  data=x$data
  mname<-ifelse(is.character(vari),vari,names(data$x)[vari])
  if(data$binpred)
  {if(!is.factor(data$x[,vari]))
   {par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
    if(full.model$distribution=="gaussian")
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim,type="response"))
    for (j in unique(data$x[,data$dirx]))
    hist(data$x[data$x[,data$dirx]==j,vari],xlab="",xlim=xlim, freq=F,main=paste(names(data$x)[data$dirx],j,sep="="))}
   else{par(mfrow=c(2,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
        if(full.model$distribution=="gaussian")
          suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter))
        else
          suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,type="response"))
     plot(data$x[,vari]~data$x[,data$dirx],ylab=mname,xlab=names(data$x)[data$dirx])}
  }
  else
  {par(mfrow=c(2,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
    if(!is.factor(data$x[,vari]))
   {if(full.model$distribution=="gaussian")
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim,type="response"))
    axis(1,at=data$x[,vari],labels=F)
    a<-marg.den(data$x[,data$dirx],data$x[,vari])
    plot(a[,2:1],ylab=names(data$x)[data$dirx],xlim=xlim,xlab="")
    lo <- loess(a[,2]~a[,1])
    lines(lo$fitted[order(lo$x)], lo$x[order(lo$x)],lwd=1)
    }
    
    #lo <- loess(data$x[,data$dirx]~data$x[,vari])
    #plot(data$x[,vari],data$x[,data$dirx],xlim=xlim,ylab=names(data$x)[data$dirx],xlab="")
    #lines(predict(lo), col='red', lwd=2)}
   else
   {if(full.model$distribution=="gaussian")
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,type="response"))
    plot(data$x[,vari],data$x[,data$dirx],ylab=names(data$x)[data$dirx],xlab="")}
  }
}
else
{full.model=x$model$model
 data=x$data
 mname<-ifelse(is.character(vari),vari,names(data$x)[vari])
 if(data$binpred)
 {if(!is.factor(data$x[,vari]))
 {par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
  if(is.null(full.model$na.action))
    data1<-full.model$data[,vari]
  else
    data1<-full.model$data[-full.model$na.action,vari]
  b<-marg.den(data1,full.model$fitted)
  plot(b,xlab=mname,ylab=paste("f(",mname,")",sep=""),xlim=xlim)
  lo1 <- loess(b[,2]~b[,1])
  lines(lo1$x[order(lo1$x)], lo1$fitted[order(lo1$x)], lwd=1)
  if(is.null(full.model$na.action))
    axis(1,at=full.model$data[,vari],labels=F)
  else
    axis(1,at=full.model$data[-full.model$na.action,vari],labels=F)
  for (j in unique(data$x[,data$dirx]))
    hist(data$x[data$x[,data$dirx]==j,vari],xlab="",xlim=xlim, freq=F,main=paste(names(data$x)[data$dirx],j,sep="="))}
 else{par(mfrow=c(2,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
      if(is.null(full.model$na.action))
        plot(full.model$y~full.model$data[,vari],ylab=paste("f(",mname,")",sep=""),xlab=mname)
      else
        plot(full.model$y~full.model$data[-full.model$na.action,vari],ylab=paste("f(",mname,")",sep=""),xlab=mname)
      plot(data$x[,vari]~data$x[,data$dirx],ylab=mname,xlab=names(data$x)[data$dirx])}
 }
 else
 {par(mfrow=c(2,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
  if(!is.factor(data$x[,vari]))
  {if(is.null(full.model$na.action))
    b<-marg.den(full.model$data[,vari],full.model$fitted)
   else
    b<-marg.den(full.model$data[-full.model$na.action,vari],full.model$fitted)
   plot(b,xlab=mname,ylab=paste("f(",mname,")",sep=""),xlim=xlim)
   lo1 <- loess(b[,2]~b[,1])
   lines(lo1$x[order(lo1$x)], lo1$fitted[order(lo1$x)], lwd=1)
   if(is.null(full.model$na.action))
     axis(1,at=full.model$data[,vari],labels=F)
   else
     axis(1,at=full.model$data[-full.model$na.action,vari],labels=F)
   a<-marg.den(data$x[,data$dirx],data$x[,vari])
   plot(a[,2:1],ylab=names(data$x)[data$dirx],xlim=xlim,xlab="")
   lo <- loess(a[,2]~a[,1])
   lines(lo$fitted[order(lo$x)], lo$x[order(lo$x)],lwd=1)}
  else
  {if(is.null(full.model$na.action))
     plot(full.model$y~full.model$data[,vari],ylab=paste("f(",mname,")",sep=""),xlab=mname)
   else
     plot(full.model$y~full.model$data[-full.model$na.action,vari],ylab=paste("f(",mname,")",sep=""),xlab=mname)
    plot(data$x[,vari],data$x[,data$dirx],ylab=names(data$x)[data$dirx],xlab="")}
 }
}
par(op)
}
  

  
  
  
  
