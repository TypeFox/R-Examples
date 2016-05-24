
skda<-function(x, y, tau, method="Bayes") {

nclass=length(unique(y))
if(length(setdiff(y,1:nclass))>0.5)
  {stop("The categorical response must be coded as 1, 2, ..., K!!!")}

y=y-1 # change category coding from 1,2,...,K to 0,1,...,(K-1)
n=length(y)
p=dim(as.matrix(x))[2]

lam=matrix(0, p,1)
phat=matrix(0,n,nclass)
maxct=40;

priorporp=switch(method,
       mle=matrix(1/nclass, nclass,1),
       Bayes=matrix(as.vector(table(y)/length(y)), nclass,1)
)

wu=.C("coordescent", as.integer(n), as.integer(p), as.integer(nclass), as.double(x),
                 as.integer(y), as.double(tau),
                 as.integer(maxct), as.double(priorporp), as.double(lam), as.double(phat));
lam=wu[[9]]
phat=matrix(wu[[10]],n,nclass)
return(list(lam=lam, phat=phat))
}
#####################################################
predprob<-function(x,y,lam,xnew, method="Bayes") {

nclass=length(unique(y))

if(length(setdiff(y,1:nclass))>0.5)
  {stop("The categorical response must be coded as 1, 2, ..., K!!!")}

y=y-1 # change category coding from 1,2,...,K to 0,1,...,(K-1)
n=length(y);
p=dim(as.matrix(x))[2];
m=dim(as.matrix(xnew))[1];


priorporp=switch(method,
       mle=matrix(1/nclass, nclass,1),
       Bayes=matrix(as.vector(table(y)/length(y)), nclass,1)
)


phat=matrix(0,m,nclass);
result=.C("condprob", as.integer(n), as.integer(p),as.integer(nclass),
             as.integer(y), as.double(x), as.double(lam),
              as.integer(m), as.double(xnew), as.double(priorporp),
              as.double(phat));
phat=matrix(result[[10]], m, nclass);
return(phat);
}

########################################################

evalonefold<-function(i,j,n,k,x,y,taus,nclass,method){
     tuind=seq(i,n,k);
     trind=setdiff(1:n,tuind)
     xtr=x[trind,];
     ytr=y[trind];
     xtu=x[tuind,];
     ytu=y[tuind];

      temp=skda(xtr, ytr, taus[j], method);
      lam=temp$lam;
      ptuhat=predprob(xtr,ytr,lam,xtu, method);
 
      fer=0
      for (m in 1:nclass) 
        {ind=which(ytu==m);
         fer=fer+sum(log(ptuhat[ind,m]));
       }
fer
}

###########################
cvskda<-function(x,y,taus,nfolds=10, method="Bayes") {
nclass=length(unique(y))

if(length(setdiff(y,1:nclass))>0.5)
  {stop("The categorical response must be coded as 1, 2, ..., K!!!")}

if(length(taus)<1.5) 
{stop("There must be more than one taus!!!")
}
n=length(y);
tmpord=sample(n);
x=x[tmpord,];
y=y[tmpord];
er=matrix(0,nfolds,length(taus));


 mc = detectCores()
 if(.Platform$OS.type=="windows") {
 cl <- makeCluster(getOption("cl.cores", mc))
 }

er=NULL
for(j in 1:length(taus)) {

if(.Platform$OS.type=="windows")  {
#a=lapply(1:nfolds, function(t)evalonefold(t,j,n,nfolds,x,y,taus,nclass,method))
a=clusterApplyLB (cl, 1:nfolds, function(t)evalonefold(t,j,n,nfolds,x,y,taus,nclass,method))
}
else {
a=mclapply(1:nfolds, function(t)evalonefold(t,j,n,nfolds,x,y,taus,nclass,method), mc.cores = mc)
}
er=cbind(er,unlist(a))

  if((j>1.5)&(sum(er[,j])<sum(er[,j-1])))
    {break;
  }
}
  
  
   if(.Platform$OS.type=="windows") {
   stopCluster(cl)
   }
cver=apply(er[,1:min(j,length(taus))],2,sum);
jbest=which.max(cver);

  temp=skda(x, y, taus[jbest], method);
lam=temp$lam
lam
}
 

