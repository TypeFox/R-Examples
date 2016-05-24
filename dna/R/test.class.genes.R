test.class.genes=function(X1,X2,genelist=NULL,scores="PLS",distance="abs",
num.permutations=1000,check.networks=TRUE,...){
 X1=as.matrix(X1)
 X2=as.matrix(X2)
 if (check.networks==TRUE){
  PairNetworks=new("pairOfNetworks",network1=X1,network2=X2)
  CommonNetworks=get.common.networks(PairNetworks)
  X1=CommonNetworks$network1
  X2=CommonNetworks$network2
 }
 if (is.null(genelist))
  genelist=1:ncol(X1)
 if (is.character(genelist)){
  genelist=which(colnames(X1)%in%genelist)
 }
 if (is.function(scores)==TRUE)
  R.test.class.genes(X1,X2,genelist,scores,distance,num.permutations,...)
 else if (scores=="PLS")
       PLSnet.test.class.genes(X1,X2,genelist,distance,num.permutations,...)
 else if (scores=="PC")
       PCnet.test.class.genes(X1,X2,genelist,distance,num.permutations,...)
 else if (scores=="cor")
       cornet.test.class.genes(X1,X2,genelist,distance,num.permutations,...)
 else if (scores=="RR")
       RRnet.test.class.genes(X1,X2,genelist,distance,num.permutations,...)
 else
  stop("Error: scores invalid!")
}

PLSnet.test.class.genes=function(X1,X2,genelist,distance="abs",
num.permutations=1000,ncom=3,rescale.data=TRUE,symmetrize.scores=TRUE,
rescale.scores=FALSE){
 if (distance=="abs")
  distancetype=1
 else if (distance=="sqr")
  distancetype=2
 n1=as.integer(nrow(X1))
 n2=as.integer(nrow(X2))
 p=as.integer(ncol(X1))
 n.geneslist=length(genelist)
 out=.C("tdcclassPLS",as.double(X1),as.double(X2),as.integer(genelist),as.integer(n.geneslist),pval=double(1),dlt=double(1),n1,n2,p,as.integer(ncom),as.integer(num.permutations),as.integer(rescale.data),as.integer(symmetrize.scores),as.integer(rescale.scores),as.integer(distancetype),PACKAGE="dna") 
 new("resultsClassTest",p.value=out$pval,delta=out$dlt,class.genes=colnames(X1)[genelist])
}

PCnet.test.class.genes=function(X1,X2,genelist,distance="abs",
num.permutations=1000,ncom=3,rescale.data=TRUE,symmetrize.scores=TRUE,
rescale.scores=FALSE){
 if (distance=="abs")
  distancetype=1
 else if (distance=="sqr")
  distancetype=2
 n1=as.integer(nrow(X1))
 n2=as.integer(nrow(X2))
 p=as.integer(ncol(X1))
 n.geneslist=length(genelist)
 out=.C("tdcclassPC",as.double(X1),as.double(X2),as.integer(genelist),as.integer(n.geneslist),pval=double(1),dlt=double(1),n1,n2,p,as.integer(ncom),as.integer(num.permutations),as.integer(rescale.data),as.integer(symmetrize.scores),as.integer(rescale.scores),as.integer(distancetype),PACKAGE="dna") 
 new("resultsClassTest",p.value=out$pval,delta=out$dlt,class.genes=colnames(X1)[genelist])
}

RRnet.test.class.genes=function(X1,X2,genelist,distance="abs",
num.permutations=1000,lambda=1,rescale.data=TRUE,symmetrize.scores=TRUE,
rescale.scores=FALSE){
 if (distance=="abs")
  distancetype=1
 else if (distance=="sqr")
  distancetype=2
 n1=as.integer(nrow(X1))
 n2=as.integer(nrow(X2))
 p=as.integer(ncol(X1))
 n.geneslist=length(genelist)
 out=.C("tdcclassRR",as.double(X1),as.double(X2),as.integer(genelist),as.integer(n.geneslist),pval=double(1),dlt=double(1),n1,n2,p,as.double(lambda),as.integer(num.permutations),as.integer(rescale.data),as.integer(symmetrize.scores),as.integer(rescale.scores),as.integer(distancetype),PACKAGE="dna") 
 new("resultsClassTest",p.value=out$pval,delta=out$dlt,class.genes=colnames(X1)[genelist])
}

cornet.test.class.genes=function(X1,X2,genelist,distance="abs",
num.permutations=1000,rescale.scores=FALSE){
 if (is.function(distance)==TRUE){
  dist.f=distance 
 }
 else if (distance=="abs")
       dist.f=function(score1,score2){
               abs(score1-score2)
              }
 else if (distance=="sqr")
       dist.f=function(score1,score2){
               (score1-score2)^2
              }
 else
  stop("Error: distance invalid!")

 s1=cornet(X1,rescale.scores)
 s2=cornet(X2,rescale.scores)
 n1=nrow(X1)
 n2=nrow(X2) 

 X=rbind(X1,X2)
 n=n1+n2
 nG=ncol(X)
 n.geneslist=length(genelist)
 temp=dist.f(s1[genelist,genelist],s2[genelist,genelist])
 delta=(sum(temp)-sum(diag(temp)))/(n.geneslist*(n.geneslist-1))
 count.overall=0
 cat("Starting permutation test:\n")
 for (i in 1:num.permutations){
   cat("permutation",i,"out of",num.permutations,"\n")
   i1=as.vector(sample(1:n,n1))
   perm.X1=X[i1,]
   perm.X2=X[-i1,]
   perm.s1=cornet(perm.X1,rescale.scores)
   perm.s2=cornet(perm.X2,rescale.scores)
   temp=dist.f(perm.s1[genelist,genelist],perm.s2[genelist,genelist])
   perm.delta=(sum(temp)-sum(diag(temp)))/(n.geneslist*(n.geneslist-1))
   if (perm.delta>=delta)
    count.overall=count.overall+1
 }

 p.value.overall=count.overall/num.permutations

 new("resultsClassTest",p.value=p.value.overall,delta=delta,class.genes=colnames(X1)[genelist])
}

R.test.class.genes=function(X1,X2,genelist,f,distance="abs",
num.permutations=1000,...){
 if (is.function(distance)==TRUE){
  dist.f=distance 
 }
 else if (distance=="abs")
       dist.f=function(score1,score2){
               abs(score1-score2)
              }
 else if (distance=="sqr")
       dist.f=function(score1,score2){
               (score1-score2)^2
              }
 else
  stop("Error: distance invalid!")

 s1=f(X1,...)
 s2=f(X2,...)
 n1=nrow(X1)
 n2=nrow(X2) 

 X=rbind(X1,X2)
 n=n1+n2
 nG=ncol(X)
 n.geneslist=length(genelist)
 temp=dist.f(s1[genelist,genelist],s2[genelist,genelist])
 delta=(sum(temp)-sum(diag(temp)))/(n.geneslist*(n.geneslist-1))
 count.overall=0
 cat("Starting permutation test:\n")
 for (i in 1:num.permutations){
   cat("permutation",i,"out of",num.permutations,"\n")
   i1=as.vector(sample(1:n,n1))
   perm.X1=X[i1,]
   perm.X2=X[-i1,]
   perm.s1=f(perm.X1,...)
   perm.s2=f(perm.X2,...)
   temp=dist.f(perm.s1[genelist,genelist],perm.s2[genelist,genelist])
   perm.delta=(sum(temp)-sum(diag(temp)))/(n.geneslist*(n.geneslist-1))
   if (perm.delta>=delta)
    count.overall=count.overall+1
 }

 p.value.overall=count.overall/num.permutations

 new("resultsClassTest",p.value=p.value.overall,delta=delta,class.genes=colnames(X1)[genelist])
}
