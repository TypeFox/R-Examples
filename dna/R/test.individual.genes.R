test.individual.genes=function(X1,X2,scores="PLS",distance="abs",
num.permutations=1000,check.networks=TRUE,...){
 X1=as.matrix(X1)
 X2=as.matrix(X2)
 if (check.networks==TRUE){
  PairNetworks=new("pairOfNetworks",network1=X1,network2=X2)
  CommonNetworks=get.common.networks(PairNetworks)
  X1=CommonNetworks$network1
  X2=CommonNetworks$network2
 }
 if (is.function(scores)==TRUE)
  R.test.individual.genes(X1,X2,scores,distance,num.permutations,...)
 else if (scores=="PLS")
       PLSnet.test.individual.genes(X1,X2,distance,num.permutations,...)
 else if (scores=="PC")
       PCnet.test.individual.genes(X1,X2,distance,num.permutations,...)
 else if (scores=="cor")
       cornet.test.individual.genes(X1,X2,distance,num.permutations,...)
 else if (scores=="RR")
       RRnet.test.individual.genes(X1,X2,distance,num.permutations,...)
 else
  stop("Error: scores invalid!")
}

PLSnet.test.individual.genes=function(X1,X2,distance="abs",
num.permutations=1000,ncom=3,rescale.data=TRUE,symmetrize.scores=TRUE,
rescale.scores=FALSE){
 if (distance=="abs")
  distancetype=1
 else if (distance=="sqr")
  distancetype=2
 n1=as.integer(nrow(X1))
 n2=as.integer(nrow(X2))
 p=as.integer(ncol(X1))
 gene.names=colnames(X1)
 out=.C("tdcindPLS",as.double(X1),as.double(X2),pval=double(p),d=double(p),n1,n2,p,as.integer(ncom),as.integer(num.permutations),as.integer(rescale.data),
as.integer(symmetrize.scores),as.integer(rescale.scores),
as.integer(distancetype),PACKAGE="dna") 
 names(out$pval)=gene.names
 names(out$d)=gene.names
 new("resultsIndTest",p.values=out$pval,d=out$d)
}

PCnet.test.individual.genes=function(X1,X2,distance="abs",
num.permutations=1000,ncom=3,rescale.data=TRUE,symmetrize.scores=TRUE,
rescale.scores=FALSE){
 if (distance=="abs")
  distancetype=1
 else if (distance=="sqr")
  distancetype=2
 n1=as.integer(nrow(X1))
 n2=as.integer(nrow(X2))
 p=as.integer(ncol(X1))
 gene.names=colnames(X1)
 out=.C("tdcindPC",as.double(X1),as.double(X2),pval=double(p),d=double(p),n1,n2,
p,as.integer(ncom),as.integer(num.permutations),as.integer(rescale.data),
as.integer(symmetrize.scores),as.integer(rescale.scores),
as.integer(distancetype),PACKAGE="dna") 
 names(out$pval)=gene.names
 names(out$d)=gene.names
 new("resultsIndTest",p.values=out$pval,d=out$d)
}

RRnet.test.individual.genes=function(X1,X2,distance="abs",
num.permutations=1000,lambda=1,rescale.data=TRUE,symmetrize.scores=TRUE,
rescale.scores=FALSE){
 if (distance=="abs")
  distancetype=1
 else if (distance=="sqr")
  distancetype=2
 n1=as.integer(nrow(X1))
 n2=as.integer(nrow(X2))
 p=as.integer(ncol(X1))
 gene.names=colnames(X1)
 out=.C("tdcindRR",as.double(X1),as.double(X2),pval=double(p),d=double(p),n1,n2,
p,as.double(lambda),as.integer(num.permutations),as.integer(rescale.data),
as.integer(symmetrize.scores),as.integer(rescale.scores),
as.integer(distancetype),PACKAGE="dna") 
 names(out$pval)=gene.names
 names(out$d)=gene.names
 new("resultsIndTest",p.values=out$pval,d=out$d)
}

cornet.test.individual.genes=function(X1,X2,distance="abs",
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
 di=rep(0,nG)
 for (g in 1:nG)
  di[g]=sum(dist.f(s1[g,],s2[g,]))
 di=di/(nG-1)

 perm.di=rep(0,nG)
 count.ind=rep(0,nG)
 cat("Starting permutation test:\n")
 for (i in 1:num.permutations){
  cat("permutation",i,"out of",num.permutations,"\n")
  i1=as.vector(sample(1:n,n1))
  perm.X1=X[i1,]
  perm.X2=X[-i1,]
  perm.s1=cornet(perm.X1,rescale.scores)
  perm.s2=cornet(perm.X2,rescale.scores)
  for (g in 1:nG)
   perm.di[g]=sum(dist.f(perm.s1[g,],perm.s2[g,]))
  perm.di=perm.di/(nG-1)
  for (g in 1:nG)
  if (perm.di[g]>=di[g])
    count.ind[g]=count.ind[g]+1
 }

 p.value.ind=rep(0,nG)
 for (g in 1:nG){
   p.value.ind[g]=count.ind[g]/num.permutations
 }
 names(p.value.ind)=colnames(X1)
 names(di)=colnames(X1)
 new("resultsIndTest",p.values=p.value.ind,d=di)
}


R.test.individual.genes=
function(X1,X2,f,distance="abs",num.permutations=1000,...){
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
 di=rep(0,nG)
 for (g in 1:nG)
  di[g]=sum(dist.f(s1[g,],s2[g,]))
 di=di/(nG-1)

 perm.di=rep(0,nG)
 count.ind=rep(0,nG)
 cat("Starting permutation test:\n")
 for (i in 1:num.permutations){
  cat("permutation",i,"out of",num.permutations,"\n")
  i1=as.vector(sample(1:n,n1))
  perm.X1=X[i1,]
  perm.X2=X[-i1,]
  perm.s1=f(perm.X1,...)
  perm.s2=f(perm.X2,...)
  for (g in 1:nG)
   perm.di[g]=sum(dist.f(perm.s1[g,],perm.s2[g,]))
  perm.di=perm.di/(nG-1)
  for (g in 1:nG)
  if (perm.di[g]>=di[g])
    count.ind[g]=count.ind[g]+1
 }

 p.value.ind=rep(0,nG)
 for (g in 1:nG){
   p.value.ind[g]=count.ind[g]/num.permutations
 }
 names(p.value.ind)=colnames(X1)
 names(di)=colnames(X1)
 new("resultsIndTest",p.values=p.value.ind,d=di)
}
