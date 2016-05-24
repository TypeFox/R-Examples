test.modular.structure=function(X1,X2,scores="PLS",min.module.size=5,epsilon=.5,num.permutations=1000,check.networks=TRUE,...){
 X1=as.matrix(X1)
 X2=as.matrix(X2)
 if (check.networks==TRUE){
  PairNetworks=new("pairOfNetworks",network1=X1,network2=X2)
  CommonNetworks=get.common.networks(PairNetworks)
  X1=CommonNetworks$network1
  X2=CommonNetworks$network2
 }
 if (is.function(scores)==TRUE)
  R.test.modular.structure(X1,X2,scores,min.module.size,epsilon,num.permutations,...)
 else if (scores=="PLS")
       PLSnet.test.modular.structure(X1,X2,min.module.size,epsilon,num.permutations,...)
 else if (scores=="PC")
       PCnet.test.modular.structure(X1,X2,min.module.size,epsilon,num.permutations,...)
 else if (scores=="cor")
       cornet.test.modular.structure(X1,X2,min.module.size,epsilon,num.permutations,...)
 else if (scores=="RR")
       RRnet.test.modular.structure(X1,X2,min.module.size,epsilon,num.permutations,...)
 else
  stop("Error: scores invalid!")
}

PLSnet.test.modular.structure=function(X1,X2,min.module.size=5,epsilon=.5,num.permutations=1000,ncom=3,rescale.data=TRUE,symmetrize.scores=TRUE,rescale.scores=FALSE){
 n1=as.integer(nrow(X1))
 n2=as.integer(nrow(X2))
 p=as.integer(ncol(X1))
 gene.names=colnames(X1) 
 if (is.null(gene.names))
  gene.names=paste("Gene",1:p)
 out=.C("tdmsPLS",as.double(X1),as.double(X2),as.integer(min.module.size),as.double(epsilon),pval=double(1),sN=double(1),module1=integer(p),module2=integer(p),n1,n2,p,as.integer(ncom),as.integer(num.permutations),as.integer(rescale.data),as.integer(symmetrize.scores),as.integer(rescale.scores),PACKAGE="dna")
 names(out$module1)=gene.names
 names(out$module2)=gene.names
 new("resultsModTest",p.value=out$pval,N=out$sN,modules1=new("modules",module=as.factor(out$module1)),modules2=new("modules",module=as.factor(out$module2)))
}

PCnet.test.modular.structure=function(X1,X2,min.module.size=5,epsilon=.5,num.permutations=1000,ncom=3,rescale.data=TRUE,symmetrize.scores=TRUE,rescale.scores=FALSE){
 n1=as.integer(nrow(X1))
 n2=as.integer(nrow(X2))
 p=as.integer(ncol(X1))
 gene.names=colnames(X1) 
 if (is.null(gene.names))
  gene.names=paste("Gene",1:p)
 out=.C("tdmsPC",as.double(X1),as.double(X2),as.integer(min.module.size),as.double(epsilon),pval=double(1),sN=double(1),module1=integer(p),module2=integer(p),n1,n2,p,as.integer(ncom),as.integer(num.permutations),as.integer(rescale.data),as.integer(symmetrize.scores),as.integer(rescale.scores),PACKAGE="dna")
 names(out$module1)=gene.names
 names(out$module2)=gene.names
 new("resultsModTest",p.value=out$pval,N=out$sN,modules1=new("modules",module=as.factor(out$module1)),modules2=new("modules",module=as.factor(out$module2)))
}

RRnet.test.modular.structure=function(X1,X2,min.module.size=5,epsilon=.5,num.permutations=1000,lambda=1,rescale.data=TRUE,symmetrize.scores=TRUE,rescale.scores=FALSE){
 n1=as.integer(nrow(X1))
 n2=as.integer(nrow(X2))
 p=as.integer(ncol(X1))
 gene.names=colnames(X1) 
 if (is.null(gene.names))
  gene.names=paste("Gene",1:p)
 out=.C("tdmsRR",as.double(X1),as.double(X2),as.integer(min.module.size),as.double(epsilon),pval=double(1),sN=double(1),module1=integer(p),module2=integer(p),n1,n2,p,as.double(lambda),as.integer(num.permutations),as.integer(rescale.data),as.integer(symmetrize.scores),as.integer(rescale.scores),PACKAGE="dna")
 names(out$module1)=gene.names
 names(out$module2)=gene.names
 new("resultsModTest",p.value=out$pval,N=out$sN,modules1=new("modules",module=as.factor(out$module1)),modules2=new("modules",module=as.factor(out$module2)))
}

cornet.test.modular.structure=function(X1,X2,min.module.size=5,epsilon=.5,num.permutations=1000,rescale.scores=FALSE){
s1=cornet(X1,rescale.scores)
s2=cornet(X2,rescale.scores)
n1=nrow(X1)
n2=nrow(X2)
X=rbind(X1,X2)
p=ncol(X)
n=n1+n2

modF1=network.modules(abs(s1),min.module.size,epsilon)
modF2=network.modules(abs(s2),min.module.size,epsilon)
F1=get.modules(modF1)
F2=get.modules(modF2)
G0=(F1!=0)|(F2!=0)
sNc=0
for (g in which(G0))
  if ((F1[g]!=0)&(F2[g]!=0))
  sNc=sNc+(F1[g]!=0)*(F2[g]!=0)*length(intersect(which(F1==F1[g]),which(F2==F2[g])))/length(union(which(F1==F1[g]),which(F2==F2[g])))
if (sum(G0)==0)
  sN=0
else
  sN=1-sNc/sum(G0)
names(sN)=NULL
permutation.list=sample(1:n,n1)
for (i in 2:num.permutations)
 permutation.list=rbind(permutation.list,sample(1:n,n1))
num.perm=nrow(permutation.list)
perm.sN=rep(0,num.perm)
cat("Starting permutation test:\n")
cat("permutation 1 out of",num.permutations,"\n")

i=1
perm.sG0=rep(0,num.perm)

perm.mF1=rep(0,num.perm)
perm.mF2=rep(0,num.perm)
perm.sF1=rep(0,num.perm)
perm.sF2=rep(0,num.perm)
while (i<=num.perm){
  cat("permutation",i,"out of",num.permutations,"\n")
  i1=as.vector(permutation.list[i,])
  perm.X1=X[i1,]
  perm.X2=X[-i1,]
  perm.s1=cornet(perm.X1,rescale.scores)
  perm.s2=cornet(perm.X2,rescale.scores)
  perm.F1=get.modules(network.modules(abs(perm.s1),min.module.size,epsilon))
  perm.F2=get.modules(network.modules(abs(perm.s2),min.module.size,epsilon))
  perm.G0=(perm.F1!=0)|(perm.F2!=0)
  sNc=0
  for (g in which(perm.G0))
  if ((perm.F1[g]!=0)&(perm.F2[g]!=0))
    sNc=sNc+(perm.F1[g]!=0)*(perm.F2[g]!=0)*length(intersect(which(perm.F1==perm.F1[g]),which(perm.F2==perm.F2[g])))/length(union(which(perm.F1==perm.F1[g]),which(perm.F2==perm.F2[g])))
  if (sum(perm.G0)==0)
  perm.sN[i]=0
  else
  perm.sN[i]=1-sNc/sum(perm.G0)
  i=i+1
}
p.value=mean(perm.sN>=sN)
new("resultsModTest",p.value=p.value,N=sN,modules1=modF1,modules2=modF2)
}

R.test.modular.structure=function(X1,X2,f,min.module.size=5,epsilon=.5,num.permutations=1000,...){
s1=f(X1,...)
s2=f(X2,...)
n1=nrow(X1)
n2=nrow(X2)
X=rbind(X1,X2)
p=ncol(X)
n=n1+n2

modF1=network.modules(abs(s1),min.module.size,epsilon)
modF2=network.modules(abs(s2),min.module.size,epsilon)
F1=get.modules(modF1)
F2=get.modules(modF2)
G0=(F1!=0)|(F2!=0)
sNc=0
for (g in which(G0))
  if ((F1[g]!=0)&(F2[g]!=0))
  sNc=sNc+(F1[g]!=0)*(F2[g]!=0)*length(intersect(which(F1==F1[g]),which(F2==F2[g])))/length(union(which(F1==F1[g]),which(F2==F2[g])))
if (sum(G0)==0)
  sN=0
else
  sN=1-sNc/sum(G0)
names(sN)=NULL
permutation.list=sample(1:n,n1)
for (i in 2:num.permutations)
 permutation.list=rbind(permutation.list,sample(1:n,n1))
num.perm=nrow(permutation.list)
perm.sN=rep(0,num.perm)
cat("Starting permutation test:\n")
cat("permutation 1 out of",num.permutations,"\n")

i=1
perm.sG0=rep(0,num.perm)

perm.mF1=rep(0,num.perm)
perm.mF2=rep(0,num.perm)
perm.sF1=rep(0,num.perm)
perm.sF2=rep(0,num.perm)
while (i<=num.perm){
  cat("permutation",i,"out of",num.permutations,"\n")
  i1=as.vector(permutation.list[i,])
  perm.X1=X[i1,]
  perm.X2=X[-i1,]
  perm.s1=f(perm.X1,...)
  perm.s2=f(perm.X2,...)
  perm.F1=get.modules(network.modules(abs(perm.s1),min.module.size,epsilon))
  perm.F2=get.modules(network.modules(abs(perm.s2),min.module.size,epsilon))
  perm.G0=(perm.F1!=0)|(perm.F2!=0)
  sNc=0
  for (g in which(perm.G0))
  if ((perm.F1[g]!=0)&(perm.F2[g]!=0))
    sNc=sNc+(perm.F1[g]!=0)*(perm.F2[g]!=0)*length(intersect(which(perm.F1==perm.F1[g]),which(perm.F2==perm.F2[g])))/length(union(which(perm.F1==perm.F1[g]),which(perm.F2==perm.F2[g])))
  if (sum(perm.G0)==0)
  perm.sN[i]=0
  else
  perm.sN[i]=1-sNc/sum(perm.G0)
  i=i+1
}
p.value=mean(perm.sN>=sN)
new("resultsModTest",p.value=p.value,N=sN,modules1=modF1,modules2=modF2)
}
