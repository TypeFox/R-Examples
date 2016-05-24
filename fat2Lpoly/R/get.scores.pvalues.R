
# fonction qui prend comme argument le resultat de la fonction scores.covs.multiSNP et retourne un tableau des valeur-p 
# (et scores des tests individuels) des differents tests pour tous les SNPs testes.  

# Jordie Croteau
# 16 aout 2012

# correspond au fichier get_scores_and_p-values_function_v3.R dans le dossier "programmes"

# modifie le 27 aout

# correction apportee le 18 mars pour que ca fonctionne meme dans le cas d'un seul SNP teste.

get.scores.pvalues=function(test,joint.tests){

snp.names.mat=test$snp.names.mat
test=test$scores.covs.all.SNPs

# mettre les indices en ordre croissant au cas ou ca n'aurait pas ete fourni dans cet ordre.
if(!is.null(joint.tests)) joint.tests=lapply(joint.tests,function(x) x[order(x)])

# tableau des scores des tests individuels
scores.indiv.mat=array(NA,c(length(test),dim(test[[1]][[1]])[2]))

# tableaux de valeur-p
p.values.indiv.mat=array(NA,c(length(test),dim(test[[1]][[1]])[2]))
p.value.global.vec=rep(NA,length(test))
if(!is.null(joint.tests)) p.value.joint.mat=array(NA,c(length(test),length(joint.tests)))

for(j in 1:length(test))
 {
  scores=apply(test[[j]]$scores.mat,2,sum)
  Sigma=apply(test[[j]]$cov.mat,2:3,sum,na.rm=TRUE)
  scores.std=scores/sqrt(diag(Sigma))
  scores.indiv.mat[j,]=round(scores.std,3)
  p.values.indiv=2*(1-pnorm(abs(scores.std)))	
  scores.mat=matrix(scores,ncol=1)
  p.value.global=try(1-pchisq(as.numeric(t(scores.mat)%*%chol2inv(chol(Sigma))%*%scores.mat),df=dim(Sigma)[1]),silent=TRUE)
  p.values.indiv.mat[j,]=signif(p.values.indiv,3)
  p.value.global.vec[j]=try(signif(p.value.global,3),silent=TRUE)
  
  if(!is.null(joint.tests))
   {
    for(k in 1:length(joint.tests))
     {
      scores.mat.joint=matrix(scores[joint.tests[[k]]],ncol=1)
      p.value.joint=try(1-pchisq(as.numeric(t(scores.mat.joint)%*%chol2inv(chol(Sigma[joint.tests[[k]],joint.tests[[k]]]))%*%scores.mat.joint),df=length(joint.tests[[k]])),silent=TRUE)
      p.value.joint.mat[j,k]=try(signif(p.value.joint,3),silent=TRUE)
     }
   }
 }		

p.value.global.vec=as.numeric(p.value.global.vec)
res=data.frame(snp.names.mat,p.value.global.vec,scores.indiv.mat,p.values.indiv.mat)
col.names.tmp=c("snp.test","global_p",paste("param",1:ncol(p.values.indiv.mat),"score",sep="_"),paste("param",1:ncol(p.values.indiv.mat),"p",sep="_"))
if(ncol(snp.names.mat)==2) col.names.tmp=c("snp.cond",col.names.tmp)
colnames(res)=col.names.tmp

if(!is.null(joint.tests))
 {
  p.value.joint.mat=matrix(apply(p.value.joint.mat,2,as.numeric),nrow=nrow(p.value.joint.mat),ncol=ncol(p.value.joint.mat))
  res=data.frame(snp.names.mat,p.value.global.vec,p.value.joint.mat,scores.indiv.mat,p.values.indiv.mat)
  col.names.tmp=c("snp.test","global_p",paste("params.joint",unlist(lapply(joint.tests,paste,collapse="-")),"p",sep="_"),paste("param",1:ncol(p.values.indiv.mat),"score",sep="_"),paste("param",1:ncol(p.values.indiv.mat),"p",sep="_"))
  if(ncol(snp.names.mat)==2) col.names.tmp=c("snp.cond",col.names.tmp)
  colnames(res)=col.names.tmp
 }
res
}




