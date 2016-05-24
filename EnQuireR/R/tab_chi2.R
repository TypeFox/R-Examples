#Tableau avec les résultats des test du chi-2
"tab_chi2"=function(don,X,Y,seuil,print=TRUE){
#don:jeu de données
#X et Y : numéros des 2 groupes de variables à tester
#seuil:pour colorier les cellules
mat_pc=mat_chi=matrix(0,length(X),length(Y))
b=1
options(warn=-1)
for (i in Y){
a=1
  for (j in X){
  mat_pc[a,b]= chisq.test(don[,i],don[,j],correct=FALSE)[[3]]
  mat_chi[a,b]= chisq.test(don[,i],don[,j],correct=FALSE)[[1]]
  #options(warn=-1)
  a=a+1
  }
b=b+1
}

rownames(mat_chi)=rownames(mat_pc)=colnames(don)[X]
colnames(mat_chi)=colnames(mat_pc)=colnames(don)[Y]

mat_chi=mat_chi

if(print==TRUE){
coltable(mat_chi,col.mat=mat_pc,main.title="",level.lower = seuil)
tpolice=par("cex")
title("Chi-2 test",cex=tpolice)
}
options(warn=1)
return(list(mat_chi,mat_pc))
}
