#global script for HCPC for dataframe2
if(is.data.frame(x)==TRUE){
  quanti=c()
  quali=c()
  posi=c()
  for ( i in 1:dim(x)[2]){
    if (is.numeric(x[,i])==TRUE){
      quanti=c(quanti,colnames(x)[i])
    }
    else{
      quali=c(quali,colnames(x)[i])
    }
  }
  VariableChoices=quanti
  nom=rownames(x)
  nbindiv=length(nom)
  num=c(1:length(nom))
  QualiChoice=quali
  IdChoices=c(1:length(VariableChoices))
  Idqualisup=c(1:length(QualiChoice))
  
  nomData=nomData
  consolidf=FALSE
  metricdf="euc"
  drawdf=FALSE
  df=FALSE
  centerdf=FALSE
  numdf=60
  nb1df=1
  nb2df=2
  title1="Hierarchical clustering on the factor map"
  title2="Factor map"
  title3="Hierarchical Clustering"
}

if(is.data.frame(x)==FALSE){

if(inherits(x, "HCPCshiny")){
  nomData=x$nomData
  clustdf=x$clust
  consolidf=x$consoli
  metricdf=x$metric
  drawdf=x$drawtree
  df=x$nom3D
  centerdf=x$center
  numdf=x$num
  nb1df=x$nb1
  nb2df=x$nb2
  x=x$data
  title1=x$title1
  title2=x$title2
  title3=x$title3
}

if(inherits(x, "HCPC")){
  nomData=x$call$call[2]
  clustdf=x$call$t$nb.clust
  consolidf=FALSE
  if(x$call$t$tree[7]=="euclidean"){
    metricdf="euc"
  }
  if(x$call$t$tree[7]=="manhattan"){
    metricdf="manh"
  }
  drawdf=FALSE
  df=FALSE
  centerdf=FALSE
  numdf=60
  nb1df=1
  nb2df=2
  title1="Hierarchical clustering on the factor map"
  title2="Factor map"
  title3="Hierarchical Clustering"
}


quanti=c()
quali=c()
posi=c()
for ( i in 1:dim(x)[2]){
  if (is.numeric(x[,i])==TRUE){
    quanti=c(quanti,colnames(x)[i])
  }
  else{
    quali=c(quali,colnames(x)[i])
  }
}

VariableChoices=quanti
nom=rownames(x)
nbindiv=length(nom)
num=c(1:length(nom))
QualiChoice=quali
IdChoices=c(1:length(VariableChoices))
Idqualisup=c(1:length(QualiChoice))
}