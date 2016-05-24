#global script for HCPC2
if((inherits(x, "PCA") | inherits(x, "MCA") | inherits(x, "CA")| inherits(x, "MFA"))){
  results=x
  anafact=lignecode
  baba=HCPC(results,nb.clust=-1,graph=FALSE)$call$t$nb.clust
  nbindiv=dim(results$ind$coord)[1]
  nomData=nomData
  clustdf=baba
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

if(!((inherits(x, "PCA") | inherits(x, "MCA") | inherits(x, "CA")| inherits(x, "MFA")))){
if(inherits(x, "HCPCshiny")){
anafact=x$anafact
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
results=x$data
title1=x$title1
title2=x$title2
title3=x$title3
}

if(inherits(x, "PCAshiny") | inherits(x, "CAshiny") | inherits(x, "MCAshiny")){
  results=x$anafact
  anafact=x$code1
  baba=HCPC(results,nb.clust=-1,graph=FALSE)$call$t$nb.clust
  nbindiv=dim(results$ind$coord)[1]
  nomData=nomData
  clustdf=baba
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

baba=HCPC(results,nb.clust=-1,graph=FALSE)$call$t$nb.clust
nbindiv=dim(results$ind$coord)[1]
}