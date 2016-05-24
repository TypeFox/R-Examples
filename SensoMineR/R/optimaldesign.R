optimaldesign <- function(nbPanelist,nbProd,nbProdByPanelist=nbProd,seed=NULL){

if (is.null(seed)) seed <- sample(1:1000,1)
set.seed(seed)

coupe = nbProd%/%nbProdByPanelist
if (coupe*nbProdByPanelist!=nbProd) print("Warning: the design in the output is unbalanced")

if ((nbProd%%2)==0) nb.will = nbPanelist%/%(nbProd*coupe)
if ((nbProd%%2)==1) nb.will = nbPanelist%/%(2*nbProd*coupe)
if (nb.will*nbProd*coupe<nbPanelist) nb.will=nb.will+1

if ((nb.will*(nbProd*coupe))!=nbPanelist){
 print(paste("The design could be balanced with",(nb.will+1)*(nbProd*coupe),"panelists"))
 nb.will=nb.will+1
}
plan = NULL

for (k in 1:nb.will){
  zz=WilliamsDesign(nbProd, seed= sample(1:100001,1))
  aux = zz[,1:nbProdByPanelist]
  if (coupe>1) {
    for (i in 2:coupe){
      aux = rbind(aux,zz[,(nbProdByPanelist*(i-1)+1):(nbProdByPanelist*i)])
      for (k in 1:nbProd) aux[aux==k]=nbProd+zz[1,k]
      aux = aux-nbProd
    }
  }

  if (is.null(plan)){ plan =  aux
  } else{ plan = rbind(plan,aux)}
}

plan = plan[1:nbPanelist,]
rownames(plan) = paste("panelist",1:nrow(plan),sep=".")
colnames(plan) = paste("rank",1:ncol(plan),sep=".")

valid.rang = matrix(NA,nbProd,nbProdByPanelist)
for (k in 1:ncol(plan)){ for (i in 1:nbProd) valid.rang[i,k]=sum(plan[,k]==i)}
rownames(valid.rang)=paste("Prod",1:nrow(valid.rang))
colnames(valid.rang)=paste("Rank",1:ncol(valid.rang))

valid.succ = matrix(0,nbProd,nbProd)
for (i in 1:nrow(plan)){
 for (k in 1:(nbProdByPanelist-1)){
  valid.succ[plan[i,k],plan[i,k+1]]=valid.succ[plan[i,k],plan[i,k+1]]+1
 }
}
rownames(valid.succ)=paste("Prod",1:nrow(valid.succ))
colnames(valid.succ)=paste("Prod",1:ncol(valid.succ))

return(list(design=plan,valid.rank=valid.rang,valid.succ=valid.succ))
}