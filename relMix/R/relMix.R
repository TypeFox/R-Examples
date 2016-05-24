relMix <-
function(pedigrees,locus,datamatrix,kinship=0){
if(class(locus)!= "FamiliasLocus")
  stop("first argument locus should be a Familias locus")
if(!is.data.frame(datamatrix))
  stop("second argument datamatrix should be a dataframe")
for (i in 1:length(pedigrees)){
  if(class(pedigrees[[i]])!="pedigree") 
  stop("Third argument pedigrees should be a list of Familias pedigree-s")
  }
n=dim(datamatrix)[2]/2
myloci=vector("list",n)
for (i in 1:n) {
  myloci[[i]]=locus
  myloci[[i]]$locusname=paste("Term",i,sep="")
}
alt=FamiliasPosterior(pedigrees, myloci, datamatrix,ref=2,kinship=kinship)
result=apply(alt$likelihoodsPerSystem,2,sum)
terms=alt$likelihoodsPerSystem
terms=terms[terms[,1]>0|terms[,2]>0,]
list(likIsFather=result[1],likUnrelated=result[2],
  LR=result[1]/result[2],terms=terms)
}
