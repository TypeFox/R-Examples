cccvc<-
function(dataset,ry,rind,rmet,covar=NULL,int=FALSE,cl=0.95){
if (int==TRUE)cccvc2(dataset,ry,rind,rmet,covar,cl) else cccvc1(dataset,ry,rind,rmet,covar,cl)}