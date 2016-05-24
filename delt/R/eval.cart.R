eval.cart<-function(dendat,leaf,minobs=NULL)
{
n<-dim(dendat)[1]
if (is.null(minobs)) minobs<-ceiling(sqrt(n)/2)

bt<-densplit(dendat,minobs)
treeseq<-prune(bt)
aplnum<-roundlnum(treeseq$leafs,leaf)
eval<-eval.pick(treeseq,aplnum)  
eval$lnum<-aplnum

return(eval)
}



