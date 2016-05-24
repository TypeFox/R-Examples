prune<-function(et)
{

ini<-initial(et$ssr,et$left,et$right)
et$left<-ini$left
et$right<-ini$right
treeseq<-pruseq(et)

treeseq$support<-et$support
return(treeseq)
}
