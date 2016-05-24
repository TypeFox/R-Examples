scaspa<-function(treeseq,bind,eind)
{

alkm<-eind-bind+1
modelkm<-matrix(0,alkm,1)

j<-1
for (i in bind:eind){
     leafnum<-treeseq$leafs[i]
     tree<-eval.pick(treeseq,leafnum)
     pv<-partition(tree)
     values<-pv$values  
     recs<-pv$recs
     if (length(values==1)) modelkm[j]<-1
     else{
       pg<-profgene(values,recs)
       parents<-pg$parent
       mlkm<-moodilkm(parents)
       modelkm[j]<-mlkm$lkm
     }
     j<-j+1
}

leafnums<-treeseq$leafs[bind:eind]
alfas<-treeseq$alfa[bind:eind]
return(list(moodilkm=t(modelkm),alfas=alfas,leafnums=leafnums))
}
