pickseq<-function(treeseq,suppo,lsets=FALSE,invalue=FALSE,parvec=NULL)
{

if (is.null(parvec)) leaf<-treeseq$leafs
else leaf<-treeseq$leafs[parvec]
alfa<-treeseq$alfa[parvec]

alkm<-length(parvec)
for (inds in alkm:1){  # start with the oversmoothed estimate 
     leafnum<-leaf[inds]
   
     tree<-eval.pick(treeseq,leafnum)
     #pv<-partition(tree,suppo)
     #curtree<-profgene(pv$values,pv$recs,frekv=F,cvol=T,ccen=T,cfre=F)
     curtree<-proftree(tree)

     if (inds==alkm){
        if (alkm==1){
            treelist<-curtree
        }
        else{
           treelist=list(curtree)
        }
     }
     else{
        treelist=c(treelist,list(curtree))
     }
}

return(list(treelist=treelist,alfa=alfa,leaf=leaf))
}
