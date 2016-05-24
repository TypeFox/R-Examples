profdelt<-function(treeseq,leafnum,suppo,frekv=NULL,cvol=FALSE,ccen=FALSE,cfre=FALSE){
#Profiles an adaptive histogram

tree<-NULL
#tree<-picktree(treeseq,leafnum)

# Tehdaan binaaripuusta paloittain vakio
pv<-partition(tree,suppo)
recs<-pv$recs
values<-pv$values
#                               
pg<-profgene(values,recs,frekv,cvol,ccen,cfre)
#
return(pg)
}


