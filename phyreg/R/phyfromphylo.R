phyfromphylo <-
function(phylo) {     

     Nnode<-phylo$Nnode
     edge<-phylo$edge
     Ntip<-nrow(edge)-Nnode+1
     edge.length<-phylo$edge.length
     
     worklist<-c(Ntip+1)
     oldnames<-c()
     newnames<-c()
     newname<-Ntip+Nnode
     while (length(worklist)>0L) {
        for (kk in 1:length(worklist)) {oldnames<-union(oldnames,worklist[kk]) ; newnames<-union(newnames,newname); newname<-newname-1}
        nextworklist<-c();
        for (kk in 1:length(worklist)) for (jj in 1:nrow(edge)) if (edge[jj,1] %in% worklist) if (edge[jj,2]>Ntip+0.5) nextworklist<-union(nextworklist,edge[jj,2])
        worklist<-nextworklist
       }
       
       # Now we have the renaming vectors, we create a subscripting vector whose ii'th element is the new name of old node ii
       
       trans<-rep(0,Ntip+Nnode)
       trans[1:Ntip]<-1:Ntip;
       trans[oldnames]<-newnames
       
       # Then we just take trans of each element of edge
       
       edge<-trans[edge]
       edge<-array(edge,c(Ntip+Nnode-1,2))
       
       # and we have translated so each node's parent has a higher number than it does. Then we order by the daughter node (2nd column)
       
       edge<-edge[ order(edge[,2]) , ]
       
       # and then the first column is what we want
       
       outphy<-c(edge[,1])
       
       if (!all(edge[,2]==seq(1,Ntip+Nnode-1))) {cat(paste(" edge not as it should be as col2 is not 1:length ")); browser()}
       
       # Then we need also to convert edge.length if it exists-- using the old untranslated edge as stored in phylo$edge
       
       if (!is.null(edge.length)) {
       edge.length[trans[phylo$edge[,2]]]<-phylo$edge.length
       
       # and then turn into heights  - for convenience first with the root as zero and each node as the distance to the root
       
       hts<-rep(0,Ntip+Nnode)
       tempdata<-c(edge.length, 0)
       for ( ii in length(outphy):1) hts[ii]<-tempdata[ii]+hts[outphy[ii]] 

       # Now turn them round
       
       hts<-max(hts)-hts
       
       }
       
       if (is.null(edge.length)) return(list(phy=outphy)) else  return(list(phy=outphy , hts=hts))
         
}
