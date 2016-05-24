findbnodes<-function(lst,modenum=1,num=NULL)
{
# prunes from a level set tree "lst" the modes with "num" 
# smallest excess masses 
# or the modes with smaller excess mass than "exmalim"

if (is.null(num)){
    curmodenum<-moodilkm(lst$parent)$lkm
    num<-curmodenum-modenum
}

len<-length(lst$parent)
child.frekve<-matrix(0,len,1)
for (i in 1:len){
    if (lst$parent[i]>0) 
    child.frekve[lst$parent[i]]<-child.frekve[lst$parent[i]]+1
}

ml<-moodilkm(lst$parent)
mode.list<-ml$modloc
roots.of.modes<-matrix(0,length(mode.list),1)
for (aa in 1:length(mode.list)){
    node<-mode.list[aa]
    while ((lst$parent[node]>0) && (child.frekve[lst$parent[node]]==1)){ 
         node<-lst$parent[node]
    }
    roots.of.modes[aa]<-node
}

em<-excmas(lst)
or<-order(em[roots.of.modes],decreasing=TRUE)
#nodes<-ml$modloc[or[1:modenum]]
nodes<-roots.of.modes[or[1:modenum]]

return(nodes=nodes)
}


