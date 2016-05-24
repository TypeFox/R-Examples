findneighbor<-function(lst,node)
{
mu<-multitree(lst$parent)

no<-lst$parent[node]
while ((no!=0) && (mu$sibling[mu$child[no]]==0)){
     no<-lst$parent[no]
}


return(no)
}

