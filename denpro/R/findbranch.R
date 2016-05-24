findbranch<-function(parent,colo="red1",pch=22)
{
# finds the nodes which make the tree of the branches

#pch=19: solid circle, pch=20: bullet (smaller circle), 
#pch=21: circle, pch=22: square, 
#pch=23: diamond, pch=24: triangle point-up, 
#pch=25: triangle point down. 

len<-length(parent)
frekve<-matrix(0,len,1)

for (i in 1:len){
   if (parent[i]>0) frekve[parent[i]]<-frekve[parent[i]]+1
}

tulos<-matrix(0,len,1)
colovec<-matrix("black",len,1)
pchvec<-matrix(21,len,1)

for (i in 1:len){
    if (parent[i]==0){ #root node
             tulos[i]<-1  
             colovec[i]<-colo
             pchvec[i]<-pch
 
    }  
    else if (frekve[parent[i]]>1){ #result of a branching
                 tulos[i]<-1
                 colovec[i]<-colo
                 pchvec[i]<-pch  
    }
}

return(list(indicator=tulos,colovec=colovec,pchvec=pchvec))
}    
