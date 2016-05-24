findbranch.pare<-function(parent)
{
# finds the nodes who have more than 1 child

len<-length(parent)
frekve<-matrix(0,len,1)

for (i in 1:len){
   if (parent[i]>0) frekve[parent[i]]<-frekve[parent[i]]+1
}

tulos<-matrix(0,len,1)

for (i in 1:len){
     #if (parent[i]==0) tulos[i]<-1
     #else 
     if ((parent[i]!=0) && (frekve[parent[i]]>1)){ #result of a branching
                 tulos[parent[i]]<-1
    }
}

if (sum(tulos)==0) ans<-NULL else ans<-which(tulos==1)

return(ans)
}    
