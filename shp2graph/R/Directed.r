# edgelist: result from reading a shp file, the list of edges
# direction.v: a binary vector to specify if an edge is directed or bidirected, 1 means directed, 0 otherwise
Directed<-function(edgelist, direction.v=rep(0,length(edgelist[,1])), eadf=NULL)
{
   newElist<-c()
   ne<-length(edgelist[,1])
   newEID<-0
   neweadf<-c()
   if(!is.null(eadf))
   {
      anames<-names(eadf)
   }
   for (i in 1:ne)
   {
     if (direction.v[i]==0)
     {
       newEID=newEID+1
       edge<-c(newEID, edgelist[i,2],edgelist[i,3])
       newElist<-rbind(newElist, edge) 
       neweadfi<-eadf[i,]
       neweadf<-rbind(neweadf, neweadfi)
       newEID=newEID+1
       edge<-c(newEID, edgelist[i,3],edgelist[i,2])
       newElist<-rbind(newElist, edge)
       neweadf<-rbind(neweadf, neweadfi)
     }
     else
     {
       newEID=newEID+1
       edge<-c(newEID, edgelist[i,2],edgelist[i,3])
       newElist<-rbind(newElist, edge)
       neweadfi<-eadf[i,]
       neweadf<-rbind(neweadf, neweadfi)  
     }
   }
   if(!is.null(eadf))
   {
      neweadf<-data.frame(neweadf) 
      neweadf<-cbind(newElist[,1],neweadf)
      names(neweadf)<-c("NEdgeID", anames)
   }
   rownames(neweadf)<-NULL
   rownames(newElist)<-NULL
   res<-list(newElist, neweadf)
   res 
}
