#Returns a membership matrix from a cluster
#NA values will be unassigned objects (0 for all)
as.memb<-function(cluster){
   cln =levels(as.factor(cluster))
   k = length(cln)
   u = matrix(0,nrow=length(cluster),ncol=k)    
   for(i in 1:k) {
   	u[,i]=ifelse(cluster==cln[i],1,0)
   } 
   u[is.na(u)]<-0
   if(!is.null(names(cluster))) rownames(u) = names(cluster)
   colnames(u) = cln
   return(u)
}