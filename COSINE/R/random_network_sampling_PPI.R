random_network_sampling_PPI <-
function(size,PPI,all_genes){
   a<-0
   while(sum(a)==0){     
    current<-sample(all_genes,size)   
    a<- (PPI[,1] %in% current & PPI[,2] %in% current)
   }
   return(current)
}

