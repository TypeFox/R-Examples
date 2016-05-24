get_components_PPI <-
function(gene_names,vector,PPI,minsize){

  selected_node<-gene_names[vector==1]

  # Remove nodes with degree 0

  related_PPI <- PPI[PPI[,1] %in% selected_node & PPI[,2] %in% selected_node,]
  connected_node <- union(related_PPI[,1],related_PPI[,2])
  component_list <- vector(mode="list",length=1)
  t<-0

  while(length(connected_node)>0){
    seed <- connected_node[1]
    current_component <- seed
    genes_to_add<-NULL
    previous_length<-0
    current_length<-1
    flag <- 1
    iter <- 0
    while(flag==1){
      iter <- iter+1
      for(i in (previous_length + 1):current_length ){
    #     print(i)
         current_gene<-current_component[i]
         a1<-which(related_PPI[,1]==current_gene)
         b1<-which(related_PPI[,2]==current_gene)
         genes<-union(related_PPI[a1,2],related_PPI[b1,1])    
         genes_to_add <- union(genes,genes_to_add)
         genes_to_add <-setdiff(genes_to_add,current_component)
     }
     if(length(genes_to_add)==0) flag<-0  else{
         previous_length<-current_length
         current_component<-union(current_component,genes_to_add)
         current_length<-length(current_component)
         }
   #  print(paste("Finished iteration:",iter))
    } 

    if(length(current_component)>=minsize){
       t <- t+1
       component_list[[t]]<-current_component
   #    print(paste("Finished component",t))
    }
    connected_node<-setdiff(connected_node,current_component)
  }

  return (component_list)
}

