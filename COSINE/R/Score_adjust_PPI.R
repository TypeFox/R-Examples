Score_adjust_PPI <-
function(scaled_node_score,scaled_edge_score,PPI,lam,subnet,num_random_sampling,best_score)
{


    all_genes<-names(scaled_node_score)
    node_num<-length(subnet)
    genes_selected<-all_genes[subnet]
    edges_selected<- PPI[,1] %in% genes_selected & PPI[,2] %in% genes_selected
    num_edges_selected<-sum(edges_selected)

    #Random sampling
   
    random_score<-rep(0,num_random_sampling) 
    
    for(i in 1:num_random_sampling){

       sampled_edges <- sample(1:dim(PPI)[1],num_edges_selected)
       edge_score<-sum(scaled_edge_score[sampled_edges])/sqrt(num_edges_selected)
       sampled_nodes <- sample(1:length(all_genes),node_num)
       node_score<-sum(scaled_node_score[sampled_nodes])/sqrt(node_num)
       random_score[i]<- lam*edge_score + (1-lam)*node_score
       print(i)
    }

    mean<-mean(random_score)
    sd<-sd(random_score)
    adjusted_score<-(best_score-mean)/sd  
    return (adjusted_score)
}

