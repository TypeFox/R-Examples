get_quantiles_PPI <-
function(scaled_node_score,scaled_edge_score,PPI,klist,pop_size){


# The function to get the "node_score_term" and "edge_score_term" of a sub-network denoted by "vector"
# "vector" is a binary vector with length equal to the size of the whole network. 
# An element of value "1" indicates the inclusion of that gene in the selected sub-network.


node_edge<-function(sub){    
    n<-length(sub)    
    node_score<-sum(scaled_node_score[sub])/sqrt(n)
    edges<- PPI[,1] %in% sub & PPI[,2] %in% sub
    m<-sum(edges)
    edge_score<-sum(scaled_edge_score[edges])/sqrt(m)
    return(c(node_score,edge_score))
}

n<-length(scaled_node_score)
all_genes<-names(scaled_node_score)
edge_score_term<-vector(length=length(klist),mode="list")
node_score_term<-vector(length=length(klist),mode="list")
save_sub<-NULL
for(i in 1:length(klist)){
     k<-klist[i]
     for(j in 1:pop_size){         
         sub<-random_network_sampling_PPI(k,PPI,all_genes)
         node_edge_score<-node_edge(sub)
         node_score_term[[i]][j]<-node_edge_score[1]
         edge_score_term[[i]][j]<-node_edge_score[2]
         if(is.na(edge_score_term[[i]][j])) save_sub<-sub
         print(c(i,j))
     }   
     
}


     log_abs_edge_node_ratio <- vector(length=length(klist),mode="list")
     for(i in 1:length(klist)){
          log_abs_edge_node_ratio[[i]] <- log10(abs(edge_score_term[[i]]/node_score_term[[i]]))
     }
     b <- NULL
     for(i in 1:length(klist)){
         a <- summary(log_abs_edge_node_ratio[[i]])
         b <- rbind(b,a)
     }

     ratio <- apply(b,2,mean)[-4]
     lambda <- sort(1/(1+10^ratio))
     names(lambda)<-names(ratio)
     return(list(ratio,lambda))
}

