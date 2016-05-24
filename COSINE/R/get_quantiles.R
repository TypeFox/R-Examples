get_quantiles <-
function(diff_expr,diff_coex,klist,pop_size){


##### The function to get the "node_score_term" and "edge_score_term" of a sub-network denoted by "vector"

my.fun<-function(vector){
return(diff_coex[vector[1],vector[2]])
}


# "vector" is a binary vector with length equal to the size of the whole network. 
# An element of value "1" indicates the inclusion of that gene in the selected sub-network.


node_edge<-function(vector){    
    selected_subset<-which(vector==1)
    n<-length(selected_subset)    
    node_score<-sum(diff_expr[selected_subset])/sqrt(n)
    edges<-combn(selected_subset,2)
    edge_score<-sum(apply(edges,2,my.fun))/sqrt(choose(n,2))
    return(c(node_score,edge_score))
}

n<-length(diff_expr)

edge_score_term<-vector(length=length(klist),mode="list")

node_score_term<-vector(length=length(klist),mode="list")

for(i in 1:length(klist)){
     k<-klist[i]
     for(j in 1:pop_size){
         sub<-sample(n,k)
         vector<-rep(0,n)
         vector[sub]<-1
         node_edge_score<-node_edge(vector)
         node_score_term[[i]][j]<-node_edge_score[1]
         edge_score_term[[i]][j]<-node_edge_score[2]
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

