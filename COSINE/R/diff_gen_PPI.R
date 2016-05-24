diff_gen_PPI <-
function(data1,data2,PPI){
    num_sample_1 <- dim(data1)[1]
    num_sample_2 <- dim(data2)[1]
    num_gene <- dim(data1)[2]
    PPI <-PPI[PPI[,1]%in% colnames(data1) & PPI[,2]%in% colnames(data1),]
    num_edge <- dim(PPI)[1]
    node_score <- rep(0,num_gene)
    edge_score <- rep(0,num_edge)
    type <- c(rep(0,num_sample_1),rep(1,num_sample_2))

    # calculate the statistics measuing the differential expression of each gene between the 2 groups

    for(i in 1:num_gene){
       data <- c(data1[,i],data2[,i])
       node_score[i] <- f.test(data,type)
       print(i)
    }

    # calculate the statistics measuing the differential co-expression of each gene-pair between the 2 groups

    for(i in 1:num_edge){
              gene1 <- as.character(PPI[i,1])
              gene2 <- as.character(PPI[i,2])
              data.x <- c(data1[,gene1],data2[,gene1])
              data.y <- c(data1[,gene2],data2[,gene2])
              edge_score[i] <- cov(data.x,data.y)
              print(i)
    }
    
    scaled_node_score <- (node_score - mean(node_score))/sd(node_score)
    scaled_edge_score <- (edge_score - mean(edge_score))/sd(edge_score)
    names(scaled_node_score) <- colnames(data1)
    return(list(scaled_node_score, scaled_edge_score))
}

