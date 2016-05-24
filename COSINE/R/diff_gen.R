diff_gen <-
function(data1,data2){
    num_sample_1 <- dim(data1)[1]
    num_sample_2 <- dim(data2)[1]
    num_gene <- dim(data1)[2]
    diff_expr <- rep(0,num_gene)
    diff_coex <- matrix(0,ncol=num_gene,nrow=num_gene)
    type <- c(rep(0,num_sample_1),rep(1,num_sample_2))

    # calculate the statistics measuing the differential expression of each gene between the 2 groups

    for(i in 1:num_gene){
       data <- c(data1[,i],data2[,i])
       diff_expr[i] <- f.test(data,type)
    }

    # calculate the statistics measuing the differential co-expression of each gene-pair between the 2 groups

    for(i in 1:(num_gene-1)){
         for(j in (i+1):num_gene){
              data.x <- c(data1[,i],data2[,i])
              data.y <- c(data1[,j],data2[,j])
              cond_fyx <- cond.fyx(data.y,data.x,type)
              cond_fxy <- cond.fyx(data.x,data.y,type)
              diff_coex[i,j] <- (cond_fyx+cond_fxy)/2
              diff_coex[j,i] <- diff_coex[i,j]
        }
    }

    return(list(diff_expr,diff_coex))
}

