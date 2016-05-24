score_scaling <-
function(diff_expr,diff_coex){

      gene_num<-length(diff_expr)
      n<-gene_num
      unique_edge<-rep(0,choose(n,2)) 
      k=1   
      for(i in 1:(n-1)){
          for(j in (i+1):n){
          unique_edge[k]<-diff_coex[i,j]
          k<-k+1
          }
      }


      hist(diff_expr,breaks=100,main="Diff_expr_before_scaling")
      hist(unique_edge,breaks=100,main="Diff_coex_before_scaling")


      # Standardize the node score and edge score

      mean_node<-mean(diff_expr)
      sigma_node<-sd(diff_expr)
      diff_expr<-(diff_expr-mean_node)/sigma_node

      mean_edge<-mean(unique_edge)
      sigma_edge<-sd(unique_edge)
      diff_coex<-(diff_coex-mean_edge)/sigma_edge   

      hist(diff_expr,breaks=100,main="Diff_expr_after_scaling")
      hist((unique_edge-mean_edge)/sigma_edge,breaks=100,main="Diff_coex_after_scaling")

      return(list(diff_expr,diff_coex))
}

