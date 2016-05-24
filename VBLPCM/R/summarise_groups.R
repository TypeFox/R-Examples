vblpcmgroups<-function(v.params, colours=palette())
  {
  V_lambda<-v.params$V_lambda
  col_vec<-apply(V_lambda,2,which.max)
  cat("The ", v.params$G, " Groups are:\n", sep="")
  
  for (g in 1:v.params$G)
    {
    cat("\n*****************", colours[g], "**********************\n")
    for (i in 1:sum(col_vec==g))
      cat(network.vertex.names(v.params$net)[col_vec==g][i], "\n") 
    }
  }
