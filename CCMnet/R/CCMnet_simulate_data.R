
NS_Multinomial <- function(g, Network_stats, mean_inflate = 0, var_inflate = 1, covPattern = NULL) {
  if (g$gal$bipartite > 0) {
    print("Not implemented")
  } else {
    NS_Multinomial_uni(g, Network_stats, mean_inflate, var_inflate, covPattern = NULL)    
  }
}


NS_Multinomial_uni <- function(g, Network_stats, mean_inflate = 0, var_inflate = 1, covPattern = NULL) {
  
  error = 0
  
  if (Network_stats == "DegreeDist") {
    mean_NS = tabulate(degree(g, gmode = 'graph')+1) + mean_inflate
    sample_size = sum(mean_NS)
    
    p = mean_NS / sum(mean_NS)
    var_inflate_num = network.size(g)
  } else if (Network_stats == "DegMixing") {
    nedges = c(network.edgecount(g),0,0)
    max_degree_g = max(degree(g, gmode='graph'))
    g_dmm = matrix(0,  nrow = max_degree_g, ncol = max_degree_g)
    edge_list = unlist(g$mel)
    dim(edge_list) = c(3,nedges[1])
    g_degree = degree(g, gmode = "graph")
    for (num_edge in c(1:nedges[1])) {
      deg1 = g_degree[edge_list[1,num_edge]]
      deg2 = g_degree[edge_list[2,num_edge]]
      if ((deg1 <= max_degree_g) && (deg2 <= max_degree_g)) {
        g_dmm[deg1, deg2] = g_dmm[deg1,deg2] + 1
        if (deg1 != deg2) {
          g_dmm[deg2, deg1] = g_dmm[deg2,deg1] + 1
        }
      }
    }
    mean_NS = g_dmm[upper.tri(g_dmm, diag = TRUE)] + mean_inflate
    sample_size = sum(mean_NS)
    
    p = mean_NS / sum(mean_NS)
    var_inflate_num = network.edgecount(g)
  } else {
    print("Error: No such network statistic currently implemented")
    error = 1
    return(NULL)
  }
  
  if (error == 0){
    var_x = matrix(data = 0, nrow = length(p), ncol = length(p))
    
    for (i in c(1:length(p))) {
      for (j in c(1:i)) {
        if (i == j) {
          var_x[i,j] = p[i]*(1-p[i])
        } else {
          var_x[j,i] = var_x[i,j] = -p[i]*p[j]
        }
      }
    }
    diag(var_x) = var_inflate*diag(var_x)
  }
  return(list(mean_NS, (var_x*var_inflate_num) ))
}




