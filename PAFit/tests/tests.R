#testing codes
result <- NULL
data   <- NULL
stats  <- NULL
library(PAFit)
for (prob_m in c("TRUE", "FALSE"))
   for (inc in c("TRUE","FALSE"))
      for (log in c("TRUE", "FALSE"))             
              for (i in 1:3) {
                  data  <- GenerateNet(N = 50, m = 5,prob_m = prob_m, increase = inc, log = log,
                                      mode = i, shape = 1, rate = 1)
                  for (bin in c("TRUE","FALSE"))   
                    for (deg_thresh in c(0,2)) {  
                        stats <- GetStatistics(data$graph,deg_threshold = deg_thresh, Binning = bin, G = 10) 
                        #check stats
                        if (sum(stats$m_t) != sum(stats$Sum_m_k))
                            print("wrong at m_t and sum_m_k")
                        if (sum(abs(colSums(stats$m_tk) - stats$Sum_m_k)) != 0)
                            print("wrong at m_tk and sum_m_k")
                        temp <- sapply(1:(stats$T-1),function(x) sum(stats$node_degree[x,] != -1))
                        if (sum(abs(rowSums(stats$n_tk) - (rowSums(stats$offset_tk) + temp))))
                            print("wrong at node_degree, n_tk, offset_tk") 
                        if (sum(stats$z_j) > sum(stats$m_t))
                            print("wrong at z_j")   
                        result <- PAFit(stats, mode_f ="Constant_PA",weight_PA_mode = 0,stop_cond = 10^-1)
                        plot(result,stats,plot = "A")
                        plot(result,stats,plot = "f")
                        plot(result,data = stats,true_f = data$fitness,plot = "true_f")
                  }
    }
