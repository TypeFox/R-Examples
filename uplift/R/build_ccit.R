######################################################################
# Build causal conditional inference tree 
######################################################################

  build_ccit <- function(b.dframe, 
                         mtry, 
                         split_method,
                         interaction.depth,
                         pvalue,
                         bonferroni,
                         minsplit, 
                         minbucket_ct0, 
                         minbucket_ct1, 
                         nr_vars,
                         var_names,
                         nr_in_samples,
                         nr_nodes,
                         ...) {
  
  # create variables used to store results at each iteration
  s_curr_node <- rep(NA, nr_nodes)
  s_node_status <- rep(NA, nr_nodes)
  s_n.ct1 <- rep(NA, nr_nodes)
  s_n.ct0 <- rep(NA, nr_nodes)
  s_pr.y1_ct1 <- rep(NA, nr_nodes)
  s_pr.y1_ct0 <- rep(NA, nr_nodes)
  s_bs.var <- rep(NA, nr_nodes)
  s_bs.x.value <- vector("list", nr_nodes)
  s_bs.s.value <- rep(NA, nr_nodes)
                                    
  #initialize iterations and observation index
  iter <- 1 
  curr_node <- 1
  obs_node <- rep(1, nr_in_samples) #initilize observation index (which observations belong to which node at the current iteration)

  ### Start building

  repeat{
  ### check split requirements are met

    repeat {
  
    node_status <- 2 #indicates what type a node is: 1: parent node, 2: has not been processed yet, -1: terminal node
    obs_curr_node.ind <- which(obs_node == curr_node) # which observations will be split at the current iteration

    ### compute misc. arguments passed to findBestSplit
    n.ct1 <- sum(b.dframe[obs_curr_node.ind, ]$ct)
    n.ct0 <- sum(b.dframe[obs_curr_node.ind, ]$ct == 0)  
    pr.ct1 <- (n.ct1 + 1) / (n.ct1 + n.ct0 + 2)
    pr.ct0 <- 1 - pr.ct1
    pr.y1_ct1 <- sum(b.dframe[obs_curr_node.ind, ]$y & b.dframe[obs_curr_node.ind, ]$ct) / n.ct1
    pr.y1_ct0 <- sum(b.dframe[obs_curr_node.ind, ]$y & !b.dframe[obs_curr_node.ind, ]$ct) / n.ct0 
  
   # check which variables have less than 2 levels at the current node, 
   # and exclude them as candidates for mtry selection
    single.lev <- logical(nr_vars)
    for (i in 1:nr_vars) {
      single.lev[i] <- length(unique(b.dframe[obs_curr_node.ind, i])) == 1
    }

    ok_vars <- (1:nr_vars)[!single.lev]

    if (length(ok_vars) >= mtry) {
      mtry.ind <- sample(ok_vars, mtry) #select mtry variables
    } 
    
    ### termination conditions
    last_node <- curr_node == max(obs_node)
    split_cond <- all(pr.y1_ct1 > 0 & pr.y1_ct0 > 0 & pr.y1_ct1 < 1 & pr.y1_ct0 < 1 #pure node?
                      & n.ct1 > minsplit & n.ct0 > minsplit # minsplit condtion satisfied?
                      & length(ok_vars) >= mtry # enough ok variables in node?
                      & (is.null(interaction.depth) ||  curr_node < 2 ^ interaction.depth))

    
    ### Independence test
    if (split_cond) {
      
      indTest <- vapply(mtry.ind, 
                        function(i) 
                          pvalue(independence_test(as.formula(paste('z ~', 
                                                                    paste(var_names[i]))),
                                                   data = b.dframe[obs_curr_node.ind, ], ...)), FUN.VALUE = 0);
      
      if (bonferroni) {
        indTest_pass <- any(indTest <= pvalue / length(mtry.ind)) } else {
          indTest_pass <- any(indTest <= pvalue)
        }
      
      if (!indTest_pass) split_cond <- FALSE
      
    }
    
    
    if (!split_cond) node_status <- -1 #update node status to terminal
  
    ### store node stats
    s_curr_node[iter] <- curr_node
    s_node_status[iter] <- node_status
    s_n.ct1[iter] <- n.ct1
    s_n.ct0[iter] <- n.ct0
    s_pr.y1_ct1[iter] <- pr.y1_ct1
    s_pr.y1_ct0[iter] <- pr.y1_ct0

    # go to next node
    if (!last_node)  {
      curr_node_temp <- curr_node
      curr_node <- min(unique(obs_node)[curr_node < unique(obs_node)])
      iter <- iter + 1
    }

    if (last_node | split_cond)  break 

}

  if (!last_node) {
  curr_node <- curr_node_temp
  iter <- iter - 1
  }

  if (split_cond) { # note that from the loop above, I can have last_node = T and split_cond = F 
    
    
    bs.var <- mtry.ind[which(indTest ==  min(indTest))] ### best var
    
    ### break ties randomly
    if (length(bs.var) > 1) {
      bs.var <- sample(bs.var, 1)
    }
     
    bs.mat <- findBestSplit(n_data = b.dframe[obs_curr_node.ind, ], 
                            var_ind = bs.var, 
                            split_method,
                            n.ct1,
                            n.ct0, 
                            pr.ct1,
                            pr.ct0, 
                            pr.y1_ct1, 
                            pr.y1_ct0, 
                            minbucket_ct0, 
                            minbucket_ct1);
    
    bs.s.value <- bs.mat$s.value
    bs.x.value <- bs.mat$x.value
    num.var <- is.numeric(b.dframe[, bs.var])
    
    if (bs.s.value > 0) {
      
      node_status = 1
      
      ### update obs index based on best partition
      if (num.var) {
        obs_node[obs_curr_node.ind] <- ifelse(b.dframe[obs_curr_node.ind, bs.var] <= bs.x.value, 2 * curr_node,
                                              2 * curr_node + 1)
      } else {
        obs_node[obs_curr_node.ind] <- ifelse(b.dframe[obs_curr_node.ind, bs.var] %in%
                                                names(bs.x.value[bs.x.value == TRUE]), 2 * curr_node,
                                              2 * curr_node + 1)
      }      
      
    } else {
      node_status <- -1
      bs.var <- bs.x.value <- bs.s.value <- NA
    } 
    
    ### store split stats 
 
    s_node_status[iter] <- node_status
    s_bs.var[iter] <- bs.var
    s_bs.x.value[[iter]] <- bs.x.value
    s_bs.s.value[iter] <- bs.s.value

    ### go to next iteration
    iter <- iter + 1
    ### give me next node for which to attempt a partition
    if (max(obs_node) > curr_node)  {curr_node <- min(unique(obs_node)[curr_node < unique(obs_node)])} else break;
  } else break; # if last node = T and split conditions are not satisfied
}
  total_nr_nodes <- sum(!is.na(s_curr_node))
  res.tree <- list(total_nr_nodes = total_nr_nodes,
                    s_curr_node = s_curr_node[1:total_nr_nodes],
                    s_node_status = s_node_status[1:total_nr_nodes],
                    s_n.ct1 = s_n.ct1[1:total_nr_nodes],
                    s_n.ct0 = s_n.ct0[1:total_nr_nodes],
                    s_pr.y1_ct1 = round(s_pr.y1_ct1[1:total_nr_nodes], 4),
                    s_pr.y1_ct0 = round(s_pr.y1_ct0[1:total_nr_nodes], 4),
                    s_bs.var = s_bs.var[1:total_nr_nodes],
                    s_bs.x.value = s_bs.x.value[1:total_nr_nodes][],
                    s_bs.s.value = s_bs.s.value[1:total_nr_nodes])
                   
}
### END FUN
  
  
  

                      









