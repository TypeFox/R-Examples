######################################################################
# Find best split
######################################################################

  findBestSplit <- function(n_data, 
                            var_ind, 
                            split_method,
                            n.ct1,
                            n.ct0, 
                            pr.ct1,
                            pr.ct0, 
                            pr.y1_ct1, 
                            pr.y1_ct0, 
                            minbucket_ct0, 
                            minbucket_ct1) {
  
  # temporary correction in case probabilities passed to node equal zero
  critmax <- 1E-11
  pr.y1_ct1 <- ifelse(pr.y1_ct1 == 0, critmax, pr.y1_ct1)
  pr.y1_ct0 <- ifelse(pr.y1_ct0 == 0, critmax, pr.y1_ct0)
  
  # Create matrix with variables required later in the computations 
  if (is.numeric(n_data[, var_ind])) {
    
    n_data_small <- data.frame(n_data[, var_ind], ct = n_data$ct, y = n_data$y)
    s.n_data_small <- n_data_small[sort.list(n_data_small[, 1], na.last = NA, method = "quick"), ]
    s.n_mat <- cbind(.var = s.n_data_small[, 1], 
                      cs.y.ct1 = cumsum(s.n_data_small$y & s.n_data_small$ct), 
                      cs.y.ct0 = cumsum(s.n_data_small$y & !s.n_data_small$ct),
                      ncs.y.ct1 = sum(s.n_data_small$y & s.n_data_small$ct) - cumsum(s.n_data_small$y & s.n_data_small$ct),
                      ncs.y.ct0 = sum(s.n_data_small$y & !s.n_data_small$ct) - cumsum(s.n_data_small$y & !s.n_data_small$ct),
                      cs.ct1 = cumsum(s.n_data_small$ct),
                      cs.ct0 = cumsum(s.n_data_small$ct == 0), 
                      ncs.ct1 = sum(s.n_data_small$ct) - cumsum(s.n_data_small$ct), 
                      ncs.ct0 = sum(s.n_data_small$ct == 0) - cumsum(s.n_data_small$ct == 0))
    
    } else {  # if not numeric
      
    n_data_agg <- cbind(tapply(n_data$y, list(factor(n_data[, var_ind]), n_data$ct), sum),
                        tapply(n_data$ct == 0, factor(n_data[, var_ind]), sum),
                        tapply(n_data$ct == 1, factor(n_data[, var_ind]), sum))
    n_data_agg[is.na(n_data_agg)] <- 0  #I might have levels of factors for which I have records in 
                                        # treatment but not in control. 
    s.n_data <- n_data_agg[order(n_data_agg[, 2] / n_data_agg[, 4] - 
                                 n_data_agg[, 1] / n_data_agg[, 3]), ]
    s.n_mat <- cbind(.var = 1:nrow(s.n_data), 
                      cs.y.ct1 = cumsum(s.n_data[, 2]),
                      cs.y.ct0 = cumsum(s.n_data[, 1]),
                      ncs.y.ct1 = sum(s.n_data[, 2]) - cumsum(s.n_data[, 2]),
                      ncs.y.ct0 = sum(s.n_data[, 1]) - cumsum(s.n_data[, 1]),
                      cs.ct1 = cumsum(s.n_data[, 4]),
                      cs.ct0 = cumsum(s.n_data[, 3]),
                      ncs.ct1 = n.ct1 - cumsum(s.n_data[, 4]), 
                      ncs.ct0 = n.ct0 - cumsum(s.n_data[, 3]))       
 
  }
  
  # Is minbucket constraint satisfied?
  minbucket.ok_temp <- (s.n_mat[, 6] >= minbucket_ct1 & 
                        s.n_mat[, 7] >= minbucket_ct0 &
                        s.n_mat[, 8] >= minbucket_ct1 &
                        s.n_mat[, 9] >= minbucket_ct0)
  dups <- duplicated(s.n_mat[, 1], fromLast = TRUE)
  minbucket.ok <- (minbucket.ok_temp & !dups)
  
  if (any(minbucket.ok)) {
  
    # Randomized assignment implies pr.l_ct1 = pr.l_ct0 for all possible splits
    pr.l_ct1 <- (s.n_mat[, 6] + 1) / (n.ct1 + 2) # added Laplace correction
    pr.l_ct0 <- (s.n_mat[, 7] + 1) / (n.ct0 + 2)
    # Both versions below should work out about the same for numeric
    #pr.l <- (1:nrow(s.n_mat) + 1) / (nrow(s.n_mat) + 2)
    #pr.r <- 1 - pr.l
    pr.l <- pr.l_ct1 * pr.ct1 + pr.l_ct0 * pr.ct0
    pr.r <- 1 - pr.l
  
    # Add Laplace correction to probablities
    pr.y1_l.ct1 <- (s.n_mat[, 2] + 1) / (s.n_mat[, 6] + 2)
    pr.y1_l.ct0 <- (s.n_mat[, 3] + 1) / (s.n_mat[, 7] + 2)
    pr.y1_r.ct1 <- (s.n_mat[, 4] + 1) / (s.n_mat[, 8] + 2)
    pr.y1_r.ct0 <- (s.n_mat[, 5] + 1) / (s.n_mat[, 9] + 2) 
    
    # Number of treatment/control observations at left and right child nodes
    cs.ct1 <- s.n_mat[, 6]
    cs.ct0 <- s.n_mat[, 7]
    ncs.ct1 <- s.n_mat[, 8]
    ncs.ct0 <- s.n_mat[, 9]
    
    # This is an alternative way of computing these. Better pass as an argument to the fun as
    # it is the same for all variables 
    #pr.y1_ct1 <- pr.y1_l.ct1 * pr.l_ct1 + pr.y1_r.ct1 * (1 - pr.l_ct1)
    #pr.y1_ct0 <- pr.y1_l.ct0 * pr.l_ct0 + pr.y1_r.ct0 * (1 - pr.l_ct0)
  
    bs.res <- switch(split_method, 
                     KL = klDivergence(n_data,
                                       s.n_mat,
                                       var_ind,
                                       pr.y1_ct1, 
                                       pr.y1_ct0,
                                       pr.l, 
                                       pr.r,
                                       pr.y1_l.ct1,
                                       pr.y1_l.ct0,
                                       pr.y1_r.ct1,
                                       pr.y1_r.ct0,
                                       pr.ct1,
                                       pr.ct0,
                                       pr.l_ct1,
                                       pr.l_ct0,
                                       minbucket.ok),
                                     
                   
                     Chisq =  chisq(n_data,
                                    s.n_mat,
                                    var_ind,
                                    pr.y1_ct1, 
                                    pr.y1_ct0,
                                    pr.l, 
                                    pr.r,
                                    pr.y1_l.ct1,
                                    pr.y1_l.ct0,
                                    pr.y1_r.ct1,
                                    pr.y1_r.ct0,
                                    pr.ct1,
                                    pr.ct0,
                                    pr.l_ct1,
                                    pr.l_ct0,
                                    minbucket.ok),
          
                     ED = eucliDist(n_data,
                                    s.n_mat,
                                    var_ind,
                                    pr.y1_ct1, 
                                    pr.y1_ct0,
                                    pr.l, 
                                    pr.r,
                                    pr.y1_l.ct1,
                                    pr.y1_l.ct0,
                                    pr.y1_r.ct1,
                                    pr.y1_r.ct0,
                                    pr.ct1,
                                    pr.ct0,
                                    pr.l_ct1,
                                    pr.l_ct0,
                                    minbucket.ok),
                     
                     L1 = L1norm(n_data,
                                 s.n_mat,
                                 var_ind,
                                 pr.y1_ct1, 
                                 pr.y1_ct0,
                                 pr.l, 
                                 pr.r,
                                 pr.y1_l.ct1,
                                 pr.y1_l.ct0,
                                 pr.y1_r.ct1,
                                 pr.y1_r.ct0,
                                 pr.ct1,
                                 pr.ct0,
                                 pr.l_ct1,
                                 pr.l_ct0,
                                 minbucket.ok),
                     
                     Int = intSplit(n_data,
                                    s.n_mat,
                                    var_ind,
                                    pr.y1_ct1, 
                                    pr.y1_ct0,
                                    pr.l, 
                                    pr.r,
                                    pr.y1_l.ct1,
                                    pr.y1_l.ct0,
                                    pr.y1_r.ct1,
                                    pr.y1_r.ct0,
                                    pr.ct1,
                                    pr.ct0,
                                    pr.l_ct1,
                                    pr.l_ct0,
                                    cs.ct1, 
                                    cs.ct0,
                                    ncs.ct1,
                                    ncs.ct0,
                                    minbucket.ok))
                   
    } else {
    
   
    bs.res <- list(s.value = -Inf, 
                   x.value = NA)
      
  }
  
  return(bs.res)
  
}

### END FUN

