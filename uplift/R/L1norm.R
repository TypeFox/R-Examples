######################################################################
# Euclidean distance split criteria
######################################################################

     L1norm <- function(n_data,
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
                        minbucket.ok) {
                         
  ### L1 gain
  L1.node <- abs(pr.y1_ct1 - pr.y1_ct0) + abs((1 - pr.y1_ct1) - (1 - pr.y1_ct0)) 
  L1.l <- abs(pr.y1_l.ct1 - pr.y1_l.ct0) + abs((1 - pr.y1_l.ct1) - (1 - pr.y1_l.ct0)) 
  L1.r <- abs(pr.y1_r.ct1 - pr.y1_r.ct0) + abs((1 - pr.y1_r.ct1) - (1 - pr.y1_r.ct0)) 
  L1.lr <- pr.l * L1.l + pr.r * L1.r
  L1.gain <- L1.lr - L1.node
 
  ### L1 Normalization factor
  gini.ct <- 2 * pr.ct1 * (1 - pr.ct1) 
  L1.ct <- abs(pr.l_ct1 - pr.l_ct0) + abs((1 - pr.l_ct1) - (1 - pr.l_ct0)) 
  gini.ct1 <- 2 * pr.l_ct1 * (1 - pr.l_ct1)
  gini.ct0 <- 2 * pr.l_ct0 * (1 - pr.l_ct0)
  L1.norm <- gini.ct * L1.ct + gini.ct1 * pr.ct1  + gini.ct0 * pr.ct0 + 0.5
     
  ### Output
  s.value.t <- L1.gain / L1.norm
  s.value <- max(s.value.t[minbucket.ok])
  wh.max <- which(s.value.t == s.value)
  wh.max <- wh.max[minbucket.ok[wh.max]] #Ensures the max selected satisfies the constraint (in case of duplicates)
  
  ### break ties randomly
  
  if (length(wh.max) > 1) {
    wh.max <- sample(wh.max, 1)
  }
    
  if (is.numeric(n_data[, var_ind])) {
    x.value = s.n_mat[wh.max, 1] 
  } else x.value =  s.n_mat[, 1] <= wh.max
      
  criteria.res <- list(s.value = s.value, 
                       x.value = x.value)
  return(criteria.res)
}
    
  
### END FUN