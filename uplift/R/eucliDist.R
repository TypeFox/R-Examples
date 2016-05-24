######################################################################
# Euclidean distance split criteria
######################################################################

  eucliDist <- function(n_data,
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
                         
  ### Euclidean gain
  eucli.node <- (pr.y1_ct1 - pr.y1_ct0) ^ 2 + ((1 - pr.y1_ct1) - (1 - pr.y1_ct0)) ^ 2
  eucli.l <- (pr.y1_l.ct1 - pr.y1_l.ct0) ^ 2 + ((1 - pr.y1_l.ct1) - (1 - pr.y1_l.ct0)) ^ 2
  eucli.r <- (pr.y1_r.ct1 - pr.y1_r.ct0) ^ 2 + ((1 - pr.y1_r.ct1) - (1 - pr.y1_r.ct0)) ^ 2
  eucli.lr <- pr.l * eucli.l + pr.r * eucli.r
  eucli.gain <- eucli.lr - eucli.node
 
  ### Euclidean Normalization factor
  gini.ct <- 2 * pr.ct1 * (1 - pr.ct1) 
  eucli.ct <- (pr.l_ct1 - pr.l_ct0) ^ 2 + ((1 - pr.l_ct1) - (1 - pr.l_ct0)) ^ 2
  gini.ct1 <- 2 * pr.l_ct1 * (1 - pr.l_ct1)
  gini.ct0 <- 2 * pr.l_ct0 * (1 - pr.l_ct0)
  eucli.norm <- gini.ct * eucli.ct + gini.ct1 * pr.ct1  + gini.ct0 * pr.ct0 + 0.5
     
  ### Output
  s.value.t <- eucli.gain / eucli.norm
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