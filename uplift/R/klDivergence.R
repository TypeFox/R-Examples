######################################################################
# Kullback-Leibler divergence split criteria
######################################################################

  klDivergence <- function(n_data,
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
                         
  ### KL Gain
  kl.node <- pr.y1_ct1 * log(pr.y1_ct1/pr.y1_ct0, 2) +
             (1 - pr.y1_ct1) * log((1 - pr.y1_ct1) / (1 - pr.y1_ct0), 2)
  kl.l <- pr.y1_l.ct1 * log(pr.y1_l.ct1 / pr.y1_l.ct0, 2) +
          (1 - pr.y1_l.ct1) * log((1 - pr.y1_l.ct1) / (1 - pr.y1_l.ct0), 2)
  kl.r <- pr.y1_r.ct1 * log(pr.y1_r.ct1 / pr.y1_r.ct0, 2) +
          (1 - pr.y1_r.ct1) * log((1 - pr.y1_r.ct1) / (1 - pr.y1_r.ct0), 2)
  kl.lr <- pr.l * kl.l + pr.r * kl.r
  kl.gain <- kl.lr - kl.node
    
  ### KL Normalization factor
  ent.ct <- -(pr.ct1 * log(pr.ct1, 2) + pr.ct0 * log(pr.ct0, 2))
  kl.ct <- pr.l_ct1 * log(pr.l_ct1 / pr.l_ct0, 2) +
           (1 - pr.l_ct1) * log ((1 - pr.l_ct1) / (1 - pr.l_ct0), 2)
  ent.ct1 <- -(pr.l_ct1 * log(pr.l_ct1, 2) + (1 - pr.l_ct1) * log((1 - pr.l_ct1), 2))
  ent.ct0 <- -(pr.l_ct0 * log(pr.l_ct0, 2) + (1 - pr.l_ct0) * log((1 - pr.l_ct0), 2))
    
  norm <- kl.ct * ent.ct + ent.ct1 * pr.ct1 + ent.ct0 * pr.ct0 + 0.5
    
  ### Output
  s.value.t <- kl.gain / norm
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

