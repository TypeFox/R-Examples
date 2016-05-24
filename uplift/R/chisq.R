######################################################################
# Chisq split criteria
######################################################################

  chisq <- function(n_data,
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
                         
  ### Chi-squared gain
  chisq.node <- ((pr.y1_ct1 - pr.y1_ct0) ^ 2) / pr.y1_ct0  +
                (((1 - pr.y1_ct1) - (1 - pr.y1_ct0)) ^ 2) / (1 - pr.y1_ct0) 
  chisq.l <- ((pr.y1_l.ct1 - pr.y1_l.ct0) ^ 2) / pr.y1_l.ct0  +
             (((1 - pr.y1_l.ct1) - (1 - pr.y1_l.ct0)) ^ 2) / (1 - pr.y1_l.ct0)
  chisq.r <- ((pr.y1_r.ct1 - pr.y1_r.ct0) ^ 2) / pr.y1_r.ct0  +
             (((1 - pr.y1_r.ct1) - (1 - pr.y1_r.ct0)) ^ 2) / (1 - pr.y1_r.ct0)
  chisq.lr <- pr.l * chisq.l + pr.r * chisq.r
  chisq.gain <- chisq.lr - chisq.node
 
  ### Chi-squared Normalization factor
  gini.ct <- 2 * pr.ct1 * (1 - pr.ct1) 
  chisq.ct <- ((pr.l_ct1 - pr.l_ct0) ^ 2) / pr.l_ct0  + 
              (((1 - pr.l_ct1) - (1 - pr.l_ct0)) ^ 2) / (1 - pr.l_ct0)
  gini.ct1 <- 2 * pr.l_ct1 * (1 - pr.l_ct1)
  gini.ct0 <- 2 * pr.l_ct0 * (1 - pr.l_ct0)
  chisq.norm <- gini.ct * chisq.ct + gini.ct1 * pr.ct1  + gini.ct0 * pr.ct0 + 0.5
     
  ### Output
  s.value.t <- chisq.gain / chisq.norm
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
    
  
