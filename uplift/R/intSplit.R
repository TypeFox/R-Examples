######################################################################
# Interaction split criteria
# Based on Radcliffe and Surry "Real-world uplift modeling with 
# significance based uplift trees" Technical Report (Eq. 17)
######################################################################

intSplit <- function(n_data,
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
                     minbucket.ok) {
  
  
  ### Compute elements for split formula
  C44 <- 1/cs.ct1 + 1/cs.ct0 + 1/ncs.ct1 + 1/ncs.ct0
  
  UR <- pr.y1_r.ct1 - pr.y1_r.ct0
  UL <- pr.y1_l.ct1 - pr.y1_l.ct0
  
  SSE <- cs.ct1 * pr.y1_l.ct1 * (1 - pr.y1_l.ct1)  + 
         ncs.ct1 * pr.y1_r.ct1 * (1 - pr.y1_r.ct1) +
         cs.ct0 * pr.y1_l.ct0 * (1 - pr.y1_l.ct0)  +
         ncs.ct0 * pr.y1_r.ct0 * (1 - pr.y1_r.ct0)
         
  n.node <- nrow(n_data)       
  
  ### Interaction split
  s.value.t <- ((n.node - 4) * (UR - UL)^2) / (C44 * SSE)
  
  
  ### Output
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