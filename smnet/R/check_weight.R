
check_weight <- function(adj, wgt, silent = F){
  wgt <- try(as.numeric(as.matrix(wgt)))
  if(class(wgt) == "numeric" & (!anyNA(wgt))){
    # pull out the bid from the adjacency object
    bid  <- adj$rid_bid[,2]
    nbid <- nchar(bid) 
    min_weight <- min(wgt)
    # the idea here is that for a network weight, the sum of the weights at the i^th 
    # level of the network heirarchy should equal the number of segments in the i-1^th level
    # of the network, minus the number of 'source' segments (dead_ends) that occur at i-1^th level
    # for an additive weight, the idea is similar, the sum of the weights a the i^th level
    # equal the sum of the weights at the i-1^th level, minus the number of dead_ends at i-1
    # since all dead_ends will have order 1 (under Shreve - if the lowest order is different, I'll need to fix this)
    dis_additive     <- 0
    dis_network      <- 0
    for(i in 2:max(nbid)){
      which_lower     <- which(nbid == i-1)
      which_upper     <- which(nbid == i)
      n_lower         <- length(which_lower)
      sum_up_weight   <- sum(wgt[which_upper])
      sum_dn_weight   <- sum(wgt[which_lower])
      dead_ends       <- sum(colSums(adj$adjacency)[which_lower] == 0)
      dis_network     <- dis_network + (n_lower - sum_up_weight - dead_ends)^2
      dis_additive    <- dis_additive + (sum_dn_weight - sum_up_weight - dead_ends*min_weight)^2
    }
    if(round(dis_additive, 6) == 0){
      weight.type <- "additive"
      if(!silent) cat("Provided weight passes additivity check... \n")
    } else if(round(dis_network, 6)  == 0){
      weight.type <- "network"
      if(!silent) cat("Provided weight passes network check... \n")
    } else {weight.type <- "unrecognised"}
  } else {
    weight.type <- "unrecognised"
  }
  weight.type
}


