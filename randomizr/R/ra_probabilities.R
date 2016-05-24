#' Probabilties of assignment: Simple Random Assignment
#'
#' @param N The total number of units in the experimental sample (required).
#' @param prob The probability of assignment to treatment. If specified, a two-group design is assumed.
#' @param num_arms The total number of treatment arms. If unspecified, num_arms will be determined from the length of m_each or condition_names.
#' @param prob_each A numeric giving the probability of assignment to each treatment arm. Must sum to 1. Please note that due to rounding, these probabilities are approximate. For finer control, please use m_each.
#' @param condition_names A character vector giving the names of the treatment groups. If unspecified, the treatment groups will be named T1, T2, T3, etc. An execption is a two-group design in which N only or N and m are specified, in which the condition names are 0 and 1.
#'
#' @return A matrix of probabilities of assignment
#' 
#' @examples 
#' # Two Group Designs
#' simple_ra_probabilities(N=100)
#' simple_ra_probabilities(N=100, prob=0.5)
#' simple_ra_probabilities(N=100, prob_each = c(0.3, 0.7), 
#'                         condition_names = c("control", "treatment"))
#' # Multi-arm Designs
#' simple_ra_probabilities(N=100, num_arms=3)
#' simple_ra_probabilities(N=100, prob_each=c(0.3, 0.3, 0.4))
#' simple_ra_probabilities(N=100, prob_each=c(0.3, 0.3, 0.4), 
#'                         condition_names=c("control", "placebo", "treatment"))
#' simple_ra_probabilities(N=100, condition_names=c("control", "placebo", "treatment"))
#' 
#' @export
simple_ra_probabilities <- function(N, prob = NULL, num_arms = NULL, prob_each = NULL, condition_names = NULL){
  
  # Setup: obtain number of arms and condition_names
  
  if(is.null(num_arms)){
    num_arms <- 2
    if(!is.null(prob_each)){num_arms <- length(prob_each)}
    if(!is.null(condition_names)){num_arms <- length(condition_names)}
  }
  
  if(is.null(condition_names)){
    if(num_arms==2){
      condition_names = c(0,1)
    }else{
      condition_names <- paste0("T", 1:num_arms)    
    }
  }
  if(is.null(prob) & is.null(prob_each)){
    prob_each <- rep(1/num_arms, num_arms)
  }
  if(!is.null(prob)){
    condition_probabilities <- c(1-prob, prob)
  }
  if(!is.null(prob_each)){
    condition_probabilities <- prob_each
  }
  
  # Build prob_mat
  prob_mat <- matrix(rep(condition_probabilities, N), 
                     byrow=TRUE, ncol=length(condition_probabilities), 
                     dimnames = list(NULL,  paste0("prob_",condition_names)))
  return(prob_mat)
  
}

#' Probabilties of assignment: Complete Random Assignment
#'
#' @param N The total number of units in the experimental sample (required).
#' @param m If specified, a two-group design is assumed. m is the total number of units to be assigned to treatment. Should only be specified for a two group design in which exactly m of N units are assigned to treatment. If not specified, half of the sample (N/2) will be assigned to treatment (if N is odd, m will be set to either floor(N/2) or ceiling(N/2) with equal probability. m is NULL by default. 
#' @param prob The probability of assignment to treatment. If specified, a two-group design is assumed.
#' @param num_arms The total number of treatment arms. If unspecified, num_arms will be determined from the length of m_each, prob_each, or condition_names.
#' @param m_each A numeric vector giving the size of each treatment group. Must sum to N. If unspecified, equally sized (rounded) groups will be assumed.
#' @param prob_each A numeric giving the probability of assignment to each treatment arm. Must sum to 1. Please note that due to rounding, these probabilities are approximate. For finer control, please use m_each.
#' @param condition_names A character vector giving the names of the treatment groups. If unspecified, the treatment groups will be named T1, T2, T3, etc. An execption is a two-group design in which N only or N and m are specified, in which the condition names are 0 and 1.
#'
#' @return A matrix of probabilities of assignment
#' 
#' @examples 
#' # 2-arm designs
#' complete_ra_probabilities(N=100)
#' complete_ra_probabilities(N=100, m=50)
#' complete_ra_probabilities(N=100, prob = .3)
#' 
#' complete_ra_probabilities(N=100, m_each = c(30, 70), 
#'                           condition_names = c("control", "treatment"))
#' 
#' # Multi-arm Designs
#' complete_ra_probabilities(N=100, num_arms=3)
#' complete_ra_probabilities(N=100, m_each=c(30, 30, 40))
#' 
#' complete_ra_probabilities(N=100, m_each=c(30, 30, 40), 
#'                           condition_names=c("control", "placebo", "treatment"))
#' 
#' complete_ra_probabilities(N=100, condition_names=c("control", "placebo", "treatment"))
#' complete_ra_probabilities(N=100, prob_each = c(.2, .7, .1))
#' 
#' @export
complete_ra_probabilities <- function(N, m = NULL, prob = NULL, num_arms = NULL, m_each = NULL, prob_each = NULL, condition_names = NULL){
  
  # Setup: obtain number of arms and condition_names
  
  if(is.null(num_arms)){
    num_arms <- 2
    if(!is.null(m_each)){num_arms <- length(m_each)}
    if(!is.null(prob_each)){num_arms <- length(prob_each)}
    if(!is.null(condition_names)){num_arms <- length(condition_names)}
  }
  
  if(is.null(condition_names)){
    if(num_arms==2){
      condition_names = c(0,1)
    }else{
      condition_names <- paste0("T", 1:num_arms)    
      }
  }
  
  # Case 0: Two Arms and N = 1
  if(is.null(m_each) & is.null(prob_each) & num_arms ==2 & N ==1) {
    prob <- 0.5
    prob_mat <- matrix(rep(c(1-prob, prob), N), byrow=TRUE, ncol=2, dimnames = list(NULL,  paste0("prob_",condition_names)))
    return(prob_mat)
  }
  
  # Case 1: Two Arms and N > 1
  if(is.null(m_each) & is.null(prob_each) & num_arms==2 & N > 1){
    m_floor <- m
    m_ceiling <- m
    
    if(is.null(m)){
      m_floor <- floor(N/2)
      m_ceiling <- ceiling(N/2)
    }
    
    if(!is.null(prob)){
      m_floor <- floor(N*prob)
      m_ceiling <- ceiling(N*prob)
    }
    
    prob <- 0.5*(m_floor/N) + 0.5*(m_ceiling/N)
    prob_mat <- matrix(rep(c(1-prob, prob), N), byrow=TRUE, ncol=2, 
                       dimnames = list(NULL,  paste0("prob_",condition_names)))
    return(prob_mat)
  }
  
  # Case 2: We need to obtain "condition_probabilities" then make a matrix.
  
  # 2a: If m_each is specified
  if(!is.null(m_each) & is.null(prob_each)){
    remainder <-  N%%num_arms
    condition_probabilities <- (m_each/N)
  }
  
  # 2b: if neither m_each nor prob_each is specified
  if(is.null(m_each) & is.null(prob_each)){
    m_each <- rep(N%/%num_arms, num_arms)
    remainder <-  N%%num_arms
    condition_probabilities <- 
      (1-(remainder/num_arms))*(m_each/N) +
      (remainder/num_arms)*((m_each +1)/N)
  }
  
  # 2c: if prob_each is specified
  if(!is.null(prob_each)){
    m_each <- floor(N*prob_each)
    remainder <- N - sum(m_each)
    condition_probabilities <- 
      (1-(remainder/length(prob_each)))* (m_each/N) +
      (remainder/length(prob_each))* ((m_each +1)/N)
  } 
  
  # 2d: if N is smaller than number of arms, we just flip coins
  if(N < num_arms){
    condition_probabilities <- rep(N/num_arms, num_arms)
  }
  
  # Build prob_mat
  prob_mat <- matrix(rep(condition_probabilities, N), 
                     byrow=TRUE, ncol=length(condition_probabilities), 
                     dimnames = list(NULL,  paste0("prob_",condition_names)))
  return(prob_mat)
  
}

#' Probabilties of assignment: Block Random Assignment
#'
#' @param block_var A vector of length N indicating which block each unit belongs to.
#' @param num_arms The total number of treatment arms. If unspecified, will be determined from the number of columns of block_m, the length of prob_each, or the length of condition_names.
#' @param block_m A matrix of arm sizes whose number of rows is equal to the number of blocks and whose number of columns is equal to the number of treatment arms. The rows should respect the alphabetical ordering of the blocks as determined by sort(unique(block_var). The columns should be in the order of condition_names, if specified.
#' @param prob_each A numeric vector whose length is equal to the number of treatment conditions. When specified, prob_each assigns the same (within rounding) proportion of each block to each treatment condition, using complete random assignment. prob_each must sum to 1.
#' @param condition_names A character vector giving the names of the treatment conditions. If unspecified, the treatment conditions. will be named T1, T2, T3, etc.
#'
#' @return A matrix of probabilities of assignment
#' 
#' @examples 
#' 
#' block_var <- rep(c("A", "B","C"), times=c(50, 100, 200))
#' block_ra_probabilities(block_var=block_var)
#' 
#' block_m <- rbind(c(25, 25),
#'                  c(50, 50),
#'                  c(100, 100))
#' 
#' block_ra_probabilities(block_var=block_var, block_m=block_m)
#' 
#' block_m <- rbind(c(10, 40),
#'                  c(30, 70),
#'                  c(50, 150))
#' 
#' block_ra_probabilities(block_var=block_var, block_m=block_m, 
#'                        condition_names=c("control", "treatment"))
#' 
#' block_ra_probabilities(block_var=block_var, num_arms=3)
#' 
#' block_m <- rbind(c(10, 20, 20),
#'                  c(30, 50, 20),
#'                  c(50, 75, 75))
#' block_ra_probabilities(block_var = block_var, block_m = block_m)
#' 
#' block_ra_probabilities(block_var=block_var, block_m=block_m, 
#'                        condition_names=c("control", "placebo", "treatment"))
#' 
#' block_ra_probabilities(block_var=block_var, prob_each=c(.1, .1, .8))
#' 
#' 
#' @export
block_ra_probabilities <- function(block_var, num_arms = NULL, block_m=NULL, prob_each = NULL, condition_names = NULL){
  
  # Setup: obtain number of arms and condition_names
  
  if(is.null(num_arms)){
    num_arms <- 2
    if(!is.null(block_m)){num_arms <- dim(block_m)[2]}
    if(!is.null(prob_each)){num_arms <- length(prob_each)}
    if(!is.null(condition_names)){num_arms <- length(condition_names)}
  }
  
  if(is.null(condition_names)){
    if(num_arms==2){
      condition_names = c(0,1)
    }else{
      condition_names <- paste0("T", 1:num_arms)    
    }
  }
  
  blocks <- sort(unique(block_var))
  prob_mat <- matrix(NA, 
                     nrow = length(block_var), 
                     ncol = length(condition_names),
                     dimnames = list(NULL,  paste0("prob_",condition_names)))
  
  # Case 1: Assume (approximately) equal probabilities for all blocks and conditions.
  if(is.null(block_m) & is.null(prob_each) & is.null(prob_each)){
    for(i in 1:length(blocks)){
      N_block <- sum(block_var==blocks[i])
      prob_mat[block_var==blocks[i],] <- complete_ra_probabilities(N = N_block, condition_names=condition_names)
    }
    return(prob_mat)
  }
  
  # Case 2: block_m is specified
  if(!is.null(block_m)){
    for(i in 1:length(blocks)){
      N_block <- sum(block_var==blocks[i])
      prob_mat[block_var==blocks[i],] <- complete_ra_probabilities(N = N_block, 
                                                                   m_each = block_m[i,], 
                                                                   condition_names=condition_names)
    }
    return(prob_mat)
  }
  
  # Case 3: prob_each is specified
  if(!is.null(prob_each)){
    for(i in 1:length(blocks)){
      N_block <- sum(block_var==blocks[i])
      prob_mat[block_var==blocks[i],] <- complete_ra_probabilities(N = N_block, 
                                                                   prob_each = prob_each, 
                                                                   condition_names=condition_names)
    }
    return(prob_mat)
  }
  
  # Case 4: prob_each is specified
  if(!is.null(prob_each)){
    for(i in 1:length(blocks)){
      N_block <- sum(block_var==blocks[i])
      prob_mat[block_var==blocks[i],] <- complete_ra_probabilities(N = N_block, 
                                                                   prob_each = prob_each[i,], 
                                                                   condition_names=condition_names)
    }
    return(prob_mat)
  }
}

#' Probabilties of assignment: Cluster Random Assignment
#'
#' @param clust_var A vector of length N that indicates which cluster each unit belongs to.
#' @param m The total number clusters to be treated. Should only be specified for a two group design in which exactly m of N clusters is assigned to treatment. If not specified, half of the clusters will be assigned to treatment. Is NULL by default. 
#' @param num_arms The total number of treatment arms. If unspecified, will be determined from the length of m_each or condition_names.
#' @param m_each A numeric vector giving the number of clusters to be assigned to each treatment group. Must sum to the total number of clusters. If unspecified, equally sized (rounded) groups will be assumed.
#' @param prob_each A numeric vector giving the probability of assignment to each treatment arm. Must sum to 1. Please note that due to rounding, these probabilities are approximate. For finer control, please use m_each.
#' @param condition_names A character vector giving the names of the treatment groups.  If unspecified, the treatment groups will be named T1, T2, T3, etc. 
#'
#' @return A matrix of probabilities of assignment
#' 
#' @examples 
#' 
#' # Two Group Designs
#' clust_var <- rep(letters, times=1:26)
#' cluster_ra_probabilities(clust_var=clust_var)
#' 
#' cluster_ra_probabilities(clust_var=clust_var, m=10)
#' 
#' cluster_ra_probabilities(clust_var=clust_var, m_each = c(9, 17), 
#'                          condition_names = c("control", "treatment"))
#' 
#' # Multi-arm Designs
#' cluster_ra_probabilities(clust_var=clust_var, num_arms=3)
#' cluster_ra_probabilities(clust_var=clust_var, m_each=c(7, 7, 12))
#' 
#' cluster_ra_probabilities(clust_var=clust_var, m_each=c(7, 7, 12), 
#'                          condition_names=c("control", "placebo", "treatment"))
#' 
#' cluster_ra_probabilities(clust_var=clust_var, 
#'                          condition_names=c("control", "placebo", "treatment"))
#' 
#' cluster_ra_probabilities(clust_var=clust_var, prob_each = c(.1, .2, .7))
#' 
#' 
#' 
#' @export
cluster_ra_probabilities <- function(clust_var, m=NULL, num_arms = NULL, m_each = NULL, prob_each = NULL, condition_names = NULL){
  unique_clus <- unique(clust_var)
  n_clus <- length(unique_clus)
  probs_clus <- complete_ra_probabilities(N = n_clus, m = m, num_arms = num_arms, m_each = m_each, prob_each = prob_each, condition_names = condition_names)
  merged <- merge(x = data.frame(clust_var, init_order = 1:length(clust_var)), 
                  data.frame(clust_var=unique_clus, probs_clus), by="clust_var")
  merged <- merged[order(merged$init_order),]
  prob_mat <- as.matrix(merged[,colnames(probs_clus)])
  return(prob_mat)
}

#' Probabilties of assignment: Blocked and Clustered Random Assignment
#'
#' @param clust_var A vector of length N that indicates which cluster each unit belongs to.
#' @param block_var A vector of length N that indicates which block each unit belongs to.
#' @param num_arms The total number of treatment arms. If unspecified, will be determined from the number of columns of block_m or the length of condition_names. 
#' @param block_m A matrix of arm sizes whose number of rows is equal to the number of blocks and whose number of columns is equal to the number of treatment arms. The rows should respect the alphabetical ordering of the blocks as determined by sort(unique(block_var). The columns should be in the order of condition_names, if specified.
#' @param prob_each A vector whose length is equal to the number of treatment assignments. When specified, prob_each assigns the same (within rounding) proportion of each block to each treatment condition, using complete random assignment. prob_each must sum to 1.
#' @param condition_names A character vector giving the names of the treatment conditions. If unspecified, the treatment conditions. will be named T1, T2, T3, etc.
#'
#' @return A matrix of probabilities of assignment
#' 
#' @examples 
#' 
#' clust_var <- rep(letters, times=1:26)
#' block_var <- rep(NA, length(clust_var))
#' block_var[clust_var %in% letters[1:5]] <- "block_1"
#' block_var[clust_var %in% letters[6:10]] <- "block_2"
#' block_var[clust_var %in% letters[11:15]] <- "block_3"
#' block_var[clust_var %in% letters[16:20]] <- "block_4"
#' block_var[clust_var %in% letters[21:26]] <- "block_5"
#' 
#' 
#' block_and_cluster_ra_probabilities(clust_var = clust_var, 
#'                                    block_var = block_var)
#' block_and_cluster_ra_probabilities(clust_var = clust_var, 
#'                                    block_var = block_var, 
#'                                    num_arms = 3)
#' block_and_cluster_ra_probabilities(clust_var = clust_var, 
#'                                    block_var = block_var, 
#'                                    prob_each = c(.2, .5, .3))
#' 
#' block_m <- rbind(c(2, 3),
#'                  c(1, 4),
#'                  c(3, 2),
#'                  c(2, 3),
#'                  c(5, 1))
#' 
#' block_and_cluster_ra_probabilities(clust_var = clust_var, block_var = block_var, block_m = block_m)
#' 
#' 
#' @export
block_and_cluster_ra_probabilities <- 
  function(clust_var, block_var, num_arms = NULL, block_m=NULL, prob_each=NULL, condition_names = NULL){
    unique_clus <- unique(clust_var)
    
    ## get the block for each cluster
    clust_blocks <- rep(NA, length(unique_clus))
    for(i in 1:length(unique_clus)){
      clust_blocks[i] <- unique(block_var[clust_var==unique_clus[i]])  
    }
    
    probs_clus <- block_ra_probabilities(block_var = clust_blocks, 
                                         block_m = block_m,
                                         num_arms = num_arms,
                                         prob_each = prob_each,
                                         condition_names = condition_names)
    
    merged <- merge(x = data.frame(clust_var, init_order = 1:length(clust_var)), 
                    data.frame(clust_var=unique_clus, probs_clus), by="clust_var")
    merged <- merged[order(merged$init_order),]
    prob_mat <- as.matrix(merged[,colnames(probs_clus)])
    return(prob_mat)
  }



