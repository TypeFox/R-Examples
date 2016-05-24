#' Blocked and Clustered Random Assignment.
#' 
#' Random assignment where units are assigned as clusters and clusters are nested within blocks. 
#'
#' @param clust_var A vector of length N that indicates which cluster each unit belongs to.
#' @param block_var A vector of length N that indicates which block each unit belongs to.
#' @param num_arms The total number of treatment arms. If unspecified, will be determined from the number of columns of block_m, the length of prob_each, or the length of condition_names.
#' @param block_m A matrix of arm sizes with blocks in the rows and treatment conditions in the columns. The rows should respect the alphabetical ordering of the blocks as determined by sort(unique(block_var). The columns should be in the order of condition_names, if specified.
#' @param prob_each A vector whose length is equal to the number of treatment assignments. When specified, prob_each assigns the same (within rounding) proportion of each block to each treatment condition, using complete random assignment. prob_each must sum to 1.
#' @param condition_names A character vector giving the names of the treatment conditions. If unspecified, the treatment conditions. will be named T1, T2, T3, etc.
#'
#' @return A vector of length N that indicates the treatment condition of each unit.
#' 
#' @examples 
#' clust_var <- rep(letters, times=1:26)
#' block_var <- rep(NA, length(clust_var))
#' block_var[clust_var %in% letters[1:5]] <- "block_1"
#' block_var[clust_var %in% letters[6:10]] <- "block_2"
#' block_var[clust_var %in% letters[11:15]] <- "block_3"
#' block_var[clust_var %in% letters[16:20]] <- "block_4"
#' block_var[clust_var %in% letters[21:26]] <- "block_5"
#' 
#' table(block_var, clust_var)
#' 
#' Z <- block_and_cluster_ra(clust_var = clust_var, block_var = block_var)
#' 
#' table(Z, clust_var)
#' table(Z, block_var)
#' 
#' Z <- block_and_cluster_ra(clust_var = clust_var, block_var = block_var, num_arms = 3)
#' 
#' table(Z, clust_var)
#' table(Z, block_var)
#' 
#' Z <- block_and_cluster_ra(clust_var = clust_var, block_var = block_var, prob_each = c(.2, .5, .3))
#' 
#' block_m <- rbind(c(2, 3),
#'                  c(1, 4),
#'                  c(3, 2),
#'                  c(2, 3),
#'                  c(5, 1))
#' 
#' Z <- block_and_cluster_ra(clust_var = clust_var, block_var = block_var, block_m = block_m)
#' 
#' table(Z, clust_var)
#' table(Z, block_var)
#' @export
block_and_cluster_ra <- 
  function(clust_var, block_var, num_arms = NULL, block_m=NULL, prob_each=NULL, condition_names = NULL) {
    
    # confirm that all units within clusters are in the same block
    # is there a computationally faster way to confirm this (possible c++ loop?)
    
    if(!all(rowSums(table(clust_var, block_var) != 0)==1)){
      stop("All units within a cluster must be in the same block.")
    }
    
    # Setup: obtain unique clusters
    unique_clust <- unique(clust_var)
    
    # get the block for each cluster
    clust_blocks <- rep(NA, length(unique_clust))
    for(i in 1:length(unique_clust)){
      clust_blocks[i] <- unique(block_var[clust_var==unique_clust[i]])  
    }
    
    # Conduct random assignment at cluster level
    z_clust <- block_ra(block_var = clust_blocks, 
                        num_arms = num_arms,
                        block_m = block_m, 
                        prob_each = prob_each,
                        condition_names = condition_names)
    
    # Merge back up to the individual level, maintaining original ordering
    merged <- merge(x = data.frame(clust_var, init_order = 1:length(clust_var)), 
                    y = data.frame(clust_var=unique_clust, z_clust), by="clust_var")
    merged <- merged[order(merged$init_order),]
    return(merged$z_clust)
  }
