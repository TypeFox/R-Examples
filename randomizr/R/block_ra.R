#' Block Random Assignment
#'
#' Random assignment where complete random assignment is employed within experimental blocks. For example, imagine that 50 of 100 men are assigned to treatment and 75 of 200 women are assigned to treatment.
#' 
#' @param block_var A vector of length N that indicates which block each unit belongs to.
#' @param num_arms The total number of treatment arms. If unspecified, will be determined from the number of columns of block_m, the length of prob_each, or the length of condition_names.
#' @param block_m A matrix of arm sizes with blocks in the rows and treatment conditions in the columns. The rows should respect the alphabetical ordering of the blocks as determined by sort(unique(block_var)). The columns should be in the order of condition_names, if specified.
#' @param prob_each A numeric vector whose length is equal to the number of treatment conditions. When specified, prob_each assigns the same (within rounding) proportion of each block to each treatment condition, using complete random assignment. prob_each must sum to 1.
#' @param condition_names A character vector giving the names of the treatment conditions. If unspecified, the treatment conditions will be named T1, T2, T3, etc.
#' @return A vector of length N that indicates the treatment condition of each unit.
#' @export
#' @examples
#' block_var <- rep(c("A", "B","C"), times=c(50, 100, 200))
#' Z <- block_ra(block_var=block_var)
#' table(block_var, Z)
#' 
#' block_m <- rbind(c(25, 25),
#'                  c(50, 50),
#'                  c(100, 100))
#' 
#' Z <- block_ra(block_var=block_var, block_m=block_m)
#' table(block_var, Z)
#' 
#' block_m <- rbind(c(10, 40),
#'                  c(30, 70),
#'                  c(50, 150))
#' 
#' Z <- block_ra(block_var=block_var, block_m=block_m, 
#'               condition_names=c("control", "treatment"))
#' table(block_var, Z)
#' 
#' # Multi-arm Designs
#' Z <- block_ra(block_var=block_var, num_arms=3)
#' table(block_var, Z)
#' 
#' block_m <- rbind(c(10, 20, 20),
#'                  c(30, 50, 20),
#'                  c(50, 75, 75))
#' Z <- block_ra(block_var=block_var, block_m=block_m )
#' table(block_var, Z)
#' 
#' Z <- block_ra(block_var=block_var, block_m=block_m, 
#'               condition_names=c("control", "placebo", "treatment"))
#' table(block_var, Z)
#' 
#' Z <- block_ra(block_var=block_var, prob_each=c(.1, .1, .8))
#' table(block_var, Z)
block_ra <- function(block_var, num_arms= NULL, block_m=NULL, prob_each=NULL, condition_names = NULL){
  
  if(!is.null(block_m) & !is.null(prob_each)){
    stop("Do not specify both block_m and prob_each at the same time.")      
  }

  if(all(!is.null(num_arms), !is.null(prob_each), num_arms != length(prob_each))){
    stop("If both num_arms and prob_each are specified, num_arms must be equal to the length of prob_each")
  }
  
  if(all(!is.null(num_arms), !is.null(block_m), num_arms != ncol(block_m))){
    stop("If both num_arms and block_m are specified, num_arms must be equal to the number of columns of block_m")
  }
  
  
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
  assign <- rep(NA, length(block_var))
  
  if(is.null(block_m) & is.null(prob_each) & is.null(num_arms)){
    for(i in 1:length(blocks)){
      N_block <- sum(block_var==blocks[i])
      assign[block_var==blocks[i]] <- complete_ra(N = N_block, condition_names=condition_names)
    }
    return(assign)
  }
  
  if(is.null(block_m) & is.null(prob_each) & !is.null(num_arms)){
    for(i in 1:length(blocks)){
      N_block <- sum(block_var==blocks[i])
      assign[block_var==blocks[i]] <- complete_ra(N = N_block, num_arms=num_arms, condition_names=condition_names)
    }
    return(assign)
  }

  if(!is.null(block_m)){
    for(i in 1:length(blocks)){
      if(nrow(block_m)!=length(unique(blocks))){
        stop("block_m should have the same number of rows as there are unique blocks in block_var")
      }
      N_block <- sum(block_var==blocks[i])
      assign[block_var==blocks[i]] <- complete_ra(N = N_block, m_each = block_m[i,], condition_names=condition_names)
    }
    return(assign)
  }
  
  if(!is.null(prob_each)){
    
    for(i in 1:length(blocks)){
      if(sum(prob_each)!=1){
        stop("prob_each must sum to 1.")
      }
      N_block <- sum(block_var==blocks[i])
      assign[block_var==blocks[i]] <- complete_ra(N = N_block, prob_each = prob_each, condition_names=condition_names)
    }
    return(assign)
  }
}
