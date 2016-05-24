#
# determine adjusted pvalues taking into account the logical structure of the region hypotheses
#

setClass("region",
         representation(
           alpha = "numeric",           # stores chosen alpha for testing
           allpvalues = "matrix",       # stores (adjusted) p-values
                                        #  (has value NA if adjusted p > alpha)
                                        #  (0 x 0 matrix if no p-values calculated)
           implications = "matrix",     # stores implications including (adjusted) pvalues at chosen alpha
           isadjusted = "logical",      # T if adjusted p-values are calculated, F otherwise
           weights = "numeric",         # stores chosen weights
           totalrejected = "numeric"    # stores total number of rejected region hypotheses
         )
)

LEFT        <- 1
RIGHT       <- 2
RATIO       <- 3
PVALUE      <- 4
CAND_FIELDS <- 4
ADJ_PVALUE  <- 3
IMPL_FIELDS <- 3


# region method determines which region hypotheses can be rejected on a chosen alpha-level while adjusting for multiple testing
# test is a function that accepts two parameters, left and right.
# test(left,right) calculates a pvalue by testing the hyps from left up to and including right.
#
# TODO: default weights: impossible -> user has to specify num_hyps in that case
# TODO: pvalues op 1, iets met alpha_max = NA ? For now: alpha_max always set to 1 if a higher value is chosen
# hypotheses that cannot get rejected on this level are given the value NA (which thus means, not rejected on alpha_max)
regionmethod <- function(weights, test, alpha_max = 0.05, all_pvalues = FALSE, isadjusted = FALSE, verbose = FALSE)
{

  if(alpha_max > 1) {
    alpha_max <- 1
    warning("The alpha-level exceeds 1. We have set it to 1.")
  }
  
  num_hyps <- length(weights)
  if(num_hyps <= 0)
    stop("No weights have been given!")
  if(any(weights <= 0))
    stop("Not all weights are strictly positive!")
  
  # determine partial sum of weights
  partial_weights <- rep(NA, num_hyps)
  partial_weights[1] <- weights[1]
  if(num_hyps > 1)
    for(i in 2:num_hyps)
      partial_weights[i] <- partial_weights[i-1] + weights[i]
 
  if(all_pvalues)
    adj_pvalues <- matrix(NA, num_hyps, num_hyps)
  else
    adj_pvalues <- matrix(NA, 0, 0)

  
  # add root as cand
  cands <- matrix(c(1, num_hyps, 1, test(1, num_hyps)), 1, CAND_FIELDS)
  impls <- matrix(NA, 0, IMPL_FIELDS)
  
  # record some statistics
  total_rejected  <- 0
  total_intervals <- num_hyps * (num_hyps + 1) / 2
  num_iters  <- 0
  
  # try to reject cands while possible
  alpha <- ifelse(isadjusted, 0, alpha_max)
  while(alpha <= alpha_max && nrow(cands) > 0)
  {
    # calculate ratios
    cands <- calculate_all_ratios(cands, impls, weights, partial_weights)
    
    # find new alpha
    min_new_alpha <- Inf
    for(i in 1:nrow(cands))
    {
      min_new_alpha <- min(min_new_alpha, cands[i,PVALUE] / cands[i,RATIO])
    }
    
    if(alpha < min_new_alpha)
    {
      alpha <- min_new_alpha
    }
    
    if(alpha <= alpha_max)
    {
      # try to reject cands    
      num_rejected <- 0
      for(i in 1:nrow(cands))
      {
        if(is_rejected(cands[i,],alpha))
        {
          num_rejected <- num_rejected + 1
          if(all_pvalues)
          {
            adj_pvalues[cands[i,LEFT], cands[i,RIGHT]] <- alpha
          }
          
        }
      }
      
      #NB: order matters! (update_impls uses not updated cands)
      impls <- update_impls(impls, cands,alpha, num_hyps)
      cands <- update_cands(cands,alpha, num_hyps,test)
      
      total_rejected <- total_rejected + num_rejected
      
      num_iters <- num_iters + 1
      if(verbose)
      {
        cat(sprintf("\r#rejections = %d.", total_rejected))
        flush.console()
      }
    }
    
  }
  
  # set remaining adj pvalues of unreached regions to 1
#   if(all_pvalues) {
#     for(i in 1:num_hyps)
#       for(j in i:num_hyps)
#         if(is.na(adj_pvalues[i,j]))
#           adj_pvalues[i,j] <- 1
#   }
  
  #make sure output is a matrix 
  #probably unnecessary, because it seems that implications get only modified through update_impls, which already makes sure it will be a matrix
  if(nrow(impls) > 0)
  {
    impls <- impls[(1:nrow(impls)), , drop=FALSE]
  }
  
  out <- new("region",
             alpha = alpha_max,
             allpvalues = adj_pvalues,
             implications = impls,
             isadjusted = isadjusted,
             weights = weights,
             totalrejected = total_rejected)
  
  return(out)
}


calculate_all_ratios <- function(cands, impls, weights, partial_weights) {
  
  dp_fw <- DP_FW(weights, impls)
  dp_bw <- DP_BW(weights, impls)
  
  num_cands <- nrow(cands)
  
  for(i in 1:num_cands)
  {
    ratio <- get_ratio(cands[i,], partial_weights, dp_fw, dp_bw)
    if(ratio < cands[i,RATIO])
    {
      warning("Monotonicity condition violated!")
    }
    
    cands[i, RATIO] <- ratio
  }
  return(cands)
}

#dp_fw[i] is minimal weight needed to satisfy al implications starting at or before i, while rejecting i
DP_FW <- function(weights, impls)
{

  num_hyps <- length(weights)
  
  dp_fw <- rep(NA, num_hyps)
  queue <- rep(NA, num_hyps)  
  
  num_impls <- nrow(impls)
  
  if(num_impls > 0)
  {
    index <- -Inf
    prev_min <- 0
    qhead <- 1
    qtail <- 0
    
    for(i in 1:num_impls)
    {
      if(index < impls[i,LEFT])
        index <- impls[i,LEFT]
      
      # remove indices belonging to the previous implication from consideration while not empty
      while(qhead <= qtail && queue[qhead] < impls[i,LEFT])
      {
        qhead <- qhead + 1
      }
      
      # reject current implication by rejecting hypothesis at index
      while(index <= impls[i,RIGHT])
      {  
        dp_fw[index] <- weights[index] + prev_min
        
        # remove indices that will never point to the minimum
        while(qhead <= qtail && dp_fw[queue[qtail]] >= dp_fw[index])
        {
          qtail <- qtail - 1
        }
        
        # add current index
        qtail <- qtail + 1
        queue[qtail] <- index
        
        index <- index + 1
      }
      
      # record minimum of current implication
      prev_min <- dp_fw[queue[qhead]]
    }
  }  
  
  return (dp_fw)  
  
}

#dp_bw[i] is minimal weight needed to satisfy al implications ending at or after i, while rejecting i
DP_BW <- function(weights, impls)
{

  num_hyps <- length(weights)
  
  
  dp_bw <- rep(NA, num_hyps)
  queue <- rep(NA, num_hyps)
  
  num_impls <- nrow(impls)
  
  if(num_impls > 0)
  {
    index <- Inf
    prev_min <- 0
    qhead <- 1
    qtail <- 0
    
    for(i in num_impls:1)
    {
      if(index > impls[i,RIGHT])
        index <- impls[i,RIGHT]
      
      # remove indices belonging to the previous implication from consideration while not empty
      while(qhead <= qtail && queue[qhead] > impls[i,RIGHT])
      {
        qhead <- qhead + 1
      }
      
      # reject current implication by rejecting hypothesis at index
      while(index >= impls[i,LEFT])
      {  
        dp_bw[index] <- weights[index] + prev_min
        
        # remove indices that will never point to the minimum
        while(qhead <= qtail && dp_bw[queue[qtail]] >= dp_bw[index])
        {
          qtail <- qtail - 1
        }
        
        # add current index
        qtail <- qtail + 1
        queue[qtail] <- index
        
        index <- index - 1
      }
      
      # record minimum of current implication
      prev_min <- dp_bw[queue[qhead]]
    }
  }
  
  return (dp_bw)
  
}


#
# determines minimum ratio for given candidate 
# note that the leaves directly surrounding the candidate have to be rejected!!
#
get_ratio <- function(cand, partial_weights, dp_fw, dp_bw)
{

  num_hyps <- length(partial_weights)
  
  #add rejected weights left and right from candidate, if any
  weight <- 0
  if(cand[LEFT] > 1)
    weight <- weight + dp_fw[cand[LEFT] - 1]
  if(cand[RIGHT] < num_hyps)
    weight <- weight + dp_bw[cand[RIGHT] + 1]
  
  if(cand[LEFT] > 1)
    return((partial_weights[cand[RIGHT]] - partial_weights[cand[LEFT]-1]) / (partial_weights[num_hyps] - weight))
  else
    return(partial_weights[cand[RIGHT]] / (partial_weights[num_hyps] - weight))
}  


#
# indicates whether given candidate is rejected
#
is_rejected <- function(cand,alpha)
{
  return(cand[PVALUE] / cand[RATIO] <= alpha)
}


# TODO: rewrite update candidates and implications to not reuse memory!

#
# update candidates such that the result remains sorted on the left bounds
#
update_cands <- function(cands,alpha,num_hyps,test)
{
   
  num_cands <- nrow(cands)
  new_cands <- matrix(NA, num_hyps, CAND_FIELDS)
  num_new_cands <- 0
  
  # we know that num_cands > 0
  for(i in 1:num_cands)
  {
    if(!is_rejected(cands[i, ],alpha))
    {
      # keep current candidate
      num_new_cands <- num_new_cands + 1
      new_cands[num_new_cands, ] <- cands[i, ]
    }
    else if(cands[i, LEFT] < cands[i, RIGHT]) # if no leaf
    {
      # add left child if all its parents are rejected    
      if(i == 1 || cands[i-1, RIGHT] <  cands[i, RIGHT]-1
         ||(cands[i-1, RIGHT] == cands[i, RIGHT]-1 && cands[i-1,LEFT] ==  cands[i,LEFT] -1 && is_rejected(cands[i-1,],alpha))
      )
      {
        # check if left child is not already added
        if(num_new_cands == 0 || cands[i,  LEFT]     != new_cands[num_new_cands,  LEFT] 
           || cands[i, RIGHT] - 1 != new_cands[num_new_cands, RIGHT])
        {
          num_new_cands <- num_new_cands + 1
          new_cands[num_new_cands,     LEFT] <- cands[i,  LEFT]
          new_cands[num_new_cands,    RIGHT] <- cands[i, RIGHT] - 1
          new_cands[num_new_cands,    RATIO] <- 0
          new_cands[num_new_cands,   PVALUE] <- test(new_cands[num_new_cands, LEFT], new_cands[num_new_cands, RIGHT])
        }
      }
      
      # add right child if all its parents are rejected 
      if(i == num_cands ||  cands[i+1, LEFT] >  cands[i, LEFT]+1
         || (cands[i+1,RIGHT] == cands[i,RIGHT]+1 && cands[i+1,LEFT] ==  cands[i,LEFT]+1 && is_rejected(cands[i+1,],alpha))
      )
      {
        # right child cannot have been added before, so no need to check
        num_new_cands <- num_new_cands + 1
        new_cands[num_new_cands,     LEFT] <- cands[i,  LEFT] + 1
        new_cands[num_new_cands,    RIGHT] <- cands[i, RIGHT]
        new_cands[num_new_cands,    RATIO] <- 0
        new_cands[num_new_cands,   PVALUE] <- test(new_cands[num_new_cands, LEFT], new_cands[num_new_cands, RIGHT])
      }
    }
  }
  
  if(num_new_cands == 0)
    return(matrix(NA, 0, CAND_FIELDS))
  else
    return(new_cands[(1:num_new_cands), , drop=FALSE])
}


#
# update impls such that the result remains sorted on the left bounds
# note that there is no need to check whether an impl has already been added
#
update_impls <- function(impls, cands,alpha, num_hyps)
{
  num_cands <- nrow(cands)
  num_impls <- nrow(impls)
  
  new_impls <- matrix(NA, 2*num_hyps, IMPL_FIELDS)
  
  num_new_impls <- 0
  impl_index    <- 1
  cand_index    <- 1
  
  # merge rejected cands with previous impls
  while(TRUE)
  {
    # find next rejection
    while(cand_index <= num_cands && !is_rejected(cands[cand_index,],alpha))
    {
      cand_index <- cand_index + 1
    }
    
    # stop if no new impls found, otherwise add the leftmost new impl
    if(impl_index > num_impls && cand_index > num_cands)
    {
      break
    }
    else if(impl_index > num_impls)
    {
      num_new_impls <- num_new_impls + 1
      new_impls[num_new_impls,       LEFT] <- cands[cand_index,  LEFT]
      new_impls[num_new_impls,      RIGHT] <- cands[cand_index, RIGHT]
      new_impls[num_new_impls, ADJ_PVALUE] <- alpha
      cand_index <- cand_index + 1
    }
    else if(cand_index > num_cands)
    {
      num_new_impls <- num_new_impls + 1
      new_impls[num_new_impls,       LEFT] <- impls[impl_index,       LEFT]
      new_impls[num_new_impls,      RIGHT] <- impls[impl_index,      RIGHT]
      new_impls[num_new_impls, ADJ_PVALUE] <- impls[impl_index, ADJ_PVALUE]
      impl_index <- impl_index + 1
    }
    else if(cands[cand_index, LEFT] <= impls[impl_index, LEFT])
    {
      num_new_impls <- num_new_impls + 1
      new_impls[num_new_impls,       LEFT] <- cands[cand_index,  LEFT]
      new_impls[num_new_impls,      RIGHT] <- cands[cand_index, RIGHT]
      new_impls[num_new_impls, ADJ_PVALUE] <- alpha
      cand_index <- cand_index + 1
    }
    else
    {
      num_new_impls <- num_new_impls + 1
      new_impls[num_new_impls,       LEFT] <- impls[impl_index,       LEFT]
      new_impls[num_new_impls,      RIGHT] <- impls[impl_index,      RIGHT]
      new_impls[num_new_impls, ADJ_PVALUE] <- impls[impl_index, ADJ_PVALUE]
      impl_index <- impl_index + 1
    }
  }
  
  # delete impls with rejected children
  size <- 0
  impl_index <- 1
  for(i in 1:num_new_impls)
  {
    if(
      # test whether predecessor is not your child
      (i == 1 || new_impls[i,  LEFT]   != new_impls[i-1,  LEFT]
       || new_impls[i, RIGHT]-1 != new_impls[i-1, RIGHT])
      # test whether successor is not your child
      && (i == num_new_impls || new_impls[i,  LEFT]+1 != new_impls[i+1,  LEFT]
          || new_impls[i, RIGHT]   != new_impls[i+1, RIGHT] )
    )
    {
      size <- size + 1
      new_impls[size, ] <- new_impls[i, ]
    }
  }
  
  if(size == 0)
    return(matrix(NA, 0, IMPL_FIELDS))
  else
    return(new_impls[(1:size), , drop=FALSE])
}



#given matrix with adjusted pvalues: find implications for alpha level alpha
matrix_to_implications <- function(matrix, alpha)
{
  num_hyps <- nrow(matrix)
  
  num_impls <- 0
  i <- 1
  j <- 1
  
  #initialize matrix for implications (left,right,adj pvalue)
  impls <- matrix(0, num_hyps,3)
  
  while(j <= num_hyps && i<= num_hyps)
  {
    #find first j satisfying condition, given i
    while(j <= num_hyps && (is.na(matrix[i,j]) || matrix[i,j] > alpha)) #NB: order of logical statments matters
    {
      j <- j + 1
    }
    if(j <= num_hyps) #implication in this column
    {
      num_impls <- num_impls + 1
      impls[num_impls,2] <- j
      
      #finds first i not satisfying condition, given j
      while(i <= num_hyps && (!is.na(matrix[i,j]) && matrix[i,j] <= alpha)) 
      {
        i <- i + 1
      }
      #implication has previous i-value
      impls[num_impls,1] <- i-1
      impls[num_impls,3] <- matrix[i-1,j]
    }
  }
  
  if(num_impls>0)
    return(impls[(1:num_impls), , drop=FALSE])
  else
    return(matrix(NA, 0, 3))
}

#function that gives back the implications corresponding to the chosen alpha-level, if these implications are constructable from the region object
#uses matrix_to_implications
get_implications <- function(region, alpha) 
{
  # check if given alpha is valid and recompute implications if necessary
  if(alpha > region@alpha) {
   stop("The region procedure has never been carried out on this alpha-level. Choose a lower level.")
  }
   
  if(alpha < region@alpha) {
   if(nrow(region@allpvalues) == 0 || !region@isadjusted)
     stop("The region object does not contain adjusted p-values for all possible region hypotheses.")
   else #implications has to be changed to this new alpha-level
     matrix_to_implications(region@allpvalues, alpha)
  }
  else
   region@implications
}


#function to find minimal number of hypotheses/minimal weight that has to be rejected in specified set, given the implications at a certain alpha
#specified set is a list containing intervals (the user can choose the intervals freely, they don't have to be ordered and can even overlap)
# TODO: alpha standaard op region@alpha zetten? 
regionpick <- function(region, intervals, alpha, silent = FALSE, ignore_weights = TRUE)
{
  if(length(intervals) == 0)
    stop("No intervals specified.")
  
  weights <- region@weights
  num_hyps <- length(weights)
  
  if(ignore_weights)
  {
    weights <- rep(1,num_hyps)
  }
  
  #if no alpha specified, take the alpha stored in the region object
  if (missing(alpha))
    alpha <- region@alpha
  
  impls <- get_implications(region, alpha)
  
  output <- function(rejected_mass) {
    if (!silent) {
      cat("At confidence level ", 1-alpha, ": ", sep="")
      if (ignore_weights)
        cat("False null-hypotheses >= ", rejected_mass, "\n", sep="")
      else
        cat("Minimal weight of false null-hypotheses >= ", rejected_mass, "\n", sep="")
      invisible(rejected_mass)
    } else
      rejected_mass
  }
  
  num_impls <- nrow(impls)
  if(num_impls == 0) {
    return(output(0))
  }
  
  
  # check for wrong intervals
  for(interval in intervals)
    if(interval[LEFT] > interval[RIGHT] || interval[LEFT] < 1 || interval[RIGHT] > num_hyps) 
      stop(sprintf("[%d, %d] is not a valid region", interval[LEFT], interval[RIGHT]))
  
  # sort intervals on increasing left bound
  intervals <- intervals[order(as.numeric(lapply(intervals, function(interval) interval[LEFT])))]
  
  # merge overlapping or touching intervals
  new_intervals <- vector("list", length(intervals))
  size <- 0
  
  index <- 1
  while(index <= length(intervals))
  {
    left <- intervals[[index]][LEFT]
    right <- intervals[[index]][RIGHT]
    
    index <- index + 1
    
    while(index <= length(intervals) && intervals[[index]][LEFT] <= (right+1)) 
    {
      right <- max(right, intervals[[index]][RIGHT])
      index <- index + 1 
    }
    
    size <- size + 1
    new_intervals[[size]] <- c(left, right)
  }
  intervals <- new_intervals[1:size]

  
  num_intervals <- length(intervals)
  num_impls_new <- 0
  impls_new <- matrix(NA,num_impls,2)
  
  
  #check for every implication whether there is an interval in which is it contained. if such an interval is found, go to next implication. 
  for(i in 1:num_impls)
  {
    for(j in 1:num_intervals)
    {
      #if impl i contained in interval j, impl i has to be satisfied within chosen interval
      if(impls[i,LEFT] >= intervals[[j]][LEFT] && impls[i,RIGHT] <= intervals[[j]][RIGHT])
      {
        num_impls_new <- num_impls_new + 1
        impls_new[num_impls_new,] <- impls[i,(LEFT:RIGHT)]
        break
      }  
    }
  }
  
  num_impls <- num_impls_new
  
  #only continue if there are still some impls left
  if(num_impls > 0)
  {
    impls <- impls_new[1:num_impls, , drop=FALSE]
    dp_fw <- DP_FW(weights,impls)
  
    min <- Inf
    
    #find smallest weight needed to reject all implications that are fully contained in the specified set
    for(i in impls[num_impls,LEFT]:impls[num_impls,RIGHT])
    {
      if(dp_fw[i] < min)
      {
        min <- dp_fw[i]
      }  
    }
    
    #min number of rejected hyps/min number of weights
    rejected_mass <- min
  }
  else
    rejected_mass <- 0
  
  
  output(rejected_mass)
}


setMethod("show", "region", function(object) {
  num_hyps <- length(object@weights)
  cat("Region method result on", num_hyps, "elementary hypotheses.\n")
  cat("There are ", object@totalrejected, " region hypotheses rejected out of a total of ", num_hyps*(num_hyps+1)/2, "\n", "at an alpha-level of ", object@alpha, ".\n", sep="")
})

#evt later veranderen. gaat over standaard summary heen? 
setMethod("summary", "region", function(object) {
  num_hyps <- length(object@weights)
  cat("Region method result on", num_hyps, "elementary hypotheses.\n")
  cat("There are ", object@totalrejected, " region hypotheses rejected out of a total of ", num_hyps*(num_hyps+1)/2, "\n", "at an alpha-level of ", object@alpha, ".\n", sep="")
})


#function that retrieves the implications from a given region object. 
#The alpha level can be varied, if the supporting information is present in the region object.
# waarom setGeneric??
setGeneric("implications", function(object, ...) standardGeneric("implications"))
setMethod("implications", "region", function(object, alpha) {
  if (missing(alpha))
    alpha <- object@alpha  
  
  impls <- get_implications(object, alpha)
#  impllist <- vector("list", nrow(impls))
#  if(nrow(impls) > 0)
#    for(i in 1:nrow(impls))
#      impllist[[i]] <- impls[i,]
#  impllist
  colnames(impls) <- c("left", "right", ifelse(object@isadjusted, "adj. pvalue", "pvalue"))
  impls
})

#function to retrieve the maximal alpha_value from the region object
setMethod("alpha", "region", function(object) {
  object@alpha
})


#function to retrieve pvalues for all possible region hypotheses (indicated with left and right bound) from the region object
#only able to return the value if the allpvalues matrix is stored in the region object
setMethod("pvalue", "region", function(object, left, right) {
  if(nrow(object@allpvalues) == 0)
    stop("The region object does not contain (adjusted) p-values for all possible region hypotheses.")
  
  # check for wrong interval
  if(left > right || left < 1 || right > length(object@weights)) 
      stop(sprintf("[%d, %d] is not a valid region", left, right))
  
  object@allpvalues[left,right]
})



#function that turns set of implications into a polygon that then can be plotted
#cannot handle empty sets of implications, but the function that calls this one should take care of that already.
impls_to_polygon <- function(num_hyps, impls)
{
  #num_impls > 0 
  num_impls <- nrow(impls)
  
  all_coords <- matrix(NA,num_impls*2+2,2)

  to_coord <- function(impl)
  {
    return(c((impl[LEFT]+impl[RIGHT])/2, (impl[RIGHT] - impl[LEFT]) + 1))
  }
  
  #left boundary. 
  intersectionleft <- c(1,impls[1,RIGHT])
  all_coords[1,] <- to_coord(intersectionleft)
  
  #first implication
  all_coords[2,] <- to_coord(impls[1,])
  
  #right boundary.   
  intersectionright <- c(impls[num_impls,LEFT],num_hyps)
  all_coords[num_impls*2+1,] <- to_coord(intersectionright)
  
  #top
  all_coords[num_impls*2+2,] <- to_coord(c(1,num_hyps))
  
  index <- 3
  
  if(num_impls>1)
  {
    for(i in 1:(num_impls-1))
    {
      intersection <- c(impls[i,LEFT], impls[i+1,RIGHT])
      all_coords[index,] <- to_coord(intersection)
      all_coords[index+1,] <- to_coord(impls[i+1,])
      index <- index + 2
    }                                    
    
  }
  
  return(list(x=all_coords[,1], y=all_coords[,2]))
  
}

#function that plots the rejected hypotheses, by drawing a polygon that follows the original graph structure.
#TODO (ever?) instead of fully coloring the polygon, draw the underlying graph (the "diamondpattern")
regionplot <- function(region, alpha, color="red")
{
  #if no alpha specified, take the alpha from the region object
  if (missing(alpha))
    alpha <- region@alpha  
  
  num_hyps <- length(region@weights)
  
  impls <- get_implications(region, alpha)
  
  plot(x=NULL,y=NULL, type="l", col="darkgrey", axes=FALSE, xlim=c(1,num_hyps), ylim=c(1,num_hyps), ylab="Region Length", xlab="" ,frame.plot=FALSE, cex.lab=1.2)
  #axis(side=2)
  axis(side=2, at=round(seq(from=1, to=num_hyps, length.out=5)), cex.axis=1.2)
  
  if(nrow(impls) == 0)
    warning("No regions got rejected at this alpha level.")
  else
    poly1 <- impls_to_polygon(num_hyps,impls)

  polygon(x=c((num_hyps+1)/2,1,num_hyps), y=c(num_hyps,1,1), lty="dashed")
  
  if(nrow(impls) > 0)
    polygon(x=poly1$x, y=poly1$y, col=color, border=color, lwd=2)
  
  par(xpd=T)
  legend("bottom", inset=c(0,-.15), legend=c("rejected hypotheses", "unrejected hypotheses"), fill=c(color, "white"), cex=1.1)
}

# NB: werkt ook voor 0 implications
regionplot2 <- function(region,alpha,color_rej="red",color_unrej="grey") 
{
  #if no alpha specified, take the alpha from the region object
  if (missing(alpha))
    alpha <- region@alpha  

  impls <- get_implications(region, alpha)
  num_impls <- nrow(impls)
  num_hyps <- length(region@weights)
  
  # TODO: dikte waarschijnlijk aanpassen!
  drawsegment <- function(left1, right1, left2, right2, color)
  {
    if(left1!=left2 || right1!=right2)      
      segments((left1+right1)/2, right1-left1+1, 
               (left2+right2)/2, right2-left2+1, col=color)        
  }

  # setup plot
  plot(x=NULL,y=NULL, type="l", col="darkgrey", axes=FALSE, xlim=c(1,num_hyps), ylim=c(1,num_hyps), ylab="Region Length", xlab="" ,frame.plot=FALSE, cex.lab=1.2)
  axis(side=2, at=round(seq(from=1, to=num_hyps, length.out=5)), cex.axis=1.2)
  
  
  # create left leaning line segments
  left <- 1
  index <- 0
  
  for(hyp in 1:num_hyps) {
    if(index+1 <= num_impls && impls[index+1, RIGHT] == hyp) {
      index <- index+1
      left <- impls[index, LEFT]
    }
    
    drawsegment(hyp, hyp, left, hyp, color_unrej)
    drawsegment(left, hyp, 1, hyp, color_rej)
  }
  
  # create right leaning line segments
  right <- num_hyps
  index <- num_impls+1
  
  for(hyp in num_hyps:1) {
    if(index-1 >= 1 && impls[index-1, LEFT] == hyp) {
      index <- index-1
      right <- impls[index, RIGHT]
    }
    
    drawsegment(hyp, hyp, hyp, right, color_unrej)
    drawsegment(hyp, right, hyp, num_hyps, color_rej)
  }
  
  par(xpd=T)
  legend("bottom", inset=c(0,-.15), legend=c("rejected hypotheses", "unrejected hypotheses"), fill=c(color_rej, color_unrej), cex=1.1)
    
}

# 
# #function that plots region object 
# setMethod("plot", "region", regionplot)




