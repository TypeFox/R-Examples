#function that carries out StructuredHolm on a DAGstructure object
#optimization = none/LP/ILP
#isadjusted T if adjusted p-values required
structuredHolm <- function(DAGstructure, test, alpha_max = 0.05, isadjusted = FALSE, optimization = "none", pvalues = NULL, verbose = FALSE)
{
  
  #check whether optimization is allowed
  if(optimization != "none" && !DAGstructure@twoway)
    stop("Your DAG structure does not allow optimization, because it has no twoway relations.")
  
  if(!(optimization %in% c("none","ILP", "LP")))
    stop("Choose a valid argument for optimization!")
  
  #choose LPsolver for later use
  if(optimization != "none")
    solveLP <- getLPsolver()
  
  sets <- DAGstructure@sets
  parents <- DAGstructure@parents
  children <- DAGstructure@children
  twoway <- DAGstructure@twoway
  
  n <- length(sets)  
  if(n<1)
    stop("no sets to test")
  
  if(length(pvalues)==0)
  {
    pvalues <- rep(NA,n)
    for(i in 1:n)
      pvalues[i] <- test(sets[[i]])    
  }
  
  #make sets in terms of leaf nodes, the so-called "leaf_based_sets", needed for implications
  leaf_based_sets <- set_in_leaves(parents,children)
  
  rejected <- rep(FALSE,n)
  adj_pvalues <- rep(Inf,n)
  
  #sort pvalues from small to large
  index <- order(pvalues) #index of sorted element in original pvalues (and sets parents and children)
  sorted <- sort(pvalues)
  
  if(sorted[1] >= alpha_max/n)
    stop("Nothing to reject on this maximal alpha-level")
  
  num_rejected <- 0
  index_p <- 1 #position smallest unrejected pvalue
  alpha <- ifelse(isadjusted,0,alpha_max) #if no adj.pvalues needed: immediately alpha_max
  implications <- rep(FALSE,n) 
  rej_ILP <- 0 #minimal number false given by ILP 
  
  #as long as new rejections possible... 
  while(num_rejected < n && alpha <= alpha_max)
  {
    
    #if possible to reject the smallest p-value on the current alpha level
    if(pvalues[index[index_p]]*(n-max(num_rejected,rej_ILP)) <= alpha)
    {
      rejected[index[index_p]] <- TRUE
      adj_pvalues[index[index_p]] <- alpha
      rej_and_adj <- propagate(c(index[index_p]), parents, rejected, adj_pvalues, implications)
      rejected <- rej_and_adj$rejected
      adj_pvalues <- rej_and_adj$adj_pvalues
      implications <- rej_and_adj$implications
      
#       for(i in 1:length(implications))
#       {
#         if(implications[i])
#         {
#           #stay implication if you are no ancestor of new rejection 
#           if(all(sets[[index[index_p]]] %in% sets[[i]])) #ancestor
#             implications[i] <- FALSE        
#         }
#       }
#       #new rejection always is implication
#       implications[index[index_p]] <- TRUE 
      
    }
      
    num_rejected <- sum(rejected)
    
    #look for smallest unrejected pvalue
    while(index_p < n && rejected[index[index_p]])
      index_p <- index_p + 1
    
    #what alpha is needed to reject this pvalue? 
    #subtract largest possible number, found with or without optimization
    #i.e find out how many unrejected hypotheses can maximally be simultaneously true
          
    adj <- pvalues[index[index_p]]*(n-num_rejected)
    
    if(num_rejected > 0 && optimization != "none" && (isadjusted || adj > alpha_max)) 
    {
      impl <- leaf_based_sets[implications]
      idx <- which(implications) #index of implications
      rej_ILP <- ILPholm(solveLP, idx, impl,parents,relaxation=(optimization=="LP"))      
      adj <- pvalues[index[index_p]]*(n-rej_ILP)
    }

    alpha <- max(alpha, adj)
    
    
    if(verbose)
    {
      cat(sprintf("\r#rejections = %d.", num_rejected))
      flush.console()
    }
    
  }
  
  
  names(leaf_based_sets) <- names(sets)
  names(adj_pvalues) <- names(sets)
  names(implications) <- names(sets)
  names(rejected) <- names(sets)
  
  #NB: implications etc here just logical vectors. functions on DAG object can return what you really want
  out <- new("DAG",
             alpha = alpha_max,
             sets = sets,
             leaf_based_sets = leaf_based_sets,
             allpvalues = adj_pvalues,
             implications = implications,
             isadjusted = isadjusted,
             rejected = rejected,
             method = "holm",
             twoway = DAGstructure@twoway)
  
  return(out)

}

#index newly rejected pvalues. reject parents if not yet rejected + give its adj.pvalue
propagate <- function(index_vector, parents, rejected, adj_pvalues, implications) 
{  
  n <- length(rejected)
  #initialize queue
  queue <- rep(0, n) 
  head <- tail <- 1
  
  #put indices of new rejections on queue
  for(i in index_vector)
  {
    queue[tail] <- i
    tail <- tail + 1 
    implications[i] <- TRUE #new rejection = implication (NB: can only be 1 new rejection each time)
  }
  
  #put parents on queue
  while(head < tail)
  {
    current <- queue[head]
    head <- head + 1
    
    #put unrejected parents on queue
    for(j in parents[[current]])
    {
      implications[j] <- FALSE #parent of new rejection can no longer be an implication
      
      if(!rejected[j])
      {
        queue[tail] <- j
        tail <- tail + 1
        
        #reject -> will only be on queue once
        rejected[j] <- TRUE
        #adj_pvalue = adj_pvalue child
        adj_pvalues[j] <- adj_pvalues[index_vector[1]] #all same adjusted p-value, take first one. 
        
      }
    }
  }
  return(list(rejected=rejected,adj_pvalues=adj_pvalues, implications=implications))
}



ILPholm <- function(solveLP, indices, implications, parents, relaxation)
{
  n <- length(parents)
  
  n_impl <- length(implications)
  #number of elements in parents is exactly the number of constraints given by these numbers..)
  n_constraints <- sum(sapply(parents, length))
  
  num_non_zero <- n_impl
  num_non_zero <- num_non_zero + sum(sapply(implications,length)) 
  num_non_zero <- num_non_zero + n_constraints*2
  
  row_index <- rep(NA,num_non_zero)
  col_index <- rep(NA,num_non_zero)
  val_index <- rep(NA,num_non_zero)
  
  counter <- 1
  
  #for each implication: sum over leaf elements >= 1
  for(i in 1:n_impl)
  {
    for(j in implications[[i]])
    {
      row_index[counter] <- i
      col_index[counter] <- j 
      val_index[counter] <- 1
      counter <- counter + 1
    }
  }
  
  index_rows <- n_impl + 1 #next row in A
  
  #for each implication: make sure value will be 1
  for(i in indices)
  {
    row_index[counter] <- index_rows
    col_index[counter] <- i
    val_index[counter] <- 1
    counter <- counter + 1
    index_rows <- index_rows + 1    
  }
  
  
  #for each child,parent pair: parent >= child
  for(i in 1:n)
  {
    if(!is.null(parents[[i]]))
    {
      for(j in parents[[i]]) #xi - x(parent j) <= 0
      {
        row_index[counter] <- index_rows
        col_index[counter] <- j
        val_index[counter] <- -1
        row_index[counter+1] <- index_rows
        col_index[counter+1] <- i
        val_index[counter+1] <- 1        
        counter <- counter + 2
        index_rows <- index_rows + 1
      }
    }
  }  
  
  #minimize sum over all nodes
  obj <- rep(1,n)
  #at least 1 rejection in each implication: >= 1, implication itself = 1, child <= parent: child - parent <= 0
  rhs <- c(rep(1, n_impl) , rep(1, n_impl), rep(0, n_constraints))
  #>= for all implications , = for each implication itself, <= for all other constraints
  dir <- c(rep('>=', n_impl) , rep('=', n_impl) , rep('<=', n_constraints))
  #minimize objective function
  minmax <- "min"
  
  result <- solveLP(col_index, row_index, val_index, obj, rhs, dir, relaxation, minmax)
  
  #in case of relaxation: allowed to round of to highest integer (given that all weights are natural numbers).
  #in case of ILP: ceiling will not alter the result
  #if ILP gives not exact numbers, rounding goes still well
  return(ifelse(isTRUE(all.equal(result, round(result))), result, ceiling(result)))
}
    