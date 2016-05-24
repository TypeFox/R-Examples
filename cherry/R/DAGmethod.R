#TODO: at loading library(cherry), give a message that gurobi would be handy 

# TODO
# potentially. add leaves to make oneway structure into twoway structure 
# later: add weight function (function that assigns weight to leaf nodes (allthough difficult because leaf nodes are unknown to user))

# TODO: allow alpha parameter at implications method


#DAG object 
setClass("DAG",
         representation(
           alpha = "numeric",           # stores chosen alpha for testing
           sets = "list",               # stores unique original sets
           leaf_based_sets = "list",    # stores sets expressed in leaves
           allpvalues = "numeric",      # stores (adjusted) p-values
                                        # (has value NA if adjusted p > alpha)
           implications = "logical",    # stores whether sets are implications at chosen alpha
           isadjusted = "logical",      # T if adjusted p-values are calculated, F otherwise
           rejected = "logical",        # stores whether sets are rejected
           method = "character",        # stores whether any-, all-parent or the structuredHolm method has been used
           twoway = "logical"           # stores whether underlying DAG structure had two-way logical relationships
         )
)


#TODO; try to remove check for length(implications)>0 

#multiple testing method for DAG structures, two variants: any-parent and all-parents method
#three options: adj/unadj, optimization: none/LP/ILP, degree: group/individual
DAGmethod <- function(DAGstructure, test, alpha_max = 0.05, method = "all", isadjusted = FALSE, optimization = "none", degree = "group", pvalues = NULL, verbose = FALSE)
{
  #check whether optimization is allowed
  if(optimization != "none" && !DAGstructure@twoway)
    stop("Your DAG structure does not allow optimization, because it has no twoway relations.")
  
  #choose LPsolver for later use
  if(optimization != "none")
    solveLP <- getLPsolver()
  
  sets <- DAGstructure@sets
  parents <- DAGstructure@parents
  children <- DAGstructure@children
  twoway <- DAGstructure@twoway
  
  n <- length(sets)
  
  #initialize where to start
  if(isadjusted)
  {
    cur_degree <- degree
    cur_optim <- optimization    
  }
  else
  {
    cur_degree <- "group"
    cur_optim <- "none"
  }
  
  #make sets in terms of their leaf nodes: so called "leaf_based_sets"
  leaf_based_sets <- set_in_leaves(parents,children)
  
  
  #
  # commence algorithm
  #
  
  
  num_rejected <- 0
  rejected <- rep(FALSE, n)
  to_be_rejected <- rep(FALSE,n) #only needed for simultaneous updating numerator and denominator (gives no adj pvalues), but has to be FALSE in other cases.. (needed in if statement)
  if(length(pvalues)==0) #not precalculated by the user
    pvalues <- rep(NA_real_, n)
  adj_pvalues <- rep(NA_real_, n) #NA if not rejected at alpha_max (also ensures that it is the numeric type)
  ratios <- rep(0,n)
  
  #optimize numerator and denominator simultaneously
  if(method == "any" && !isadjusted && degree == "individual" && optimization != "none")
  {
    max_weights <- update_max_weights(parents, children, rejected) #ingredient LPnumdenom 
  }
  
  
  #number of leaves, needed for calculating denuminator  
  num_leaves <- 0
  for(i in 1:n)
    if(length(children[[i]])==0)
      num_leaves <- num_leaves + 1
  
  implications <- rep(FALSE,n)
  
  #if adjusted p-values, update alpha incrementally, otherwise, start directly at alpha_max
  if(isadjusted)
    alpha <- 0
  else
    alpha <- alpha_max
  
  while(num_rejected < n) #still rejections possible
  {
    
    # update weight to find candidates (weight>0) and to find the numerator of the ratio 
    weights <- update_weights(parents, children, rejected, method)
    
    impl <- leaf_based_sets[implications]
    
    #without using logical relations
    if(cur_optim == "none")
    {
      unrejected_leaves <- 0
      for(i in 1:n)
        if(!rejected[i] && length(children[[i]]) == 0)
          unrejected_leaves <- unrejected_leaves + 1
    }
    else #some optimization
    {
      if(length(impl)>0)
      {
        #(I)LP per candidate group 
        if(cur_degree == "group")
          rejected_leaves_easy <- LPdenom(solveLP, impl, NULL, relaxation = (cur_optim == "LP"))      
      }
      #no implication yet
      else
        rejected_leaves_easy <- 0      
    }
    
    
    # find smallest alpha on which you can reject anything by looping over the candidates 
    min_alpha <- Inf
    
    for(i in 1:n) #for every candidate: calculate pvalue (if necessary) and ratio
    {
      if(weights[i] > 0) #candidate
      {
        
        if(is.na(pvalues[i]))
          pvalues[i] <- test(sets[[i]])
        
        #leaf_based_sets[i] or leaf_based_sets[[i]]? 
        if(pvalues[i] <= alpha_max) #only ratio for candidates with pvalue <= alpha_max (ratio always <=1)
        {
          #store ratios 
          if(cur_optim == "none") #without logical relations
          {
            ratios[i] <- weights[i]/unrejected_leaves
          }          
          else if(cur_degree == "group") #(I)LP for all candidates together
          {
            ratios[i] <- weights[i]/(num_leaves - rejected_leaves_easy)
          }          
          else if(cur_degree == "individual" && (method == "all" || isadjusted) )#(I)LP per candidate (skip if method == "any && !isadjusted)
          {
            if(length(impl)>0)
              rejected_leaves <- LPdenom(solveLP, impl,leaf_based_sets[[i]], relaxation = (cur_optim == "LP"))
            else #no impl yet
              rejected_leaves <- 0
            
            ratios[i] <- weights[i] / (num_leaves - rejected_leaves)
          }          
          
          min_alpha <- min(min_alpha, pvalues[i] / ratios[i])  #(*)
          
          if(method == "any" && !isadjusted && cur_degree == "individual" && cur_optim != "none") #optimizing numerator & denominator simultaneously, no adj.pval, only rejected yes/no
          {
            if(length(impl)>0)
            {
              par <- parents[[i]]
              num_par <- length(par)
              parentscand <- vector("list", num_par)
              counter <- 1
              for(j in par)
              {
                parentscand[[counter]] <- leaf_based_sets[[j]]
                counter <- counter + 1
              }
              
              g1 <- max_weights[i]/num_par
              palpha <- pvalues[i]/alpha_max
              res <- LPnumdenom(solveLP, impl, leaf_based_sets[[i]], parentscand, g1, palpha, relaxation = (cur_optim == "LP"))
              
              #if res >= palpha*g2, to_be_rejected! -> keep alpha at alpha_max
              #if you don't get in this loop, new_alpha will later be set to Inf
              if(res >= palpha*num_leaves)
              {
                to_be_rejected[i] <- TRUE
                min_alpha <- min(min_alpha, alpha_max)
              }
              
            }            
            else #no impl yet
            {
              rejected_leaves <- 0
              
              ratios[i] <- weights[i] / (num_leaves - rejected_leaves)
              min_alpha <- min(min_alpha, pvalues[i] / ratios[i])
            }
            
            
          }
          
          
        }
        
      }
      
    }
    
    new_alpha <- max(alpha, min_alpha)
    
    if(new_alpha > alpha_max) #no rejections possible on current degree and optimization level
    {
      if(cur_degree != degree || cur_optim != optimization)
      {
        cur_degree <- degree
        cur_optim <- optimization
        next #start while-loop again: on desired degree and optimization level
      }
      else
        break #leave while-loop: meaning termination of algorithm, report results
    }
    else #still rejections possible on current degree and optimization level
      alpha <- new_alpha
    
    
    # reject and update implications
    for(i in 1:n)
    {
      if(weights[i] > 0 && !rejected[i] && pvalues[i]<=alpha_max) #candidate with possiblity of rejection
      {
        #candidates that really have to be rejected, NB: order important because of rounding floating points, must be same as (*)
        if(alpha >= pvalues[i]/ratios[i] || to_be_rejected[i]) #second argument only used when LPnumdenom was used
        {
          if(method == "all")
          {
            rejected[i] <- TRUE
            #i implication, parents i can no longer be implications
            implications[i] <- TRUE
            for(j in parents[[i]])
              implications[j] <- FALSE
            
            num_rejected <- num_rejected + 1
            adj_pvalues[i] <- alpha
          }
          
          #in case of the any-parent method: make sure that possibly unrejected parents of rejected children get rejected themselves
          #check parents, if parent already rejected: done, else, go to their parents as well. 
          if(method == "any")
          {
            rejected[i] <- TRUE
            implications[i] <- TRUE 
            num_rejected <- num_rejected + 1
            adj_pvalues[i] <- alpha
            
            #initialize queue
            queue <- rep(0, n)
            head <- tail <- 1
            
            #put i on queue
            queue[tail] <- i
            tail <- tail + 1
            
            while(head < tail)
            {
              current <- queue[head]
              head <- head + 1
              
              #put unadjusted parents on queue
              for(j in parents[[current]])
              {
                if(!rejected[j]) #look at unrejected parents
                {
                  queue[tail] <- j
                  tail <- tail + 1
                  
                  #reject, increase num_rejected, adjust adj_pvalue
                  rejected[j] <- TRUE
                  num_rejected <- num_rejected + 1
                  adj_pvalues[j] <- alpha
                }
                else #for rejected parents, make sure it's no implication anymore (unrjected parents can not be an implication anyway)
                {
                  implications[j] <- FALSE
                }
                
              }
              
            }
          }          
          
          if(verbose)
          {
            cat(sprintf("\r#rejections = %d.", num_rejected))
            flush.console()
          }
          
        }    
        
      }
      
    }
    
    
  }
  
  if(verbose)
    cat("\n")
  
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
             method = method,
             twoway = DAGstructure@twoway)
  
  return(out)
}


#
# make weights flow up from leaf nodes upwards
#
update_weights <- function(parents, children, rejected, method)
{
  n <- length(rejected)
  queue <- rep(0, n)
  head <- tail <- 1
  
  children_left <- rep(NA,n)
  
  for(i in 1:n)
    children_left[i] <- length(children[[i]])        #initialize children_left will all children per node
  weights <- rep(0, n)
  
  # put not yet rejected leaf nodes on queue
  for(i in 1:n)
    if(!rejected[i] && children_left[i] == 0)
    {
      queue[tail] <- i
      tail <- tail + 1
      weights[i] <- 1
    }
  
  # as long as there are nodes on the queue, pass on weight if possible 
  while(head < tail)
  {
    current <- queue[head]
    head <- head + 1
    num_parents <- length(parents[[current]]) #needed for any-parent method
    
    num_unrejected <- 0
    for(i in parents[[current]])
      if(!rejected[i])
        num_unrejected <- num_unrejected + 1
    
    # give each unrejected parent its fair share 
    if(num_unrejected > 0)
    {
      for(i in parents[[current]])
        if(!rejected[i])
        {
          if(method == "all")
          {
            weights[i] <- weights[i] + weights[current] / num_unrejected
          }
          else if(method == "any")
          {
            weights[i] <- weights[i] + weights[current] / num_parents
          }
          
          children_left[i] <- children_left[i] - 1
          
          # put i on queue if it got something from all its children 
          if(children_left[i] == 0)
          {
            queue[tail] <- i
            tail <- tail + 1
          }
        }
      
      if(method == "all")
      {
        weights[current] <- 0
      }
      else if(method == "any")
      {
        #brackets needed to garantuee that that part becomes truly 1 when num_unrej=num_parents
        weights[current] <- weights[current] - weights[current] * (num_unrejected/ num_parents)
      }
      
    }
    
  }
  return(weights)
}


#function that makes leaf_based_sets. Each node has an index. Per node: give back the indices of its leaf nodes. Leaf node gets its own index.
set_in_leaves <- function(parents, children)
{
  n <- length(parents)
  queue <- rep(0, n)
  head <- tail <- 1
  
  children_left <- rep(NA,n)
  
  for(i in 1:n)
    children_left[i] <- length(children[[i]])
    
  leaf_based_sets <- vector("list",n)  
  
  # put leaf nodes on queue (it's all about their indices)
  for(i in 1:n)
    if(children_left[i] == 0)
    {
      queue[tail] <- i
      tail <- tail + 1
      leaf_based_sets[[i]] <- i
    }
  
  # as long as there are nodes on the queue, find parents 
  while(head < tail)
  {
    current <- queue[head]
    head <- head + 1
    
    # put current in leaf_based_sets of his parents (if any)
    for(i in parents[[current]])
    {
      leaf_based_sets[[i]] <- c(leaf_based_sets[[i]],leaf_based_sets[[current]])
      
      children_left[i] <- children_left[i] - 1
      
      # put i on queue if all its children are dealt with 
      if(children_left[i] == 0)
      {
        queue[tail] <- i
        tail <- tail + 1
        leaf_based_sets[[i]] <- unique(leaf_based_sets[[i]])
      }
    }
  }
  
  return(leaf_based_sets)
}

#TODO: merge with update_weights? 

#for LPnumdenomen, we need the maximal weight that a node can ever attain. This is the total incoming weight 
#from a weightflow where all leaf nodes participate in

update_max_weights <- function(parents, children, rejected)
{
  n <- length(rejected)
  queue <- rep(0, n)
  head <- tail <- 1
  
  children_left <- rep(NA,n)
  
  for(i in 1:n)
    children_left[i] <- length(children[[i]])        #initialize children_left will all children per node
  weights <- rep(0, n)
  max_weights <- rep(0,n) #needed for LPnumdenom
  
  # put unrejected leaf nodes on queue (this will be all leaf nodes, since the function is called before the algorithm starts)
  for(i in 1:n)
    if(!rejected[i] && children_left[i] == 0)
    {
      queue[tail] <- i
      tail <- tail + 1
      weights[i] <- 1
      max_weights[i] <- 1 #maximal weight of leaf nodes is just their initial weight
    }
  
  # as long as there are nodes: give weights to parent nodes
  while(head < tail)
  {
    current <- queue[head]
    head <- head + 1
    num_parents <- length(parents[[current]]) #NB LPnumdenom only needed if any parent method is chosen
    
    num_unrejected <- 0
    for(i in parents[[current]])
      if(!rejected[i])
        num_unrejected <- num_unrejected + 1
    
    # move the weight to the parents
    if(num_unrejected > 0)
    {
      for(i in parents[[current]])
        if(!rejected[i])
        {
          weights[i] <- weights[i] + weights[current] / num_parents
          
          children_left[i] <- children_left[i] - 1
          
          # put i on queue if all its children are processed
          if(children_left[i] == 0)
          {
            queue[tail] <- i
            tail <- tail + 1
            max_weights[i] <- weights[i] #maximal weight is weight just before you try to give the weight to your parents
          }
        }
      
      #brackets needed to ensure that last argument can become exactly one
      weights[current] <- weights[current] - weights[current] * (num_unrejected/ num_parents)
      
    }
    
  }
  return(max_weights)
}


#LP solver based on the package "gurobi"
solveLP_gurobi <-  function(col_index, row_index, val_index, obj, rhs, dir, relaxation, minmax)
{
  A <- sparseMatrix(i=row_index,j=col_index,x=val_index)
  
  model             <- list()  
  model$A           <- A
  model$obj         <- obj
  model$modelsense  <- minmax
  model$rhs         <- rhs
  model$sense       <- dir
  
  if(relaxation) #all parameters >= 0, <= 1 not needed to additionally specify since it's a minimization problem
    model$vtype <- 'S'
  else #all parameters binary
    model$vtype      <- 'B'    
  
  #needed by gurobi,  
  params <- list(OutputFlag=0)
  
  result <- gurobi::gurobi(model, params) #:: because gurobi function is not imported, only suggested, and not even truly loaded 
  result <- result$objval
  return(result)
  
}

#LP solver based on the package "lpSolve"
solveLP_lpsolve <- function(col_index, row_index, val_index, obj, rhs, dir, relaxation, minmax)
{
  A <- cbind(row_index,col_index,val_index)
  result <- lp(direction = minmax , objective.in = obj, const.dir = dir, const.rhs = rhs, all.bin=!relaxation, dense.const = A)
  result <- result$objval
  return(result)  
}


#function that chooses the LPsolver to work with, based on installed packages. gurobi is suggested, lpSolve is imported at installation
getLPsolver <- function() 
{
#  if(suppressWarnings(require(gurobi, quietly=TRUE)))
#    solveLP_gurobi
  if(suppressWarnings((requireNamespace("gurobi", quietly = TRUE))))
    solveLP_gurobi
  else
    solveLP_lpsolve
}


#LP solver that only optimizes the denominator 
LPdenom <- function(solveLP, implications, candidate, relaxation = TRUE)
{
  #find maximum index of leaf nodes 
  maxelem <- max(unlist(implications),unlist(candidate))
  
  n_impl <- length(implications)
  
  num_non_zero <- sum(sapply(implications,length))
  
  if(length(candidate) != 0)
    num_non_zero <- num_non_zero + length(candidate)
  
  row_index <- rep(NA,num_non_zero)
  col_index <- rep(NA,num_non_zero)
  val_index <- rep(NA,num_non_zero)
  
  counter <- 1
  
  #per implication: sum leaf nodes >= 1 (at least one rejection)
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
  
  #if necessary, add constraint candidate: sum leaf nodes in candidate = 0
  if(length(candidate) != 0)
  {
    for(j in candidate)
    {
      row_index[counter] <- n_impl + 1
      col_index[counter] <- j 
      val_index[counter] <- 1
      counter <- counter + 1      
    }
  }
  
  obj <- rep(1,maxelem)
  
  #minstens 1 verwerping in elke set: >= 1 
  if(length(candidate) != 0)
    rhs <- c(rep(1,n_impl), 0)
  else
    rhs <- rep(1, n_impl)
  
  #>= voor alles
  if(length(candidate) != 0)
    dir <- c(rep('>=', n_impl), '=')
  else
    dir <- rep('>=', n_impl)
  
  #choose ILP solver to work with and calculate minimal value objective function
  result <- solveLP(col_index=col_index, row_index=row_index, val_index=val_index, obj=obj, rhs=rhs, dir=dir, relaxation=relaxation, minmax="min")
  
  #no integer in case of relaxation, but has to be, so you can use ceiling. However rounding errorscan occur, even in case of ILP
  #use all.equal to check whether result is actually an integer, and using ceiling afterwards
  return(ifelse(isTRUE(all.equal(result, round(result))), result, ceiling(result)))
}





#optimize numerator & denominator for any-parent variant simultaneously. >= p_k/alpha*g2 means you can reject candidate
# g1*sum(parents) + p_k/alpha sum(leaves) - p_kalpha* #leaves >= 0 
#parentscand same format as implications: list
LPnumdenom <- function(solveLP, implications, candidate, parentscand, g1, palpha, relaxation = TRUE)
{
  maxelem <- max(unlist(implications),unlist(candidate),unlist(parentscand))
  
  n_impl <- length(implications)
  n_par <- length(parentscand)
  
  num_non_zero <- sum(sapply(implications,length))  
  if(n_par>0) #can happen that there are more root nodes, already an implication, but no parents
    num_non_zero <- num_non_zero + sum(sapply(parentscand,length))*2 #xj <= pi for all xj in pi 
  num_non_zero <- num_non_zero + length(candidate)
  
  row_index <- rep(NA,num_non_zero)
  col_index <- rep(NA,num_non_zero)
  val_index <- rep(NA,num_non_zero)
  
  counter <- 1
  
  #for each implication: sum of leaf nodes >= 1
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
  
  #add constraints parents: xj <= pi for all xj in pi (if node rejected, same holds for parent node): xj-pi <= 0
  if(n_par>0)
  {
    row <- n_impl + 1 #row index 
    for(i in 1:n_par)
    {
      for(j in parentscand[[i]])
      {
        row_index[counter] <- row
        col_index[counter] <- j
        val_index[counter] <- 1
        counter <- counter + 1
        row_index[counter] <- row
        col_index[counter] <- maxelem + i
        val_index[counter] <- -1
        counter <- counter + 1
        row <- row + 1
      }
    }
  }
  
  
  #add constraint candidate: sum leaves = 0
  for(j in candidate)
  {
    row_index[counter] <- n_impl + sum(sapply(parentscand,length)) + 1
    col_index[counter] <- j 
    val_index[counter] <- 1
    counter <- counter + 1      
  }
  
  
  #parents *g1, leaf nodes *pk/alpha
  obj <- c(rep(palpha,maxelem), rep(g1,n_par))
  #minstens 1 verwerping in elke set: >= 1, nothing in cand <=0, parents <=0
  rhs <- c(rep(1,n_impl),rep(0,sum(sapply(parentscand,length))+1))
  #>= voor impl, <= voor kandidaat en parents
  dir <- c(rep('>=', n_impl),rep('<=',sum(sapply(parentscand,length))+1))
  
  #choose ILP solver to work with and calculate minimal value objective function
  result <- solveLP(col_index=col_index, row_index=row_index, val_index=val_index, obj=obj, rhs=rhs, dir=dir, relaxation=relaxation, minmax="min")
  
  #does not have to be an integer, so no rouding in case of relaxation
  return(result)
}


#function that, given a number of sets, return the number of them that has to be rejected, given the rejections in the DAG object that is entered
#optimization can be LP or ILP
DAGpick <- function(DAG, indicators, optimization = "ILP")
{
  #check whether DAG has two-way property
  if(!DAG@twoway)
    stop("The DAG you work with has no two-way properties. At the moment it is not yet possible to use a DAG structure without these properties in this function.")
  
  implications <- DAG@leaf_based_sets[DAG@implications]
  sets <- DAG@leaf_based_sets[indicators] #selected sets in terms of leaf nodes
  unionsets <- Reduce(union, sets, NULL)
  implications <- Filter(function(i) all(i %in% unionsets), implications) #only look at implications that are fully embedded in the chosen sets
  
  maxelem <- max(unlist(implications))
  
  n_impl <- length(implications)
  n_sets <- length(sets)
  
  num_non_zero <- sum(sapply(implications,length)) + sum(sapply(sets,length))*2
  
  row_index <- rep(NA,num_non_zero)
  col_index <- rep(NA,num_non_zero)
  val_index <- rep(NA,num_non_zero)
  
  counter <- 1
  
  # contraint for every implication to be satisfied
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
  
  # constraints xj <= si for all j in Si
  row <- n_impl + 1 #row index
  for(i in 1:n_sets)
  {
    for(j in sets[[i]])
    {
      row_index[counter] <- row
      col_index[counter] <- j
      val_index[counter] <- 1
      counter <- counter + 1
      row_index[counter] <- row
      col_index[counter] <- maxelem + i
      val_index[counter] <- -1
      counter <- counter + 1
      row <- row + 1
    }
  }
  
  obj <- c(rep(0,maxelem), rep(1, n_sets)) #sum over chosen sets
  
  #mat least one rejection in eacht implication: >= 1 
  #for all gensets: xj - si <= 0
  rhs <- c(rep(1, n_impl), rep(0, sum(sapply(sets,length))))
  
  #>= voor alles
  dir <- c(rep('>', n_impl), rep('<', sum(sapply(sets,length))))
  
  #choose LPsolver to work with
  solveLP <- getLPsolver()

  #calculate minimal value objective function
  result <- solveLP(col_index=col_index, row_index=row_index, val_index=val_index, obj=obj, rhs=rhs, dir=dir, relaxation=(optimization == "LP"), minmax="min")
  
  #no integer in case of relaxation, but has to be, so you can use ceiling. However rounding errors can occur, even in case of ILP
  #use all.equal to check whether result is actually an integer, and using ceiling afterwards
  return(ifelse(isTRUE(all.equal(result, round(result))), result, ceiling(result)))
}


setMethod("show", "DAG", function(object) {
  meth <- switch(object@method, 
                 all = "all-parents",
                 any = "any-parent",
                 holm = "structured holm")
  
  cat("The ", meth, " method result on ", length(object@sets), " hypotheses.\n", sep="")
  cat("There are ", sum(object@rejected), " hypotheses rejected out of a total of ", length(object@sets), "\n", "at an alpha-level of ", object@alpha, ".\n", sep="")
})

setMethod("summary", "DAG", function(object) {
  meth <- switch(object@method, 
                 all = "all-parents",
                 any = "any-parent",
                 holm = "structured holm")
  
  cat("The ", meth, " method result on ", length(object@sets), " hypotheses.\n", sep="")
  cat("There are ", sum(object@rejected), " hypotheses rejected out of a total of ", length(object@sets), "\n", "at an alpha-level of ", object@alpha, ".\n", sep="")
})


setMethod("implications", "DAG", function(object) {
  #if no names in original sets, add index to this list
  allpvalues <- object@allpvalues
  if(length(names(allpvalues))==0)
    names(allpvalues) <- 1:length(allpvalues)
  impls <- allpvalues[object@implications]
  impls
})

setMethod("pvalue", "DAG", function(object, indicator) {
  object@allpvalues[indicator]
})

setMethod("alpha", "DAG", function(object) {
  object@alpha
})