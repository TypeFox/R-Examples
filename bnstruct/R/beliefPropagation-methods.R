#' @rdname belief.propagation
#' @aliases belief.propagation,InferenceEngine
setMethod("belief.propagation",
          c("InferenceEngine"),
          function(ie, observations = NULL, return.potentials = FALSE){
            {
              ###############################
              # moved inside in order to eliminate a NOTE in R CMD check
              proc.order <- function(node, from, adj)
              {
                # Recursive method to compute order of the message passing in the upward step.
                #
                # node : current node
                # from : (local) root
                # adj  : adjacency matrix
                neighbours <- setdiff(which(adj[node,] > 0), from)
                
                if (length(neighbours) > 0)
                {
                  for (n in neighbours) {
                    proc.order(n, node, adj)
                    parents.list <<- c(parents.list, node)
                  }
                }
                
                process.order <<- c(process.order, node)
              }
              
              ##############################
              
              # if (missing(net))
              # {
              net <- bn(ie)
              # }
              
              if (missing(observations))
              {
                obs <- observations(ie)
                observed.vars <- obs$observed.vars
                observed.vals <- obs$observed.vals
              }
              else
              {
                observed.vars <- observations[[1]]
                observed.vals <- observations[[2]]
                obs <- unique.observations(observed.vars, observed.vals)
                observed.vars <- obs$observed.vars
                observed.vals <- obs$observed.vals
                observations(ie) <- list(observed.vars, observed.vals)
              }
              
              num.nodes  <- num.nodes(net)
              num.cliqs  <- num.nodes(ie)
              
              # cliques contains the variables that compose each clique
              cliques    <- jt.cliques(ie)
              
              ctree      <- junction.tree(ie)
              
              cpts       <- cpts(net)
              
              variables  <- variables(net)
              
              dim.vars   <- lapply(1:num.nodes,
                                   function(x)
                                     #as.list(
                                       match(
                                         c(unlist(
                                          names(dimnames(cpts[[x]])), F, F
                                         )),
                                         variables
                                       )
                                   #)
                            )

              
              node.sizes <- node.sizes(net)
              
              # potentials is a list containing the probability tables of each clique
              potentials <- as.list(rep(as.list(c(1)), num.cliqs))
              
              # dimensions.contained contains the variables that effectively compose the cpt
              # currently contained in each node of the clique tree.
              # After last round it will match corresponding clique.
              dimensions.contained <- lapply(1:num.cliqs, function(x) as.list(c(NULL)))
              
              
              # choose as root (one among) the clique(s) whose connected edges have the highest overall sum
              root <- which.max(rowSums(ctree))
              
              # Assign factors to a cluster graph
              # Construct initial potentials:
              # initial potentials are conditional or joint probability tables, depending on the initial BN
              # and how we assign each CPT to the cliques.
              #
              # The probabilities are stored as multidimensional arrays, with the convention that
              # the variables in the tables are in alphanumerical order.
              # E.g.: the network A -> C <- D -> B, whose junction tree has two nodes (and tables) ACD and BD.
              #
              # We construct a table for a clique this way:
              # - start with an empty table (NULL)
              # - whenever a {C,J}PT is assigned to that clique, its variables are ordered
              # - then, we control if the variables of the table we're inserting are already present in the clique table:
              #   - if no, we can multiply the two tables, ensuring the variables are ordered after that
              #   - if yes, we multiply the two tables
              # If the clique is currently empty, just add the cpt.
              # We have, however, to maintain the order of the variables in the probability table.
              
              for (cpt in 1:num.nodes)
              {
                # find a suitable clique for the CPT of node `cpt` (a clique that contains node `cpt` and all of its parents)
#                 target.clique <- which.min(lapply(1:num.cliqs,
#                                                   function(x){
#                                                     length(
#                                                       which(unlist(
#                                                         is.element(
#                                                           c(unlist(dim.vars[[cpt]])),
#                                                           c(cliques[[x]])
#                                                         )
#                                                       ) == FALSE) == 0)
#                                                   }
#                 ))
                
                # TODO: PROFILING: is there anything better?
                target.clique <- which.min(lapply(1:num.cliqs,
                                                  function(x){
                                                    length(which(!is.na(match(
                                                      c(unlist(dim.vars[[cpt]])),
                                                      c(cliques[[x]])
                                                    )) == FALSE))
                                                  }))
                
                # get the variables currently contained in the selected clique
                ds <- unlist(dimensions.contained[[target.clique]], F, F)
                
                if (length(ds) == 0)
                {
                  # if current clique is empty, just insert the cpt
                  out <- sort.dimensions(cpts[[cpt]], dim.vars[[cpt]])
                  potentials[[target.clique]]           <- out$potential
                  dimensions.contained[[target.clique]] <- out$vars
                }
                else
                {
                  # multiply current prob. table for the already inserted prob. table
                  out <- mult(potentials[[target.clique]],
                              dimensions.contained[[target.clique]],
                              cpts[[cpt]],
                              dim.vars[[cpt]],
                              node.sizes)
                  potentials[[target.clique]]           <- out$potential
                  dimensions.contained[[target.clique]] <- out$vars
                }

              }
              
              # INCORPORATE EVIDENCE
              # If there are any observed variables, insert the knowledge.
              # Each observation is inserted by setting to zero all of the combinations that
              # do not match the observation. This is done this way:
              # - find a clique that contains the observed variable
              # - create a probability table for the observed variable, with only one non-zero entry,
              #   the one corresponding to the observed value
              # - multiply the table of the clique for the newly created table
              # - normalize after belief propagation
              
              if (length(observed.vars) > 0)
              {
                # observed.vars <- c(unlist(observed.vars, F, F))
                
                if (class(observed.vars) == "character") # hope that the user gives coherent input...
                  observed.vars <- sapply(observed.vars, function(x) which(x == variables))
                
                for (var in 1:length(observed.vars))
                {
                  if (observed.vals[var] <= 0 || observed.vals[var] > node.sizes[observed.vars[var]])
                  {
                    message(cat("Variable", observed.vars[var], "cannot take value", observed.vals[var], ", skipping..."))
                  }
                  else
                  {
                    # look for one clique containing the variable
                    target.clique <- which.min(lapply(1:num.cliqs,
                                                      function(x) {
                                                        which(is.element(
                                                          unlist(dimensions.contained[[x]]),
                                                          observed.vars[var]
                                                        ) == TRUE)
                                                      }
                    ))

                    tmp <- rep(0, node.sizes[observed.vars[var]])
                    tmp[observed.vals[var]] <- 1
                    out <- mult(potentials[[target.clique]],
                                dimensions.contained[[target.clique]],
                                tmp,
                                observed.vars[var],
                                node.sizes)
                    potentials[[target.clique]]           <- out$potential
                    dimensions.contained[[target.clique]] <- out$vars
                  }
                }
              }
              
              
              # compute processing order from leaves to root
              process.order <- c()
              parents.list  <- c()
              proc.order(root, c(), ctree)              
              
              # MESSAGE PASSING FROM LEAVES TO ROOT
              
              # msg.pots contains the prob.tables for the messages,
              # while msg.vars contains the corresponding variables
              msg.pots <- NULL
              msg.vars <- NULL
              for (i in 1:num.cliqs)
              {
                msg.pots[[i]] <- as.list(c(NULL))
                msg.vars[[i]] <- as.list(c(NULL))
              }
              
              # For each clique (excluding the root) compute the message by marginalizing
              # the variables not in the separator, then store the message and multiply it
              # for the cpt contained in the neighbour clique, overwriting the corresponding
              # potential and associated variables.
              for (clique in 1:(length(process.order)-1))
              {
                out <- compute.message(potentials[[process.order[clique]]],
                                       dimensions.contained[[process.order[clique]]],
                                       cliques[[process.order[clique]]],
                                       cliques[[parents.list[clique]]],
                                       node.sizes)
                msg.pots[[process.order[clique]]] <- out$potential
                msg.vars[[process.order[clique]]] <- out$vars
                
                bk <- potentials[[parents.list[clique]]]
                bkd <- dimensions.contained[[parents.list[clique]]]
                out <- mult(potentials[[parents.list[clique]]],
                            dimensions.contained[[parents.list[clique]]],
                            msg.pots[[process.order[clique]]],
                            msg.vars[[process.order[clique]]],
                            node.sizes)
                potentials[[parents.list[clique]]]           <- out$potential
                dimensions.contained[[parents.list[clique]]] <- out$vars
                
              }
                            
              # Upward step is thus completed. Now go backward from root to leaves.
              # This step is done by taking the CPT of the root node and dividing it (for each child)
              # by the message received from the corresponding child, then marginalize the variables
              # not in the separator and pass the new nessage to the child. As the messages computed
              # in the upward step are not needed anymore after the division, they can be overwritten.
              # Then multiply the cpt of the child for the message computed, and iterate by treating
              # each (internal) node as root.
              for (clique in (length(process.order)-1):1)
              {                
                out <- divide(potentials[[parents.list[clique]]],
                              dimensions.contained[[parents.list[clique]]],
                              msg.pots[[process.order[clique]]],
                              msg.vars[[process.order[clique]]],
                              node.sizes)
                msg.pots[[process.order[clique]]] <- out$potential
                msg.vars[[process.order[clique]]] <- out$vars
                
                out  <- compute.message(msg.pots[[process.order[clique]]],
                                        msg.vars[[process.order[clique]]],
                                        cliques[[parents.list[clique]]],
                                        cliques[[process.order[clique]]],
                                        node.sizes
                )
                msg.pots[[process.order[clique]]] <- out$potential
                msg.vars[[process.order[clique]]] <- out$vars
                
                out <- mult(potentials[[process.order[clique]]],
                            dimensions.contained[[process.order[clique]]],
                            msg.pots[[process.order[clique]]],
                            msg.vars[[process.order[clique]]],
                            node.sizes)
                potentials[[process.order[clique]]] <- out$potential
                dimensions.contained[[process.order[clique]]] <- out$vars
              }
              
              # Finally, normalize and add dimension names and return the potentials computed (will be all JPTs).
              for (x in 1:num.cliqs) {
                s <- sum(potentials[[x]])
                potentials[[x]] <- potentials[[x]] / s
                dmns <- list(NULL)
                for (i in length(dimensions.contained[[x]]))
                {
                  dmns[[i]] <- c(1:node.sizes[dimensions.contained[[x]][[i]]])
                }
                dimnames(potentials[[x]])        <- dmns
                names(dimnames(potentials[[x]])) <- as.list(variables[c(unlist(dimensions.contained[[x]]))])
              }
              
              if (return.potentials)
                return(potentials)
              
              ###################
              # Now create new BN with updated beliefs
              ###################
              #
              # To buil new conditional probability tables for the original network,
              # starting from the updated beliefs, we proceed this way:
              # for each node of the BN:
              # - if it is an observed variable, construct a prob.table containing only
              #   one variable (the one of the node) with one only non-zero (in fact, 1)
              #   entry, the one corresponding to the observed value
              # - if it is a non-observed variable, find a clique that contains the variables
              #   of the original cpt (it must exist, because of the moralization - we have
              #   already used this property), sum out the possible other variables introduced
              #   by the triangulation, and divide the JPT by the JPT of the parent nodes
              #   (e.g.: if we start from P(ABCD) and want P(C|A,B), we sum out D, then we
              #   obtain P(AB) by summing out C, and then we compute P(ABC)/P(AB) = P(C|A,B)).
              #   Easy peasy.
              
              
              
              nbn <- BN()
              name(nbn)         <- name(net)
              num.nodes(nbn)    <- num.nodes(net)
              variables(nbn)    <- variables(net)
              node.sizes(nbn)   <- node.sizes(net)
              discreteness(nbn) <- discreteness(net)
              dag(nbn)          <- dag(net)
              wpdag(nbn)        <- wpdag(net)
              scoring.func(nbn) <- scoring.func(net)
              struct.algo(nbn)  <- struct.algo(net)

              ncpts <- NULL # lapply(1:num.nodes, function(x) as.list(c(NULL)))
              
              for (node in 1:num.nodes)
              {
                # faster, but result does not change. While debugging, better keep this out...
#                 mpos <- match(node, observed.vars)
#                 if (!is.na(mpos))
#                 # works also when observed.vars == c()
#                 # in that case, the `else` branch will be the chosen one for every variable
#                 {
#                   ncpts[[node]] <- array(rep(0, node.sizes[observed.vars[mpos]]),
#                                          c(node.sizes[observed.vars[mpos]]))
#                   ncpts[[node]][observed.vals[mpos]] <- 1
#                   dimnames(ncpts[[node]]) <- list(c(1:node.sizes[observed.vars[mpos]]))
#                   names(dimnames(ncpts[[node]])) <- as.list(variables[observed.vars[mpos]])
#                 }
#                 else
                {
                  dnode <- unlist(dim.vars[[node]], F, F)
#                   target.clique <- which.min(lapply(1:num.cliqs,
#                                                     function(x){
#                                                       length(
#                                                         which(c(
#                                                           is.element(
#                                                             dnode,
#                                                             unlist(dimensions.contained[[x]])
#                                                           )
#                                                         ) == FALSE) == 0)
#                                                     }
#                   ))
                  
                  target.clique <- which.min(lapply(1:num.cliqs,
                                            function(x){
                                              length(which(!is.na(match(
                                                dnode,
                                                unlist(dimensions.contained[[x]])
                                              )) == FALSE))
                                            }))
                  
                  pot <- potentials[[target.clique]]
                  dms <- c(unlist(dimensions.contained[[target.clique]]))
                  vs  <- c(unlist(dim.vars[[node]]))
                  
#                   others <- setdiff(dms,vs)
#                   for (var in others)
#                   {
#                     out <- marginalize(pot, dms, var)
#                     pot <- out$potential
#                     dms <- out$vars
#                   }
                  remaining <- match(vs, dms)
                  dms <- dms[remaining]
                  pot <- apply(pot, remaining, sum)
                  pot <- pot / sum(pot)
                  
#                   cat(node, " ", dms,"\n")
#                   print(pot)
                  
                  if (length(dms) > 1)
                  {
                    pot.bak <- pot
                    dms.bak <- dms
                    # readLines(file("stdin"),1)
#                     out <- marginalize(pot, dms, node)
#                     pot <- out$potential
#                     dms <- out$vars
                    remaining <- (1:length(dms))[-which(dms == node)]
                    dms <- dms[remaining]
                    pot <- apply(pot, remaining, sum)
                    pot <- pot / sum(pot)
                    out <- divide(pot.bak,
                                  dms.bak,
                                  as.array(pot),
                                  dms, # out$vars,
                                  node.sizes)
                    
                    pot <- out$potential
                    #pot <- pot / sum(pot)
                    dms <- out$vars
                  }
                  pot <- as.array(pot)
                  dmns <- list(NULL)
                  for (i in length(dms))
                  {
                    dmns[[i]] <- c(1:node.sizes[dms[i]])
                  }
                  dimnames(pot)        <- dmns
                  names(dimnames(pot)) <- as.list(variables[dms])
                  ncpts[[node]] <- pot
#                   print(pot)
#                   readLines(file("stdin"),1)
                }
                
              }
              
              cpts(nbn) <- ncpts
              
              updated.bn(ie) <- nbn
              jpts(ie) <- potentials
              return(ie)
            }
          })




compute.message <- function(pot, dp, vfrom, vto, node.sizes)
{
  # Compute message from one node to another:
  # marginalize variables not in separator between two nodes.
  #
  # pot   : cpt to be marginalized
  # dp    : dimensions of the potential (may not contain all of the variables
  #         that have to be present in the clique, if this is performed in
  #         the upward step)
  # vfrom : variables in the sending clique
  # vto   : variables in the receiving clique
  # node.sizes : node sizes
  
  # separator is made of the shared variables between the two cliques
  vars.msg <- unlist(vfrom, F, F)
  sep      <- intersect(vars.msg, unlist(vto, F, F))
  dp       <- unlist(dp, F, F)
  
  # for all of the variables not in the separator, repeat marginalization
  # shrinking the prob.table
#   for (var in setdiff(vars.msg, sep))
#   {
#     if (length(intersect(var, dp)))
#     {
#       msg  <- marginalize(pot, dp, var)
#       pot  <- msg$potential
#       dp   <- msg$vars
#     }
#   }
  
  # othervars <- setdiff(dp, sep)
  remaining <- match(intersect(dp,sep), dp)
  if (length(remaining) > 0)
  {
    dp        <- dp[remaining] # rem.vars <- ...
    pot <- apply(pot, remaining, sum)
  }
  
#   if (sum(out) != sum(pot))
#   {
#     print(sum(out))
#     print(sum(pot))
#     readLines(file("stdin"),1)
#   }
  
  return(list("potential"=pot, "vars"=dp))
}


marginalize <- function(pot, vars, marg.var)
{
  # Marginalize a variable in a probability table.
  #
  # pot      : probability table
  # vars     : variables associated to pot
  # marg.var : variable to be marginalizes
  
  marg.dim <- which(unlist(vars, F, F) == marg.var)
#   print("marg.dim")
#   print(marg.dim)
  
  # get dimensions, compute dimensions for the soon-to-be-created prob. table
  # and number of the values that it will contain
  dims          <- dim(pot)
  new.dims      <- dims[-marg.dim]
  new.num.vals  <- prod(new.dims)
  length.of.run <- dims[marg.dim]
  num.runs      <- new.num.vals
  
  # switch dimensions in the array:
  # dimension corresponding to the variable to be marginalized goes first
  new.order <- c(marg.dim, (1:length(dims))[-marg.dim])
  
  # remove marginalized dimension name
  new.order.names <- vars[-marg.dim]
  
  # switch dimensions, make prob.table  a linear array,
  cpt  <- aperm(pot, new.order)
  # marg <- c(cpt)
  
  # the marginalization is now done by summing consecutive values
  marg <- tapply(c(cpt), rep(1:new.num.vals, each=length.of.run), sum)
  # marg <- rowSums(matrix(cpt, new.num.vals, length.of.run))
  # marg <- aggregate(as.data.frame(marg), list(rep(1:new.num.vals, each=length.of.run)), sum)[,2]
  # tapply is much faster than the conversion into a data.frame
  
  # apply new dimensions to resulting list (if needed)
  if (length(new.order.names) > 0)
  {
    marg <- array(marg, new.dims)
  }
  
  return(list("potential"=marg, "vars"=new.order.names))
}


mult <- function(cpt1, vars1, cpt2, vars2, node.sizes)
{
  # Multiply a cpt by another cpt.
  # Returns a list containing the resulting cpt and the associated variables list.
  #
  # cpt1  : first cpt
  # vars1 : variables associated to cpt1
  # cpt2  : second cpt
  # vars2 : variables associated to cpt2
  # node.sizes : sizes of the nodes
  
  # clean format
  vars1 <- unlist(vars1, F, F)
  vars2 <- unlist(vars2, F, F)
  
  # If the variables associated to cpt1 are all contained in the list of variables for cpt2, but
  # no variables of cpt2 is contained also in cpt1, swap the two cpts.
  # Handles cases such as P(AB) x P(C|AB) ==> P(C|AB) x P(AB)
  # Not really needed, but easier to understand.
  if ((length(setdiff(vars1, vars2)) == 0 &&
       length(setdiff(vars2, vars1)) > 0     )        
  )
  {
    tmp   <- vars1
    vars1 <- vars2
    vars2 <- tmp
    
    tmp  <- cpt1
    cpt1 <- cpt2
    cpt2 <- tmp
  }
  
  # For (my) simplicity, cpts are managed with the variables (and therefore dimensions) in ascending order.
  # Check this requirement, and take action if it is not met.
  out   <- sort.dimensions(cpt1, vars1)
  cpt1  <- out$potential
  vars1 <- out$vars
  
  out   <- sort.dimensions(cpt2, vars2)
  cpt2  <- out$potential
  vars2 <- out$vars
  
  # Proper multiplication starts here.
  # It works like this:
  # - look for the common variables in vars1 and vars2;
  common.vars <- intersect(vars1, vars2)
  #common1 <- match  (vars1, common.vars)
  #common1 <- which(!is.na(common1), TRUE)
  common1 <- match(common.vars, vars1)
  #common2 <- match  (vars2, common.vars)
  #common2 <- which(!is.na(common2), TRUE)
  common2 <- match(common.vars, vars2)

  
  # - if the cpts share no common variables, we can multiply them with an outer product;
  if (length(common.vars) == 0)
  {
    cpt1 <- as.vector(cpt1) %o% as.vector(cpt2)
  }
  else
    # otherwise, we have to manage the shared variables: consider P(C|A) x P(AB); unlisting the cpts we obtain
    # [ac !ac a!c !a!c], and
    # [ab !ab a!b !a!b]
    # (remember we have ordered the dimensions in ascending order). We have to handle the shared A and compute
    # ac   x ab   = acb
    # ac   x a!b  = ac!b
    # !ac  x !ab  = !acb
    # !ac  x !a!b = !ac!b
    # a!c  x ab   = a!cb
    # a!c  x a!b  = a!c!b
    # !a!c x !ab  = !a!cb
  # !a!c x !a!b = !a!c!b
  # (we will then have to reorder dimensions).
  # In order to do this, we permute dimensions for the two cpts, by putting the shared variables
  # as first dimensions of cpt1 and last dimensions of cpt2; then, we consecutively repeat every cell
  # of cpt1, and the entire sequence of cells of cpt2, the proper number of times in order
  # to reach the final number of elements (the ``proper number'' is the product of the size of
  # the variables not shared among the two cpt2); then we can finally compute the
  # element-wise product of cpt1 and cpt2;
  {
    if (length(vars1) > 1)
    {
      new.order <- c((1:length(vars1))[common1],
                     (1:length(vars1))[-common1])
      cpt1      <- aperm(cpt1, new.order)
      vars1     <- vars1[new.order]
      # common1 <- match  (vars1, common.vars)
      # common1 <- which(!is.na(common1), TRUE)
      common1 <- match(common.vars, vars1)
    }
    
    # [a b c] ==> [a a b b c c]
    if (length(vars2[-common2]) > 0)
    {
      cpt1  <- c(sapply(c(cpt1),
                        function(x){
                          rep(x, 
                              prod(node.sizes[vars2[-common2]])
                          )
                        }
      ))
    }
    
    if(length(vars2) > 1)
    {
      new.order <- c((1:length(vars2))[-common2],
                     (1:length(vars2))[common2])
      cpt2      <- aperm(cpt2, new.order)
      vars2     <- vars2[new.order]
      # common2 <- match  (vars2, common.vars)
      # common2 <- which(!is.na(common2), TRUE)
      common2 <- match(common.vars, vars2)
    }
    
    # [a b c] ==> [a b c a b c]
    cpt2  <- rep(c(cpt2),
                   prod(node.sizes[vars1[-common1]])
    )
    
    # - point-wise product
    cpt1 <- c(cpt1) * c(cpt2)
    
  }
  
  #vars1 <- unlist(vars1, F, F)
  #vars2 <- unlist(vars2, F, F)
  
  # - compute variables for the resulting cpt; if there were no shared variables, then
  #   it suffices to concatenate 
  if (length(common.vars) > 0)
  {
    # TODO: PROFILING: anything better than this?
    new.where <- which(is.na(match(vars2, common.vars)) == FALSE) # ugly, but should be more robust
    vars1     <- c(vars2[-new.where], vars1)
  }
  else
  {
    vars1 <- c(vars1, vars2)
  }
  
  cpt1 <- array(c(cpt1), c(node.sizes[vars1]))

  out  <- sort.dimensions(cpt1, vars1)

 return(list("potential"=out$potential, "vars"=out$vars))
}


divide <- function(cpt1, vars1, cpt2, vars2, node.sizes)
{
  # Divide a cpt by another cpt.
  # Returns a list containing the resulting cpt and the associated variables list.
  # cpt1  : dividend cpt
  # vars1 : variables associated to cpt1
  # cpt2  : divisor cpt
  # vars2 : variables associated to cpt2
  # node.sizes : sizes of the nodes
  
  # clean format
  vars1 <- unlist(vars1, F, F)
  vars2 <- unlist(vars2, F, F)
  
  # If the variables associated to cpt1 are all contained in the list of variables for cpt2, but
  # no variables of cpt2 is contained also in cpt1, swap the two cpts.
  # Handles cases such as P(AB) / P(C|AB) ==> P(C|AB) / P(AB)
  if ((length(setdiff(vars1, vars2)) == 0 &&
         length(setdiff(vars2, vars1)) > 0     )        
  )
  {
    tmp   <- vars1
    vars1 <- vars2
    vars2 <- tmp
    
    tmp  <- cpt1
    cpt1 <- cpt2
    cpt2 <- tmp
  }
#   print("--")
#   print(vars1)
#   print(vars2)
  
  # For (my) simplicity, cpts are managed with the variables (and therefore dimensions) in ascending order.
  # Check this requirement, and take action if it is not met.
  out   <- sort.dimensions(cpt1, vars1)
  cpt1  <- out$potential
  vars1 <- out$vars
  
  out   <- sort.dimensions(cpt2, vars2)
  cpt2  <- out$potential
  vars2 <- out$vars
  
  
  # The proper division starts here.
  # It works like this:
  # - domain of the divisor is entirely contained into the one of the dividend;
  # - look for the common variables (all of the variables in vars2, some of them in vars1);
  common.vars <- intersect(vars1, vars2)
  if (length(common.vars) == 0)
    return(list("potential"=cpt1, "vars"=vars1))
  # common1 <- match(vars1, common.vars)
  # common1 <- which(!is.na(common1), TRUE)# common1[!is.na(common1)]
  common1 <- match(common.vars, vars1)
  
  # - permute array dimensions for cpt1 putting the common variables in the first dimensions;
  if (length(vars1) > 1)
  {
    cpt1 <- aperm(cpt1, c((1:length(vars1))[common1],
                          (1:length(vars1))[-common1]
    ))
    vars1 <- c(vars1[(1:length(vars1))[common1]],
               vars1[(1:length(vars1))[-common1]])
    # common1 <- match  (vars1, common.vars)
    # common1 <- which(!is.na(common1), TRUE) # common1[!is.na(common1)]
    common1 <- match(common.vars, vars1)
  }
  
  # - unlist cpt2 and repeat it as many times as needed (product of cardinality
  #   of non-common variables of cpt1);
  cpt2      <- rep(c(cpt2),
                     prod(
                       node.sizes[vars1[-common1]]
                     )
  )
  
  # - now, every cell of cpt1 is paired with a cell of cpt2 whose variables
  #   have the same setting of it;
  # - perform element-wise division, handling the 0/0 case (conventionally set to 0 too);
  cpt1 <- sapply(1:length(cpt2),
                 function(x) {
                   if(cpt2[x] == 0) {
                     return(0)
                   } else {
                     return(cpt1[x] / cpt2[x])
                   }
                 })
  
  # - rebuild array with corresponding dimensions, and permute dimensions to reconstruct order
  #print(cpt1)
  cpt1 <- array(c(cpt1), node.sizes[vars1])
  out  <- sort.dimensions(cpt1, vars1)
  
  return(list("potential"=out$potential, "vars"=out$vars))
}

sort.dimensions <- function(cpt, new.ordering)
{
  # Permute array dimensions of the cpt accoring to the dimension names.
  # Dimensions (each corresponding to a variable) will be sorted in numerical order.
  # cpt  : conditional probability table
  # vars : dimension names
  
  # new.ordering <- unlist(vars, F, F)
  if (length(new.ordering) > 1)
  {
    # look if there is some variable out of order (preceding a variable with a lower number)
    is.ordered <- TRUE
    for (i in 1:(length(new.ordering)-1))
    {
      if (new.ordering[i] >= new.ordering[i+1])
      {
        is.ordered <- FALSE
        break
      }
    }
    
    # is.ordered <- is.na(match(FALSE,sapply(1:(length(new.ordering)-1), function(x) new.ordering[x] >= new.ordering[x+1])))
    
    if (!is.ordered) # permute
    {
      dd           <- data.frame(c1 = new.ordering,
                                 c2 = 1:length(new.ordering))
      dd           <- dd[with(dd,order(c1)),]
      new.ordering <- new.ordering[dd[,"c2"]]
      cpt          <- aperm(cpt, c(dd[,"c2"]))
    }
  }
  return(list("potential"=cpt, "vars"=new.ordering))
}
