kofnGA <- function(n,k,OF,popsize=200,keepbest=floor(popsize/10),ngen=500,
                   tourneysize=max(ceiling(popsize/10),2),mutprob=0.01,initpop=NULL, 
                   verbose=0, ...) {
# Genetic algorithm to do subset selection:  choose a subset of a fixed size k from the 
# integers 1:n, such that function OF() is minimized at that subset.  The selection step is
# done by tournament selection based on ranks.
#
# Inputs:
#   n           The number of objects to choose from.  The algorithm chooses a subset of integers
#               from 1 to n.
#   k           The number of objects to choose.
#   OF          The objective function.  The first argument of OF should be an index vector of k 
#               integers in the range [1, n]. 
#   popsize     The size of the population, the number of offspring produced each generation.
#   keepbest    This argument is used to implement elitism.  The keepbest least fit offspring each
#               generation are replaced by the keepbest most fit members of the previous generation.
#   ngen        The number of generations to run.
#   tourneysize The number of individuals to be selected for each tournament.
#   mutprob     The probability of mutation for each of the k chosen indices in each individual.
#               An index chosen for mutation simply jumps to any other unused index at random.
#   initpop     A popsize-by-k matrix of starting solutions.  Possibly useful if using this 
#               function in an adaptive, iterative, or parallel scheme.
#   verbose     An integer controlling the display of intermediate results during search. If 0 
#               (default), nothing is shown.  If verbose=c, the iteration number and best objective
#               value are displayed every c iterations.
#   ...         Used to pass other arguments to OF(subset,...) as necessary.
#
# Output is a list with the following elements:
#   bestsol     The best solution found.
#   bestobj     The objective function value for the best solution found.
#   pop         The final population, row-sorted in order of increasing objective function.
#   obj         The objective function values corresponding to each row of pop.
#   old         A list with elements giving for each generation: the best solution, the best
#               objective value, and the average objective value.
# 
# NOTES
# - Increasing tourneysize will put more "selection pressure," and will speed up convergence. 
#   Smaller tourneysize values will conversely promote better searching of the solution space.
# - Increasing the size of the elite group (keepbest) also promotes more homogeneity in the
#   population.

##=== Basic input checking =====================================================================##
stopifnot(n %% 1 == 0, n > 0, n >= k,
          k %% 1 == 0, k > 0, 
          popsize %% 1 == 0, popsize >= 2, popsize > keepbest, popsize > tourneysize,
          keepbest %% 1 == 0, keepbest >= 0,
          ngen %% 1 == 0, ngen >= 1,
          tourneysize %% 1 == 0, tourneysize >= 2,
          mutprob >= 0, mutprob <= 1,
          all(dim(initpop) == c(popsize,k)),
          verbose %% 1 == 0)
  
##=== Create needed objects ====================================================================##
indices <- 1:n                                   #-The set of possible objects to choose from.
if (keepbest>0) {
    elitespots <- 1:keepbest                     #-Portion of pop rows reserved for elites.
    newspots <- (keepbest+1):popsize             #-Portion of pop rows reserved for offspring.
}
fitness.old <- vector(mode="numeric",length=popsize)
fitness.new <- vector(mode="numeric",length=popsize)
offspring <- matrix(0,nrow=popsize,ncol=k)
old <- list()
old$best <- matrix(0,nrow=ngen+1,ncol=k)
old$obj <- vector(mode="numeric",length=ngen+1)
old$avg <- vector(mode="numeric",length=ngen+1)
Tourneys1 <- matrix(0,nrow=popsize,ncol=tourneysize)
Tourneys2 <- matrix(0,nrow=popsize,ncol=tourneysize)


##=== Initialize the population and evaluate their fitness =====================================##
# When no population supplied, initialize the population using randomly-chosen indices.
if (is.null(initpop))
    pop <- t(replicate(popsize,sample(indices,k)))
else
    pop <- initpop
for (i in 1:popsize) fitness.old[i] = OF(pop[i,], ...)
old$best[1,] = sort(pop[rank(fitness.old,ties.method="random")==1,])
old$obj[1] = min(fitness.old)
old$avg[1] = mean(fitness.old)


##=== Loop through generations =================================================================##
for (gen in 1:ngen) {

    # Choose mating pairs (selection)-------------------------------------------------------------
    # Matrix Tourneys1 is popsize-by-tourneysize.  Row i contains indices of the individuals 
    # chosen to be in the tournament to become parent 1 of the ith offspring.  Matrix Tourneys2 is 
    # created in the same way.  In each tournament, the probability of a competing individual 
    # being chosen is proportional to its fitness rank in the tournament (larger ranks given to 
    # more fit individuals).  This is done using pickfun().  The results, Parents1 and Parents2, 
    # are vectors of length popsize containing the row-indices of parents to take from pop.
    Tourneys1[,] <- t(replicate(popsize,sample(1:popsize,tourneysize)))
    Tourneys2[,] <- t(replicate(popsize,sample(1:popsize,tourneysize)))
    pickfun <- function(v) sample(v,1,prob=rank(-fitness.old[v]))
    Parents1 <- apply(Tourneys1,1,pickfun)
    Parents2 <- apply(Tourneys2,1,pickfun)
    
    # Create offspring (crossover)----------------------------------------------------------------
    # For crossover, just combine the unique elements of both parents, then keep k of them, chosen
    # at random.
    for (i in 1:popsize) {
        combo <- unique(c(pop[Parents1[i],],pop[Parents2[i],]))
        offspring[i,] <- sample(combo,k)
    }    
    
    # Perform mutation----------------------------------------------------------------------------
    # Logical matrix chosen has the same dimensions as pop, with TRUE elements indicating 
    # locations chosen for mutation.  Then each offspring in turn is mutated by replacing its
    # chosen elements with a randomly-selected unused index.
    chosen <- matrix(as.logical(rbinom(popsize*k,1,mutprob)),nrow=popsize)
    nchosen <- apply(chosen,1,sum)
    for (i in 1:popsize) {
        if (nchosen[i]>0) {
            toadd <- sample(indices[-offspring[i,]],nchosen[i])
            offspring[i,chosen[i,]] = toadd
        }
    }
    
    # Evaluate fitness----------------------------------------------------------------------------
    # Put the fitness values for the offspring in a separate vector fitness.new.  The next 
    # generation will be chosen based on both fitness.old and fitness.new, since the keepbest
    # elites from the previous generation will stay in.
    for (i in 1:popsize) fitness.new[i] = OF(offspring[i,], ...)
    
    # Form the next generation, keeping the elites------------------------------------------------
    # If keepbest=0, then no elitism.  Just set pop and fitness.old to the offspring values.
    # If keepbest>0, use logical vectors old.keep and new.keep to identify the keepbest most fit  
    # members of the previous generation and the (popsize-keepbest) most fit members of the 
    # offspring.  Then overwrite pop and fitness.old appropriately.
    if (keepbest==0) {
        pop <- offspring
        fitness.old <- fitness.new   
    }
    else {
        old.keep <- rank(fitness.old,ties.method="random") <= keepbest
        new.keep <- rank(fitness.new,ties.method="random") <= popsize-keepbest
        pop[elitespots,] = pop[old.keep,]
        pop[newspots,] = offspring[new.keep,]
        fitness.old[elitespots] = fitness.old[old.keep]
        fitness.old[newspots] = fitness.new[new.keep]   
    }
    
    # Record progress of the search---------------------------------------------------------------
    old$best[gen+1,] <- sort(pop[rank(fitness.old,ties.method="random")==1,])
    old$obj[gen+1] <- min(fitness.old)
    old$avg[gen+1] <- mean(fitness.old)
    
    # Give output to console, if requested--------------------------------------------------------
    if (verbose>0 && gen%%verbose==0) {
      cat("Finished iteration ", gen, ". Best OF value = ", old$obj[gen+1], "\n")
    }
    
}

##=== Package the outputs ======================================================================##
# Note, when returning the overall best, need to look at the history, since elitism might be 
# turned off (then final population might not contain best solution found).
out <- list()
# Return the record of all runs.
out$old <- old
# Return the final population, sorted.
ord <- order(fitness.old)
out$pop <- matrix(0,nrow=popsize,ncol=k)
for (i in 1:popsize) {
    out$pop[i,] = sort(pop[ord[i],])
}
# Return the objfun values for the final population.
out$obj <- fitness.old[ord]
# Return the best solution found, and its objfun value.
alltimebest <- which(old$obj==min(old$obj,na.rm=TRUE))
alltimebest <- tail(alltimebest,1)               #-In case multiple bests, take last one.
out$bestsol <- out$old$best[alltimebest,]
out$bestobj <- out$old$obj[alltimebest]

class(out) <- "GAsearch"
out

}
