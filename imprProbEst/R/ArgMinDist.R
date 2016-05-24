`ArgMinDist` <-
function(x,lbomega,ubomega,epsilon,ImpreciseModel){
  # x = the observations
  # it is assumed that the sample space Omega is a cuboid the k-dim.
  #   real vector space.
  #   lbomega = the vector which containes the lower bounds of the
  #      cuboid Omega.
  #   luomega = the vector which containes the upper bounds of the
  #      cuboid Omega.
  #   Hence, lbomega and luomega have k components.
  # epsilon = distance between two nodes

  PreResult <- 10  # This is a cautiuous pre-initialization because the total variation distance
                   # should never be greater than 1.
  theta <- 0  # This line is needless because the initialization of PreResult guarantees that
              # theta is defined in the first if-case of the j-loop
  k <- length(lbomega)
  x <- matrix(x,nrow=k)

  # The number of elements in Theta:
  t <- length(ImpreciseModel)

  # The number of observations:
  n <- length(x[1,])

  # The supporting lattice nodes
  omega <- BuildSupportingNodes(lbomega,ubomega,epsilon)   # the supporting lattice nodes
                                                           # in addition to the observations


  # By heuristic arguments, the following procedure makes a guess
  # which ordering of Theta might be a good one for the
  # calculation of the total variation distances.

  TotalVarIndicator <- function(FunctionsAndUpperPrevisions,observations){
       ListOfFunctions <- FunctionsAndUpperPrevisions[[1]]
       UpperPrevisions <- FunctionsAndUpperPrevisions[[2]]
       x <- observations

       EmpiricalValuesOfF <- fevaluation(x,ListOfFunctions, useApply = FALSE)
       EmpiricalIntegrals <- EmpiricalValuesOfF%*%rep(1/n,n)
       max(EmpiricalIntegrals-UpperPrevisions)
  }
  ArrayTotalVarIndicators <- matrix(unlist(lapply(ImpreciseModel, TotalVarIndicator, observations=x)), ncol=1)
  GuessOfOrdering <- sort.list(ArrayTotalVarIndicators)




  # For every parameter theta=j, the total variation distance is calculated in the following j-loop:

  runs <- 0
  for( jj in 1:t ){
       j <- GuessOfOrdering[jj]  # This line guarantees that the j-loop runs in that ordering
                                 # wich has been guessed
       if (PreResult > ArrayTotalVarIndicators[j]-0.002){ # If this was not true, there would definitely be no chance to get
                                                          # better by applying the procedure for this parameter j.
                                                          # So, if this was not true, we do not have to do anything 
                                                          # for this parameter j.
           runs <- runs+1
       
           ListOfFunctions <- ImpreciseModel[[j]][[1]]
           UpperPrevisions <- ImpreciseModel[[j]][[2]]
           s <- length(ListOfFunctions)     # Number of functions f

           # Calculate the discretized prevision
             DiscretePrevision <- DiscretizeImpreciseModel(ListOfFunctions,UpperPrevisions,x,omega,s)
       
       
       
           # Calculate the total variation distance
             totalvariation <- TotalVar(x,DiscretePrevision,n,s)


           # In order to decide if the parameter theta=j is minimizing:

           if (abs(totalvariation[[1]]-PreResult)<0.002) { # That is: If the total variation distance is equal to PreResult, ...
                theta <- cbind(theta,j)  # j is also minimizing. So, add j but do not change PreResult.
           }
           if (totalvariation[[1]]<=PreResult-0.002) { # That is: If the total variation distance is smaller than PreResult, ...
                theta <- j  # At the moment, j is minimizing and the previous parameters did not minimize. Forget them.
                PreResult <- totalvariation[[1]] # this is the smallest total variation distance at the moment
           }
       }       
       # As stated before: Do nothing ELSE in case of a FALSE if-expression.
  }

  list(theta,PreResult,runs)
}

