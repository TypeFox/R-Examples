# Randomly perturb a supplied network.
# <perturb="functions"> perturbs the functions associated with the genes directly.
# <perturb="states"> perturbs a maximum of <numStates> states 
# in the transition table resulting from the functions.
# <method="bitflip"> randomly flips up to <maxNumBits> in the functions or states.
# <method="shuffle"> randomly permutes the bits in the functions or states.
# If <simplify> is set, the perturbed network is simplified to remove irrelevant input functions.
# If <excludeFixed> is set, fixed genes are excluded from the perturbations and stay as they are.
perturbNetwork <- function(network,perturb=c("functions","transitions"),method=c("bitflip","shuffle"),
      simplify=(perturb[1]!="functions"),readableFunctions=FALSE,excludeFixed=TRUE,
      maxNumBits=1,numStates=max(1,2^length(network$genes)/100))
{
  stopifnot(inherits(network,"BooleanNetwork") | inherits(network,"ProbabilisticBooleanNetwork"))

  fixedGenes <- which(network$fixed != -1)

  if (length(perturb) == 1 && perturb == "states")
  {
    warning("perturb=\"states\" is deprecated. Use perturb=\"transitions\" instead!")
    perturb <- "transitions"
  }
  
  if (inherits(network,"BooleanNetwork"))
  # deterministic network
  {
    switch(match.arg(perturb,c("functions","transitions")),
      functions=
        switch(match.arg(method,c("bitflip","shuffle")),
        bitflip =
          {
            # choose the function to be perturbed
            if (length(fixedGenes) > 0 & excludeFixed)
              functionIdx <- sample((seq_along(network$interactions))[-fixedGenes],size=1)
            else
              functionIdx <- sample(seq_along(network$interactions),size=1)
          
            # choose the indices of the truth table to be flipped  
            flipIndices <- sample(seq_along(network$interactions[[functionIdx]]$func),
                                  size=runif(n=1,min=1,
                                  max=min(maxNumBits,
                                  length(network$interactions[[functionIdx]]$func))),
                                  replace=FALSE)
            # flip the bits
            network$interactions[[functionIdx]]$func[flipIndices] <- 
              as.integer(!network$interactions[[functionIdx]]$func[flipIndices])
            network$interactions[[functionIdx]]$expression <-
              getInteractionString(readableFunctions,
                     network$interactions[[functionIdx]]$func,
                     network$genes[network$interactions[[functionIdx]]$input])
          },    
        shuffle=
          {
            # choose the function to be perturbed
            if (length(fixedGenes) > 0 & excludeFixed)
              functionIdx <- sample((seq_along(network$interactions))[-fixedGenes],size=1)
            else
              functionIdx <- sample(seq_along(network$interactions),size=1)

            # draw a random permutation of bit positions
            flipIndices <- sample(seq_along(network$interactions[[functionIdx]]$func),
                  size=length(network$interactions[[functionIdx]]$func),
                  replace=FALSE)
                
            # permute the bits
            network$interactions[[functionIdx]]$func <- 
              network$interactions[[functionIdx]]$func[flipIndices]
            network$interactions[[functionIdx]]$expression <-
              getInteractionString(readableFunctions,
                     network$interactions[[functionIdx]]$func,
                     network$genes[network$interactions[[functionIdx]]$input])
                    
          },
        stop("'method' must be one of \"bitflip\",\"shuffle\"")),
    
    transitions =
      {  
        # turn off gene fixing - otherwise reverse-engineering of the transition table is not possible
      
        oldFixed <- network$fixed
        network$fixed <- rep(-1,length(network$genes))
        names(network$fixed) <- network$genes
      
        # calculate transition table
        table <- t(sapply(getAttractors(network)$stateInfo$table,dec2bin,length(network$genes)))
      
        # determine the states to be perturbed
        statesToChange <- sample(seq_len(nrow(table)),min(numStates,nrow(table)),replace=FALSE)
      
        lapply(statesToChange,function(state)
          {
            # choose the indices of the states that are allowed to be changed
            flipIndices <- seq_along(network$genes)
            if (length(fixedGenes) > 0 & excludeFixed)
              flipIndices <- flipIndices[-fixedGenes]
            
            switch(match.arg(method,c("bitflip","shuffle")),
            bitflip =
              {
                # choose the actual indices to be changed
                flipIndex <- sample(flipIndices,
                                    size=runif(n=1,min=1,
                                    max=min(maxNumBits,
                                    length(flipIndices))),
                                    replace=FALSE)
                    
                # flip the bits at these positions
                table[state,flipIndex] <<- 
                  as.integer(!table[state,flipIndex])
        
              },
            shuffle =
              {
                # determine a permutation of the bit indices
                flipIndex <- sample(flipIndices,
                                    size=length(flipIndices),
                                    replace=FALSE)
                    
                # permute the state
                table[state,] <<- 
                  as.integer(table[state,flipIndex])
              },
              stop("'method' must be one of \"bitflip\",\"shuffle\"")        
            )
          NULL})
        
        # restore network by assigning the columns of the state table to the corresponding genes
        network$interactions <- apply(table,2,function(gene)
                {
                  input = seq_along(network$genes)
                
                  # encoding is reversed in the transition table
                  input <- input[length(input):1]
                  list(input=input,
                       func=gene,
                       expression= getInteractionString(readableFunctions,
                       gene,
                       network$genes[input]))
                })
        # reactivate fixed genes
        network$fixed <- oldFixed
      
      
      },  
      stop("'perturb' must be one of \"functions\",\"transitions\""))
  }
  else
  # probabilistic network
  {
    if (match.arg(perturb) != "functions")
      stop("In probabilistic Boolean networks, only perturb=functions is allowed!")

      switch(match.arg(method,c("bitflip","shuffle")),
      bitflip =
        {
          # choose the gene and the function to be perturbed
      
          if (length(fixedGenes) > 0 & excludeFixed)
            geneIdx <- sample((seq_along(network$interactions))[-fixedGenes],size=1)
          else
            geneIdx <- sample(seq_along(network$interactions),size=1)
                
          functionIdx <- sample(seq_along(network$interactions[[geneIdx]]),size=1)

          # choose the indices of the truth table to be flipped  
          flipIndices <- sample(seq_along(network$interactions[[geneIdx]][[functionIdx]]$func),
                              size=runif(n=1,min=1,
                              max=min(maxNumBits,
                              length(network$interactions[[geneIdx]][[functionIdx]]$func))),
                              replace=FALSE)
          # flip the bits
          network$interactions[[geneIdx]][[functionIdx]]$func[flipIndices] <- 
            as.integer(!network$interactions[[geneIdx]][[functionIdx]]$func[flipIndices])
          network$interactions[[geneIdx]][[functionIdx]]$expression <-
            getInteractionString(readableFunctions,
                                 network$interactions[[geneIdx]][[functionIdx]]$func,
                                 network$genes[network$interactions[[geneIdx]][[functionIdx]]$input])
        },    
      shuffle=
        {
          # choose the function to be perturbed
          if (length(fixedGenes) > 0 & excludeFixed)
            geneIdx <- sample((seq_along(network$interactions))[-fixedGenes],size=1)
          else
            geneIdx <- sample(seq_along(network$interactions),size=1)
                
          functionIdx <- sample(seq_along(network$interactions[[geneIdx]]),size=1)


          # draw a random permutation of bit positions
          flipIndices <- sample(seq_along(network$interactions[[geneIdx]][[functionIdx]]$func),
                                size=length(network$interactions[[geneIdx]][[functionIdx]]$func),
                                replace=FALSE)
                
          # permute the bits
          network$interactions[[geneIdx]][[functionIdx]]$func <- 
            network$interactions[[geneIdx]][[functionIdx]]$func[flipIndices]
          network$interactions[[geneIdx]][[functionIdx]]$expression <-
            getInteractionString(readableFunctions,
                                 network$interactions[[geneIdx]][[functionIdx]]$func,
                                 network$genes[network$interactions[[geneIdx]][[functionIdx]]$input])
                  
        },
      stop("'method' must be one of \"bitflip\",\"shuffle\""))
    
  }
  # simplify the network if necessary
  if (simplify)
    network <- simplifyNetwork(network,readableFunctions)
  return(network)
}
