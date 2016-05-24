
PPSDiscrete<-function(trees, MCMC, states, N=2) {
  qa <- 1                                                                          # just a counter
  qb <- 1                                                                          # second counter
  if(class(trees) == "multiPhylo"){
    ntaxa <- length(trees[[1]]$tip)                                                  #finds the number taxa present in your trees - assumes its the same across sample
  }
  if(class(trees) == "phylo"){
    ntaxa <- length(trees$tip)
  }
  if(class(trees) == "multiPhylo"){
    tree.vec<-MCMC[,1]
    MCMC <- MCMC[, 2:ncol(MCMC)]
  }
  simdata <- data.frame(matrix(ncol = N, nrow = ntaxa))                            #initializes the data.frame that will store all the simulated data
  if(class(trees) == "multiPhylo"){
    row.names(simdata) <- trees$tips[[1]]                                              #add taxa name for each row
  }
  if(class(trees) == "phylo"){
    row.names(simdata) <- trees$tips                                              #add taxa name for each row
  }
  for(rs in 1:(length(states))){                                                     #each loop here is for the different rootstates
    for(i in 1:(ceiling(states[qb]*N))){                                             #each loop here is for an individual dataset
      y <- sample(1:(nrow(MCMC)),1)                                                  #picks the index for the randomly chosen MCMC state
      if(class(trees) == "multiPhylo"){                                              #rselects the proprer tree from posterior sample to match the state of the MCMC draw
        tree <- trees[[tree.vec[y]]]
      }
      if(class(trees) == "phylo"){
        tree <- trees
      }
      matrix <- list(structure(c(as.numeric(MCMC[y, 1:ncol(MCMC)])),
                              .Dim = c((sqrt((ncol(MCMC)))), (sqrt((ncol(MCMC))))), 
                              .Names = NULL))                                        #This line builds the matrix describing the markov model governing the trait evolution should handle traits with any number of states
      tempdata <- sim.char(tree, matrix, nsim = 1, model = "discrete", root = rs)    #performs a single simulation with root state =1
      simdata[qa] <- tempdata[,,]                                                    #feeds the simulated data into the dataframe
      qa <- qa + 1                                                                   #advances the counter for the simdata array
    }
    qb <- qb + 1                                                                     #advances counter for root states
  }
  rownames(simdata) <- dimnames(tempdata)[[1]]
  for(i in 1:ncol(simdata)){
    colnames(simdata)[i] <- paste('dataset', i, sep='.')
  }
  PPS.data <- return(simdata)
}