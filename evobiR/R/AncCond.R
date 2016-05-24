# to do:
#
#
################################################
#                                              #
#  Heath Blackmon & Richard Adams              #
#  Continuous value at nodes producing a       #
#  derived state: August 10 2015               #
#                                              #
################################################
##
## INPUT DATA
## trees: a phylo or multiPhylo object
## data: a dataframe with three collumns (tips, cont trait, disc trait)
## derived.state: a text string or numeric matching one of the 
##                entries in column 3 of data
## iterations the number of Monte Carlo simulations per tree 
##                used to calc p-value

AncCond <- function(trees, data, derived.state, iterations=1000){

  ## create named vector for disc trait for all taxa
  dt.vec <- data[, 3]
  names(dt.vec) <- data[, 1]
  
  ## create named vector for cont trait taxa not in derived state
  ct.data <- data[data[, 3] != derived.state, ] 
  ct.vec <- as.numeric(ct.data[, 2])
  names(ct.vec) <- ct.data[, 1]
  
  ## ASR for the continuous trait
  ## once for each tree since it is not stochastic
  anc.states.cont.trait <- list()
  if(class(trees) == "multiPhylo"){
    for(i in 1:length(trees)){
      cat(paste("Performing ASR of continuous character on tree:", i, "\n"))
      anc.states.cont.trait[[i]] <- anc.ML(trees[[i]], ct.vec, model="BM")
    }
  }
  if(class(trees) == "phylo"){
    anc.states.cont.trait[[1]] <- anc.ML(trees, ct.vec, model="BM")
  }
    
  ## ASR for discrete trait
  ## using stochastic mappings to nail down specific transition points
  anc.state.dt  <- list()
  if(class(trees) == "multiPhylo"){
    for(i in 1:length(trees)){
      #################################################################
      ##### pi argumentis is hard coded so derived = 0 and ancestral =1
      #################################################################
      temp.anc <- make.simmap(trees[[i]], 
                              dt.vec, 
                              model = matrix(c(0, 0, 1, 0), 2), 
                              nsim = 1, 
                              pi=c(1, 0), 
                              message=F)
      anc.state.dt[[i]] <- temp.anc
    }
  }
  if(class(trees) == "phylo"){
    temp.anc <- make.simmap(trees, 
                            dt.vec, 
                            model = matrix(c(0, 0, 1, 0), 2), 
                            nsim = 1, 
                            pi=c(1,0), 
                            message=F)
    anc.state.dt[[1]] <- temp.anc
  }
    
  ## Parse simmap to get producing nodes
  producing.nodes.list <- list()
  if(class(trees) == "multiPhylo"){
    for(i in 1:length(trees)){
      anc.state.bi <- anc.state.dt[[i]]
      ss_nodes <- anc.state.bi$mapped.edge[,1] > 0 & 
        anc.state.bi$mapped.edge[,2] > 0
      wanted_branches <- ss_nodes[ss_nodes==T]
      wanted_nodes <- names(wanted_branches)
      wanted_nodes <- gsub(",.*","",wanted_nodes)
      producing.nodes.list[[i]] <- unique(wanted_nodes)
    }
  }
  if(class(trees) == "phylo"){
    anc.state.bi <- anc.state.dt[[1]]
    ss_nodes <- anc.state.bi$mapped.edge[,1] > 0 & 
      anc.state.bi$mapped.edge[,2] > 0
    wanted_branches <- ss_nodes[ss_nodes==T]
    wanted_nodes <- names(wanted_branches)
    wanted_nodes <- gsub(",.*","",wanted_nodes)
    producing.nodes.list[[1]] <- unique(wanted_nodes)
  }
      
  ## get the mean ancestral value for the cont trait 
  ## at nodes producing the derived state marginalizing across trees
  hd.nodes <- vector()
  if(class(trees) == "multiPhylo"){
    for(i in 1:length(trees)){
      anc.states <- anc.states.cont.trait[[i]]
      producing.nodes <- producing.nodes.list[[i]]
      hd.nodes[i] <- mean(anc.states$ace[names(anc.states$ace) %in% 
                                           producing.nodes])
    }
  }
  if(class(trees) == "phylo"){
    anc.states <- anc.states.cont.trait[[1]]
    producing.nodes <- producing.nodes.list[[1]]
    hd.nodes[1] <- mean(anc.states$ace[names(anc.states$ace) %in% 
                                         producing.nodes])
  }
    hd.nodes <- mean(hd.nodes)

  ## Produce the null distribution of nodes in ancestral cond
  null.anc.nodes <- vector()
  counter <- 1
  number.of.trans <- vector()
  if(class(trees) == "multiPhylo"){
    for(i in 1:length(trees)){
      producing.nodes <- producing.nodes.list[[i]]
      number.of.trans[i] <- length(producing.nodes)
      anc.dt <- anc.state.dt[[i]]
      anc.ct <- anc.states.cont.trait[[i]]
      node.states <- describe.simmap(anc.dt)$states
      anc.cond.nodes <- anc.ct$ace[names(anc.ct$ace) %in% 
                                     names(node.states)[node.states != 
                                                          derived.state]]
      for(j in 1:iterations){
        null.anc.nodes[counter] <- mean(sample(anc.cond.nodes, 
                                               length(producing.nodes)))
        counter <- counter + 1
      }
    }
  }
  if(class(trees) == "phylo"){
      producing.nodes <- producing.nodes.list[[1]]
      number.of.trans[1] <- length(producing.nodes)
      anc.dt <- anc.state.dt[[1]]
      anc.ct <- anc.states.cont.trait[[1]]
      node.states <- describe.simmap(anc.dt)$states
      anc.cond.nodes <- anc.ct$ace[names(anc.ct$ace) %in% 
                                     names(node.states)[node.states != 
                                                          derived.state]]
      for(j in 1:iterations){
        null.anc.nodes[counter] <- mean(sample(anc.cond.nodes, 
                                               length(producing.nodes)))
        counter <- counter + 1
      }
    }
  
  ## plot results
  plot(density(null.anc.nodes), main="", xlab=colnames(data)[2], lwd=2)
  max.y <- range(density(null.anc.nodes)[2])[2]
  lines(x=c(hd.nodes,hd.nodes), y=c(0, max.y), col="red", lty=1, lwd=2)
  text(17,.4, label=paste("p-value", 
                          sum(null.anc.nodes<hd.nodes) / 
                            length(null.anc.nodes)))

  ## how many more extreme
  bigger <- (sum(null.anc.nodes >= hd.nodes)/
               length(null.anc.nodes)) * 100
  smaller <- (sum(null.anc.nodes <= hd.nodes)/
                length(null.anc.nodes)) * 100
  
  
  ## print results to terminal
  cat(paste("Derived State Mean Ancestral Cond:", hd.nodes, "\n"))
  cat(paste("Number of producing nodes:", mean(number.of.trans), "\n"))
  cat(paste("Mean of null dist:", mean(null.anc.nodes), "\n"))
  cat(paste("SD of null dist:", sd(null.anc.nodes), "\n"))
  cat(paste(bigger), "% of simulations had a derived state which arrose\n", 
      "with a mean continuous value larger than 
      inferred for the observed derived state\n")
  cat(paste(smaller), "% of simulations had a derived state which arrose\n", 
    "with a mean continuous value smaller than 
    inferred for the observed derived state\n")
  
  ## return results to user
  results <- list()
  results[[1]] <- hd.nodes
  results[[2]] <- number.of.trans
  results[[3]] <- null.anc.nodes
  results[[4]] <- bigger
  results[[5]] <- smaller
  names(results) <- c("OriginatingNodes", "NTrans", 
                      "NullDist", "bigger", "smaller")
  return(results)
}