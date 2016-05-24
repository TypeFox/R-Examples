generateJGLnetworks <-
function(compSize, compNum, popNum=3, pow=1, m=1, strUpper=.4, strLower=.1){
  # generate networks as described in Danaher et al. 2014
  #
  # the way this works is roughly as follows:
  # - we have a certain number of components (given by compNum), each of size compSize.
  # - nodes within each component follow a scale-free distribution (i.e. BA network)
  # - nodes across components are independent
  #
  # For example to recreate the simulations in Danaher et al (2014) we set compNum=10, compSize=50
  #
  # INPUT:
  #      - compSize: size of each random component
  #      - compNum: number of indenpdent (ie unconnected) components)
  #      - popNum: number of populations (or subjects). By default 3 to recreate Danaher et al. (2014) sim
  #      - pow, m: power of preferential attachment and number of edges to add at each step (from barabasi.game function in igraph)
  #      - strUpper, strLower: define interval from which to sample edge strengths
  #      -
  #
  # OUTPUT:
  #      - Adj: array with adjacency matrices
  #      - Sigma: array with covariance (sigma) matrices
  #      - 
  #
  #
  
  
  p = compSize * compNum # number of nodes (implicitly defined)
  ID = matrix(sample(1:p, size = p, replace = FALSE), ncol=compSize) # generate group IDs, each row a group ID
  wholeAdj = array(0, c(p,p,popNum))
  wholeSigma = array(0, c(p,p,popNum))
  
  compNet = vector("list", compNum) # each entry has a network and adj corresponding to connectivity for a component
  for (net in 1:compNum){
    compNet[[net]] = genSmallNet(p = compSize, pow = pow, m = m, strUpper = strUpper, strLower = strLower)      #as.matrix(get.adjacency(barabasi.game(n = p, power = pow, m = m, directed = FALSE)))
    # we give networks to all populations and then take them away!
    for (p in 1:popNum){
      wholeAdj[ID[net,], ID[net,],p] = compNet[[net]]$Adj
      wholeSigma[ID[net,], ID[net,],p] = compNet[[net]]$Sigma
    }
  }
  
  # now go through and delete components which should not be there!
  for (net in 1:(popNum-1)){
    wholeAdj[ID[net,], ID[net,], (seq(1,net))] = 0
    wholeSigma[ID[net,], ID[net,], (seq(1,net))] = 0
  }
  
  # reset diagonals to zero:
  for (i in 1:popNum) {diag(wholeSigma[,,i])=1}
  
  # convert to lists to put in same format as other method:
  wholeAdjList = lapply(vector("list", popNum), FUN=function(x){return(matrix(0, p,p))})
  wholeSigmaList = lapply(vector("list", popNum), FUN=function(x){return(matrix(0, p,p))})
  
  for (i in 1:popNum){ 
    wholeAdjList[[i]] = wholeAdj[,,i]
    wholeSigmaList[[i]] = wholeSigma[,,i]
  }
  
  return(list(Adj = wholeAdjList, SubPres=wholeSigmaList))
}
