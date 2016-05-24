gen.Network <-
function(method="cohort", p, Nobs, Nsub, sparsity, REsize, REprob, REnoise){
  # wrapper function for generating random networks
  #
  # INPUT:
  #      - method: method to be used to simulate networks. Either "cohort" for networks as described in Monti et al 2015 or "Danaher" for networks as described in Danaher et al 2013
  #        The following parameters are only relevant to "cohort":
  #          - p: number of nodes
  #          - Nsub: number of subjects
  #          - sparsity: sparsity level in population network
  #          - REsize: number of edges that have additional random effects (must be smaller than choose(p,2)!)
  #          - REprob: probability of re-wiring 
  #          - REnoise: variance of random effects for edges
  #        The following parameters are only relevant to "Danaher":
  #          - compSize: size of each random component
  #          - compNum: number of indenpdent (ie unconnected) components)
  #          - popNum: number of populations (or subjects). By default 3 to recreate Danaher et al. (2014) sim
  #          - pow, m: power of preferential attachment and number of edges to add at each step (from barabasi.game function in igraph)
  #          - strUpper, strLower: define interval from which to sample edge strengths  
  
  if (method=="cohort"){
    Networks = genREnetworksEnforceTrans(p = p, Nsub = Nsub, sparsity = sparsity, REsize=REsize, REprob=REprob, REnoise = REnoise)
  } else if (method=="Danaher"){
    compNum = 10 # arbitrarily set number of components to three following Danaher example
    compSize = floor(p/compNum)
    if (p != compNum * compSize){
      warning(paste("number of nodes p must be divisible by 10 (number of components). Number of nodes has been changed from", p,"to", compSize*compNum))
    }
    if (!missing(Nsub)){
      if (Nsub != 3){
      warning("This method only generates networks for 3 subjects.")
      }
    }
    Nsub = 3 # abritrarily set to 3 subjects as well!
    p = compSize * compNum
    Networks = generateJGLnetworks(compSize = compSize, compNum = compNum, popNum = Nsub, pow=1, m=1, strUpper =.4, strLower = .1)
  } else {
    stop("Unrecognized network generation method. \n  Networks must be generated according to the cohort model or the Danaher model")
  }
  if (!missing(Nobs)){
    # also generate data:
    Dat = lapply(Networks$SubPres, FUN=function(x){ mvrnorm(n = Nobs, mu = rep(0,p), Sigma = solve(x))})
    
    if (method=="cohort"){
      SimDat = list(Networks=Networks$SubPres, Data=Dat, PopNet=Networks$K, RanNet=Networks$RE)
      class(SimDat) = "MNS"
      return(SimDat)  
    } else {
      SimDat = list(Networks=Networks$SubPres, Data=Dat)
      class(SimDat) = "MNS"
      return(SimDat)
    }
  } else {
    
    if (method=="cohort"){
      SimNet = list(Networks=Networks$SubPres, PopNet=Networks$K, RanNet=Networks$RE)
      class(SimNet) = "MNS"
      return(SimNet)  
    } else {
      SimNet = list(Networks=Networks$SubPres)
      class(SimNet) = "MNS"
      return(SimNet)
    }
  }
}
