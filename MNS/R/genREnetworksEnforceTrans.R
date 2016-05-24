genREnetworksEnforceTrans <-
function(p, Nsub, sparsity, REsize, REprob, REnoise){
  # generate random effects networks.
  #
  # Here we build on the study on healthy subjects in ABIDE data which suggested that subject specific networks become more 
  # random ie transitivity (AKA clustering coefficient drops)
  #
  # Networks are generated as follows:
  #     1) Generate a single population network - preferably according to a "BA" structure (mimics structure found in EDA)
  #     2) define the variance network across subjects - this is variance network. Edges of this network may or may not be present. This is an ER network  by default
  #     3) for each subject, randomly chose which variance edges to be present (indep. according to prob REprob) and randomly remove an edge from their population network
  #     4) randomly sample edge strength at this edge for this subject (centered around 0, variance REnoise)
  #     5) ensure they're all pos definite
  #
  
  #
  # INPUT:
  #     - p: number of nodes
  #     - Nsub: number of subjects
  #     - sparsity: sparsity level in population network
  #     - REsize: number of edges that have additional random effects (must be smaller than choose(p,2)!)
  #     - REprob: probability of re-wiring 
  #     - REnoise: variance of random effects for edges
  #     - method: method to simulate population random network. Either Erdos Renyi (ER) or Barabasi-Albert (BA)
  #     - RElocation: one of "Random" and "Same". Defines location of random effects edges. If "Random", then random edges occur randomly. If "Same" then they occur at the same places as current edges!
  #
  # OUTPUT:
  #     -
  #     -
  #     -
  
  # generate population network:
  K = genLargerBAnetwork(p = p, sparsity = sparsity)
  TN = K$TN
  K = K$Pres
    
  # define locations of random variability across subjects!
  REii = diag(p)*0.5
  ii = sample( which(TN[upper.tri(TN)]==0)  , size=REsize, replace=FALSE) # make sure we avoid overlaps for now
  REii[upper.tri(REii)][ii] = 1
  REii = REii + t(REii)
  
  # now we're ready to simulate random networks for each subject!
  SubPres = vector("list", Nsub) # subject precision matrices
  SubPresTrue = vector("list", Nsub) # true network matrices (not guaranteed to be pos def!)
  
  for (i in 1:Nsub){
    posdef = FALSE
    while(!posdef){
      RandomLocs = diag(p)
      RandomLocs[upper.tri(RandomLocs)] = runif(n = choose(p,2), min = 0, max = 1) # generate more than we actually need but its cheap so who cares
      RandomLocs = RandomLocs + t(RandomLocs)
      RandomLocs = apply(RandomLocs, c(1,2), FUN=function(x){ifelse(x>(1-REprob),1,0)})
      
      RandomLocs = RandomLocs * REii
      varNum = sum(RandomLocs[upper.tri(RandomLocs)])
      
      subjectNoise = diag(p)*.5 #rnorm(p, 0, REnoise)
      subjectNoise[upper.tri(subjectNoise)] = runif(choose(p,2), min = REnoise/2, REnoise) #(subjectNoise %*% t(subjectNoise)) * RandomLocs
      subjectNoise = (subjectNoise+t(subjectNoise)) * RandomLocs
      
      subSign = matrix(0, ncol=p, nrow=p)
      subSign[upper.tri(subSign)] = sample(c(-1,1), size = choose(p,2), replace = TRUE)
      subSign = subSign + t(subSign) + diag(p)
      
      subjectNoise = subjectNoise * subSign
      
      subjectPop = diag(p)*.5
      subjectPop[upper.tri(subjectPop)] = TN[upper.tri(K)]
      subjectPop =  subjectPop + t(subjectPop)
      
      sudoMat = subjectPop+subjectNoise#ensure_pos_def(K+subjectNoise)
      diag(sudoMat)=1

      # normalize in same way as Danaher et al 2013:
      norm = apply(sudoMat,1,FUN=function(x){sum(abs(x))})
      sudoMat = apply(sudoMat,2, FUN=function(x){x/(1.5*norm)})
      
      sudoMat = (sudoMat + t(sudoMat))/2 # to make symmetric
      diag(sudoMat) = 1
      
      if (min(eigen(sudoMat)$values)>0){
        SubPres[[i]] = sudoMat
        posdef = TRUE
      }
    }
  }
  return(list(K=TN, Kgen=K, RE=REii, SubPres=SubPres))
}
