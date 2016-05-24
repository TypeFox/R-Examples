snphwe <- function(g) {
  gt <- table(g) ## assume sorting works, e.g. AA/AC/CC or 0/1/2 but not AA/AC/CA/CC
  if (length(gt) <= 2) return(p = 1.)
  if (length(gt) == 3) return(p = snphweCounts(gt[2], gt[1], gt[3]))
  ##warning("cannot calculate exact p-value for >3 genotypes")
  return(NA)
}

## This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
## Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
## Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76:(5) 887-893

## Original code of Wiggington, Cutler and Abecasis
## modified so that:
## 0. function name is snphweCounts instead of SNPHWE
## 1. error conditions cause stop with informative message
##    instead of returning -1
## 2. explicit return statement for computed p-value 
snphweCounts <- function(obs_hets, obs_hom1, obs_hom2) {
  if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
    stop("genotype counts must be non-negative")
  
  ## total number of genotypes
  N <- obs_hom1 + obs_hom2 + obs_hets
  
  ## rare homozygotes, common homozygotes
  obs_homr <- min(obs_hom1, obs_hom2)
  obs_homc <- max(obs_hom1, obs_hom2)
  
  ## number of rare allele copies
  rare  <- obs_homr * 2 + obs_hets
  
  ## Initialize probability array
  probs <- rep(0, 1 + rare)
  
  ## Find midpoint of the distribution
  mid <- floor(rare * ( 2 * N - rare) / (2 * N))
  if ( (mid %% 2) != (rare %% 2) ) mid <- mid + 1

  probs[mid + 1] <- 1.0
  mysum <- 1.0

  ## Calculate probablities from midpoint down 
  curr_hets <- mid
  curr_homr <- (rare - mid) / 2
  curr_homc <- N - curr_hets - curr_homr
  
  while ( curr_hets >=  2) {
    probs[curr_hets - 1]  <- probs[curr_hets + 1] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0)  * (curr_homc + 1.0))
    mysum <- mysum + probs[curr_hets - 1]
    
    ## 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
    curr_hets <- curr_hets - 2
    curr_homr <- curr_homr + 1
    curr_homc <- curr_homc + 1
  }    

  ## Calculate probabilities from midpoint up
  curr_hets <- mid
  curr_homr <- (rare - mid) / 2
  curr_homc <- N - curr_hets - curr_homr
   
  while ( curr_hets <= rare - 2) {
    probs[curr_hets + 3] <- probs[curr_hets + 1] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
    mysum <- mysum + probs[curr_hets + 3]
    
    ## add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
    curr_hets <- curr_hets + 2
    curr_homr <- curr_homr - 1
    curr_homc <- curr_homc - 1
  }    
 
  ## P-value calculation
  target <- probs[obs_hets + 1]

  ##plo <- min(1.0, sum(probs[1:obs_hets + 1]) / mysum)
  
  ##phi <- min(1.0, sum(probs[obs_hets + 1: rare + 1]) / mysum)
  
  ## This assignment is the last statement in the fuction to ensure 
  ## that it is used as the return value
  return(p = min(1.0, sum(probs[probs <= target])/ mysum))
}

