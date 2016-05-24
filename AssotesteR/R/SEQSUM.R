#'SEQSUM: Sequential Sum Score Test
#'
#'SEQSUM has been proposed by Basu and Pan (2011) as a modification of the Sum
#'test based on a model selection approach, following a similar philosophy as
#'the CARV and RARECOVER methods. Assuming that there are M variants, the main
#'idea behind the Sequential Sum test is to associate a sign to each variant
#'indicating whether it has a positive effect or a negative effect. In other
#'words, the purpose is to give signs to the variants so they reflect their
#'effect (positive or negative).
#'
#'There is no imputation for the missing data. Missing values are simply
#'ignored in the computations.
#'
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param perm positive integer indicating the number of permutations (100 by
#'default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem seqsum.stat seqsum statistic
#'@returnItem perm.pval permuted p-value
#'@returnItem signs a numeric vector with signs for the variants (1=positive,
#'-1=negative)
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{SCORE}}, \code{\link{SUM}}
#'@references Basu S, Pan W (2011) Comparison of Statistical Tests for Disease
#'Association with Rare Variants. \emph{Genetic Epidemiology}, \bold{35}:
#'606-619 \cr
#'@examples
#'
#'  \dontrun{
#'
#'  # number of cases
#'  cases = 500
#'
#'  # number of controls 
#'  controls = 500
#'
#'  # total (cases + controls)
#'  total = cases + controls
#'
#'  # phenotype vector
#'  phenotype = c(rep(1, cases), rep(0, controls))
#'
#'  # genotype matrix with 10 variants (random data)
#'  set.seed(123)
#'  genotype = matrix(rbinom(total*10, 2, 0.05), nrow=total, ncol=10)
#'
#'  # apply SEQSUM with 500 permutations
#'  myseq = SEQSUM(phenotype, genotype, perm=500)
#'  myseq
#'  }
#'
SEQSUM <- 
function(y, X, perm=100)
{
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  
  ## vector for storing stastitic and signs
  seqsum.stat = rep(0, ncol(X))
  signs = rep(-1, ncol(X))
  
  ## start with first variant
  X.new = X
  ## positive effect
  seqsum.stat = my_uni_score(y, rowSums(X.new, na.rm=TRUE))	
  ## negative effect
  X.new[,1] = (-1) * X[,1]
  aux.neg = my_uni_score(y, rowSums(X.new, na.rm=TRUE))
  ## who is the best
  if (seqsum.stat > aux.neg)
  {
    X.new[,1] = X[,1]
    signs[1] = 1	
  } else seqsum.stat = aux.neg
  
  ## continue with the other variants
  for (j in 2:ncol(X))
  {
    ## negative effect model
    X.new[,j] = (-1) * X[,j]
    aux.neg = my_uni_score(y, rowSums(X.new, na.rm=TRUE))
    # who is the best?
    if (seqsum.stat > aux.neg)
    {
      X.new[,j] = X[,j]
      signs[j] = 1
    } else seqsum.stat = aux.neg
  }
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      y.perm = y[perm.sample]
      ## start with first variant
      X.new = X
      perm.seqsum = my_uni_score(y.perm, rowSums(X.new, na.rm=TRUE))
      ## negative effect
      X.new[,1] = (-1) * X[,1]
      aux.neg = my_uni_score(y.perm, rowSums(X.new, na.rm=TRUE))
      if (perm.seqsum > aux.neg) 
      {
        X.new[,1] = X[,1]
      } else perm.seqsum = aux.neg
      ## continue with the other variants
      for (j in 2:ncol(X))
      {
        X.new[,j] = (-1) * X[,j]
        aux.neg = my_uni_score(y.perm, rowSums(X.new, na.rm=TRUE))
        if (perm.seqsum > aux.neg) 
        {
          X.new[,j] = X[,j]
        } else perm.seqsum = aux.neg
      }
      x.perm[i] = perm.seqsum
    }
    # p-value 
    perm.pval = sum(x.perm > seqsum.stat) / perm	
  }
  
  ## results
  name = "SEQSUM: Sequential Sum Test"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(seqsum.stat = seqsum.stat, 
             perm.pval = perm.pval, 
             signs = signs, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}
