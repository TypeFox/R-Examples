#'RARECOVER Algorithm
#'
#'RARECOVER is an algorithm proposed by Bhatia et al (2010) that determines the
#'set of variants in a manner of forward variable selection: starting from a
#'null model without any genetic variants, genetic variants are selected one by
#'one based on their statistical significance and then added into the model
#'
#'The applied association test statistic (denoted as \code{XCORR} in Bhatia et
#'al, 2010) is based on the Pearsons chi-square statistic \cr
#'
#'The argument \code{maf} is used to specify the threshold of the minor allele
#'frequency for rare variants. By default, only variants below \code{maf=0.05}
#'are taken into account in the analysis. However, if all variants in \code{X}
#'are considered as rare variants, setting \code{maf=1} will consider them all
#'for the analysis \cr
#'
#'There is no imputation for the missing data. Missing values are simply
#'ignored in the computations.
#'
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param maf numeric value indicating the minor allele frequency threshold for
#'rare variants (\code{maf=0.05} by default)
#'@param dif numeric value between 0 and 1 as a threshold for the decision
#'criterion in the RARECOVER algorithm (default \code{dif=0.5})
#'@param perm positive integer indicating the number of permutations (100 by
#'default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem rc.stat rarecover statistic
#'@returnItem perm.pval permuted p-value
#'@returnItem set set of selected variants
#'@returnItem args descriptive information with number of controls, cases,
#'variants, rare variants, maf, number of selected variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{WSS}}
#'@references Bhatia G, Bansal V, Harismendy O, Schork NJ, Topol EJ, Frazer K,
#'Bafna V (2010) A Covering Method for Detecting Genetic Associations between
#'Rare Variants and Common Phenotypes. \emph{PLoS Computational Biology},
#'\bold{6(10)}: e1000954
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
#'  set.seed(1234)
#'  genotype = matrix(rbinom(total*10, 2, 0.051), nrow=total, ncol=10)
#'
#'  # apply RARECOVER with dif=0.05 and 500 permutations
#'  myrc = RARECOVER(phenotype, genotype, maf=0.05, perm=500)
#'  myrc
#'  }
#'
RARECOVER <- 
function(y, X, maf=0.05, dif=0.5, perm=100)
{
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  if (mode(maf)!= "numeric" || length(maf) != 1 || maf<=0 || maf>1)
    stop("argument 'maf' incorreclty defined; must be a value between 0 and 1")
  if (mode(dif) != "numeric" || length(dif) != 1
      || dif <= 0  || dif >= 1)
    stop("argument 'dif' must be a value between 0 and 1")
  
  ## get minor allele frequencies
  MAF = colMeans(X, na.rm=TRUE) / 2   
  ## how many variants < maf
  rare = sum(MAF < maf)
  if (rare == 0)
    stop(paste("Ooops!  No rare variants detected below maf=", maf, sep=""))
  ## collapsing
  if (rare == 1) 
  {   
    # if rare variants <= 1, then NO collapse is needed
    X.new = X
  } else {
    X.new = X[,MAF < maf]	   
  }
  ## change genotype 2 into 1
  X.new[X.new == 2] = 1
  rarecov = my_rarecov_method(y, X.new, dif)
  rc.stat = rarecov$stat
  chisq.pval = rarecov$pval
  rc.sel = rarecov$sel
  names(rc.sel) = NULL
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      rarecov.perm = my_rarecov_method(y[perm.sample], X.new, dif)
      x.perm[i] = rarecov.perm$stat
    }
    ## p-value
    perm.pval = sum(x.perm > rc.stat) / perm
  }
  
  ## results
  name = "RARECOVER Algorithm"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), rare, maf, length(rc.sel), perm)
  arg.spec = as.character(arg.spec)
  names(arg.spec) = c("cases", "controls", "variants", "rarevar", "maf", "varsel", "n.perms")
  sel.names = names(X)[rc.sel] 
  if (is.null(sel.names)) 
    sel.names = paste("var", rarecov$sel, sep="")
  names(rc.sel) = sel.names
  res = list(rc.stat = rc.stat, 
             perm.pval = perm.pval, 
             set = rc.sel, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}










