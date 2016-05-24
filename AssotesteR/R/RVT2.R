#'RVT2: Rare Variant Test 2 for dichotomous traits
#'
#'RVT2 is a collapsing method developed by Morris and Zeggini (2010) based on a
#'regression framework that models the phenotype as a function of a collapsed
#'summary of the variants. In the case of RVT@, the collapsed summary consists
#'of the presence or absence of at least one minor allele at any rare variant.
#'In this sense, RVT2 is an accumulation approach that regresses phenotype on a
#'genetic score, defined as the presence of at least one minor allele at any
#'rare variant
#'
#'If no variants are below the specified \code{maf} threshold, the function
#'will stop and return an error message \cr
#'
#'There is no imputation for the missing data. Missing values are simply
#'ignored in the computations.
#'
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param maf numeric value indicating minor allele frequency threshold for rare
#'variants (\code{maf=0.05} by default)
#'@param perm positive integer indicating the number of permutations (100 by
#'default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem rvt2.stat rvt2 statistic
#'@returnItem asym.pval asymptotic p-value
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, rare variants, maf threshold, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{RVT1}}
#'@references Morris AP, Zeggini E (2010) An Evaluation of Statistical
#'Approaches to Rare Variants Analysis in Genetic Association Studies.
#'\emph{Genetic Epidemiology}, \bold{34}: 188-193 \cr
#'
#'Asimit J, Zeggini E (2010) Rare Variant Association Analysis Methods for
#'Complex Traits. \emph{Annual Review of Genetics}, \bold{44}: 293-308
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
#'  # apply RVT2 with maf=0.05 and 500 permutations
#'  myrvt2 = RVT2(phenotype, genotype, maf=0.05)
#'  myrvt2 
#'  }
#'
RVT2 <-
function(y, X, maf=0.05, perm=100)
{
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  if (mode(maf) != "numeric" || length(maf) != 1
      || maf <= 0  || maf >= 1)
    stop("argument 'maf' must be a value between 0 and 1")
  
  ## get minor allele frequencies
  MAFs = colSums(X, na.rm=TRUE) / (2*nrow(X))
  ## are there any rare variants?
  rare = sum(MAFs < maf)
  if (rare == 0)
    stop(paste("\n", "Oops: No rare variants below maf=", 
               maf, " were detected. Try a larger maf", sep=""))
  ## get only rare variants
  X.new = X[ , MAFs < maf]
  ## convert to 0 and 1 (0: no rare copies;  1: rare copies)
  X.new[X.new == 2] = 1 
  ## create indicator dummy
  x.ind = rowSums(X.new, na.rm=TRUE)
  x.ind[x.ind != 0] = 1
  # center phenotype y
  y.new = y - mean(y)
  # get score vector U
  U = sum(y.new * x.ind)
  ## V
  xv = sum((x.ind - mean(x.ind))^2)
  V = mean(y) * (1 - mean(y)) * xv
  ## get score
  score = sum(U^2 / V)
  if (is.na(score) || is.infinite(score) || is.nan(score))
    score = 0
  asym.pval = 1 - pchisq(score, 1)
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    ymean = mean(y)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      # center phenotype y
      y.perm = y[perm.sample] - ymean
      # get score vector
      U.perm = sum(y.perm * x.ind)
      x.perm[i] = sum(U.perm^2 / V)
    }
    # permuted p-value 
    perm.pval = sum(x.perm > score) / perm	
  }
  
  ## results
  name = "RVT2: Rare Variant Test 2"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), rare, maf, perm)
  arg.spec = as.character(arg.spec)
  names(arg.spec) = c("cases", "controls", "variants", "rarevar", "maf", "n.perms")	
  res = list(rvt2.stat = score, 
             asym.pval = asym.pval, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

