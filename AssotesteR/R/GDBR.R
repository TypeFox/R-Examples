#'GDBR: Genomic Distance-Based Regression
#'
#'The Genomic Distance-Based Regression has been developed by Wessel et al
#'(2006). This approach captures genotype information across multiple loci
#'through a similarity measure between any two individuals. GDBR is unique in
#'its regression analysis relating variation in the measure of genomic
#'similarity to variation in their trait values. Note that this approach is
#'computationally expensive.
#'
#'The argument \code{distance} is used to specify the similarity distance.
#'\code{"IBS"} indicates Identity-By-Share, \code{"wIBS"} indicates weighted
#'IBS.
#'
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'@param distance character string indicating the type of distance to be used.
#'Possible options are "IBS" or "wIBS" (\code{distance="IBS"} by default)
#'@param weights optional numeric vector with weights for the genetic variants
#'(\code{NULL} by default)
#'@param perm positive integer indicating the number of permutations
#'(\code{NULL} by default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem gdbr.stat gdbr statistic
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, permutations, and selected distance
#'@returnItem name name of the statistic
#'@note This method is computationally expensive
#'@author Gaston Sanchez
#'@seealso \code{\link{SKAT}}
#'@references Wessel J, Schork NJ (2006) Generalized Genomic Distance-Based
#'Regression Methodology for Multilocus Association Analysis. \emph{The
#'American Journal of Human Genetics}, \bold{79}: 792-806 \cr
#'
#'Schaid DJ (2010) Genomic Similarity and Kernel Methods I: Advancements by
#'Building on Mathematical and Statistical Foundations. \emph{The American
#'Journal of Human Heredity}, \bold{70}: 109-131
#'@examples
#'
#'  \dontrun{
#'   
#'  # number of cases
#'  cases = 250
#'
#'  # number of controls
#'  controls = 250
#'
#'  # total (cases + controls)
#'  total = cases + controls
#'
#'  # phenotype vector
#'  phenotype = c(rep(1,cases), rep(0,controls))
#'
#'  # genotype matrix with 10 variants (random data)
#'  set.seed(123)
#'  genotype = matrix(rbinom(total*10, 2, 0.05), nrow=total, ncol=10)
#'
#'  # apply GDBR with 50 permutations
#'  # (it takes some time to run the permutations!)
#'  mygdbr = GDBR(phenotype, genotype, perm=50)
#'  mygdbr
#'  }
#'
GDBR <- 
function(y, X, distance="IBS", weights=NULL, perm=NULL)
{
  ## checking number of permutations
  if (!is.null(perm))
  {
    if (mode(perm) != "numeric" || length(perm) != 1
        || perm < 0 || (perm %% 1) !=0) 
    {
      warning("argument 'perm' incorrectly defined. Value perm=100 is used")
      perm = 100
    }
  } else perm = 0
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  # distance
  if (!(distance %in% c('IBS', 'wIBS')))
    stop(paste("\n", "Sorry =(   I don't recognize distance:", distance, "\n",
               "Choose one from: 'IBS' or 'wIBS'")) 	
  if (distance == "wIBS" && is.null(weights))
    stop("when using distance 'wIBS', you must provide a vector of 'weights'")
  if (!is.null(weights))
  {
    if(!is.vector(weights) || mode(weights) != "numeric" || any(weights < 0))
      stop("argument 'weights' must be a numeric vector with non-negative values")
    if(length(weights) != ncol(X)) 
      stop("length of 'weights' does not match number of columns in 'X'")
  }
  
  # how many individuals
  n = nrow(X)
  ## Genomic similarity matrix
  Sim = switch(distance, 
               "IBS" = gdbr_IBS(X),
               "wIBS" = gdbr_wIBS(X, weights)
               )
  ## distance matrix
  D = 1 - Sim
  A = (-1/2) * D^2
  I = diag(1, n)
  ## association matrix
  A.cen = scale(A, scale=FALSE)
  G = A.cen - rowMeans(A.cen)
  
  # get Fstat ((very computer-intensive!!!)
  gdbr.stat = my_gdbr_fstat(y, G)
  
  ## permutation procedure
  perm.pval = NA
  if (perm > 0)
  {
    I = diag(1, length(y))
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      y.perm = y[perm.sample] - mean(y)
      # get projection 'hat' matrix
      H = y.perm %*% solve((t(y.perm) %*% y.perm)) %*% t(y.perm)
      # calculate F statistic
      Fstat.num = sum(diag(H %*% G %*% H))
      Fstat.denom = sum(diag((I - H) %*% G %*% t(I - H)))
      x.perm[i] = Fstat.num / Fstat.denom
    }
    # p-value 
    perm.pval = sum(x.perm > gdbr.stat) / perm
  }
  
  ## results
  name = "GBDR: Genomic Distance-Based Regression"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm, distance)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms", "distance")	
  res = list(gdbr.stat = gdbr.stat, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}
