#'SKAT: Sequence Kernel Association Test
#'
#'SKAT is a regression method to test for association between genetic variants
#'(common and rare) in a region. A score-based variance-component test.
#'
#'The argument \code{kernel} is used to specify the kernel function.
#'\code{"linear"} indicates the linear kernel, \code{"wlinear"} indicates a
#'weighted linear kernel, \code{"quadratic"} indicates the quadratic polynomial
#'kernel, \code{"IBS"} indicates Identity-By-Share, \code{"wIBS"} indicates
#'weighted IBS, and \code{"twowayx"} indicates a two-way interaction kernel.
#'\cr
#'
#'For the weighted kernels (\code{"wlinear"} and \code{"wIBS"}), there are two
#'options to get the weights. The default option (\code{weights=NULL}) involves
#'the calculation of the weights by taking into account the minor allele
#'frequency of the variants. In this case, the weights are calculated from a
#'Beta distribution with parameters \code{a} and \code{b}. The second option is
#'to specify the weights by providing a vector of weights for the variants; in
#'this case the length of the vector must equal the number of columns in
#'\code{X}. For more information see reference Wu et al (2011)
#'
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'@param kernel character string indicating the type of kernel to be used.
#'Possible options are "linear", "wlinear", "quadratic", "IBS", "wIBS",
#'"twowayx" (\code{kernel="linear"} by default)
#'@param weights optional numeric vector with weights for the genetic variants
#'(\code{NULL} by default)
#'@param a positive numeric value for the parameter \code{a} in the Beta
#'distribution (\code{a=1} by default)
#'@param b positive numeric vallue for the parameter \code{b} in the Beta
#'distribution (\code{b=25} by default)
#'@param perm positive integer indicating the number of permutations
#'(\code{NULL} by default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem skat.stat skat statistic
#'@returnItem asymp.pval asymptotic p-value of the applied statistic
#'(distributed as chi-square with df=1)
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, permutations, and selected kernel
#'@returnItem name name of the statistic
#'@note This method is computationally expensive
#'@author Gaston Sanchez
#'@seealso \code{\link{WSS}}
#'@references Wu MC, Kraft P, Epstein MP, Taylor DM, Chanock SJ, Hunter DJ, Lin
#'X (2010) Powerful SNP-Set Analysis for Case-Control Genome-wide Association
#'Studies. \emph{The American Journal of Human Genetics}, \bold{86}: 929-942
#'\cr
#'
#'Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X (2011) Rare-Variant Association
#'Testing for Sequencing Data with the Sequence Kernel Association Test.
#'\emph{The American Journal of Human Genetics}, \bold{89}: 82-93
#'@examples
#'
#'  \dontrun{
#'   
#'  # load data genodata
#'  data(genodata)
#'
#'  # phenotype (first column of genodata)
#'  pheno = genodata[,1]
#'
#'  # genotype (rest of columns of genodata)
#'  geno = genodata[,-1]
#'
#'  # apply SKAT with linear kernel 
#'  myskat.linear = SKAT(pheno, geno, kernel="linear")
#'  myskat.linear
#'
#'  # apply SKAT with weighted linear kernel
#'  # weights estimated from distribution beta(MAF, a=1, b=25)
#'  myskat.wlinear = SKAT(pheno, geno, kernel="wlinear", a=1, b=25)
#'  myskat.wlinear
#'
#'  # apply SKAT with quadratic kernel
#'  myskat.quad = SKAT(pheno, geno, kernel="quadratic")
#'  myskat.quad
#'
#'  # apply SKAT with IBS kernel
#'  myskat.ibs = SKAT(pheno, geno, kernel="IBS")
#'  myskat.ibs
#'  }
#'
SKAT <-
function(y, X, kernel="linear", weights=NULL, a=1, b=25, perm=NULL)
{
  ## checking argumetns
  if (!is.vector(y) || mode(y) != "numeric")
    stop("argument 'y' must be a numeric vector")
  if (any(is.na(y))) 
    stop("Sorry =(   No missing data allowed in argument 'y' ")	
  if (!all(y %in% c(0, 1)))
    stop("Sorry =(   argument 'y' must contain only 0 and 1")
  if(!is.matrix(X) & !is.data.frame(X))
    stop("argument 'X' must be a matrix or data.frame")
  if (nrow(X) != length(y)) 
    stop("'X' and 'y' have different lengths")
  if (!is.matrix(X)) X = as.matrix(X)
  if (!(kernel %in% c('linear', 'wlinear', 'quadratic', 'IBS', 'wIBS', 'twowayx')))
    stop(paste("\n", "Sorry =(   I don't recognize kernel:", kernel, "\n",
               "Choose one from: 'linear', 'wlinear', 'quadratic', 'IBS', 'wIBS', 'twowayx'")) 
  if (!is.null(weights))
  {
    if(length(weights) != ncol(X)) 
      stop("length of 'weights' does not match number of columns in 'X'")
    if(!is.vector(weights) || mode(weights) != "numeric" || any(weights) <= 0)
      stop("argument 'weights' must be a numeric vector with positive values")
  } 
  if (mode(a) != "numeric" || length(a) != 1 || a <= 0)
    stop("argument 'a' must be a positive number")
  if (mode(b) != "numeric" || length(b) != 1 || b <= 0)
    stop("argument 'b' must be a positive number")
  if (!is.null(perm))
  {
    if (mode(perm) != "numeric" || length(perm) != 1 || perm < 0 || (perm %% 1) !=0) 
    {
      warning("Argument 'perm' incorrectly defined. Value perm=NULL is used")
      perm = NULL
    }
  } else perm=0
  
  ## how many obs and variants
  n = nrow(X)
  p = ncol(X)
  
  ## get weights for linear-weighted or IBS-weighted
  if (kernel=="wlinear" || kernel=="wIBS" && is.null(weights))
  {   
    # MAF across cases and controls combined
    MAF <- colMeans(X, na.rm=TRUE) / 2
    # weights of the variants
    weights = dbeta(MAF, a, b)
  }	
  
  ## apply kernel
  K = switch(kernel,
             "linear" = X %*% t(X),
             "wlinear" = ((X %*% diag(weights^2)) %*% t(X) ), 
             "quadratic" = (X %*% t(X) + 1) ^ 2,
             "IBS" = kernel_IBS(X, n, p),
             "wIBS" = kernel_wIBS(X, n, p, weights),
             "twowayx" = kernel_twowayx(X, n, p)
             )
  
  ## score statistic Q (follows a Chi-square distr)
  y.new = y - mean(y)
  Q = t(y.new) %*% K %*% y.new / 2 
  
  ## parameters to estimate distr of Q    
  mu = rep(mean(y), n)
  V = diag(mu * (1-mu))
  ones = rep(1, n)
  P = V - V %*% ones %*% solve(t(ones) %*% V %*% ones) %*% t(ones) %*% V
  PK = P %*% K
  muQ = sum(diag(PK)) / 2
  itau = sum(PK * t(PK)) / 2
  itausig = sum(PK * t(P)) / 2
  isigma = sum(P^2) / 2
  iest = itau - ((itausig^2)/isigma)
  skale = iest / (2 * muQ)
  degfr = 2 * (muQ^2) / iest
  ## p-value
  skat.stat = as.numeric(Q)
  asym.pval = 1 - pchisq(skat.stat/skale, df=degfr)
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      y.perm = y[perm.sample]
      y.new = y.perm - mean(y.perm)
      Q.perm = t(y.new) %*% K %*% y.new / 2
      x.perm[i] = as.numeric(Q.perm)
    }
    # p-value 
    perm.pval = sum(x.perm > skat.stat) / perm
  } else perm = "NULL"   
  
  ## results
  name = "SKAT: Sequence Kernel Association Test"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm, kernel)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms", "kernel")
  res = list(skat.stat = skat.stat, 
             asym.pval = asym.pval, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

