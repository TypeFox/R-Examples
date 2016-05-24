#'WST: Weighted Score Test
#'
#'The WST method has been proposed by Wang and Elston (2007) and it can be seen
#'as a fixed effects method with transformed predictors based on Fourier
#'Transformations. WST is based on Fourier Transform (FT) coefficients to
#'globally test a set of correlated genetic variants (e.g. SNPs). The sequence
#'of genetic variants values is transformed into a sequence of numbers by
#'discrete FT, but only the real parts of the FT coefficients are taken into
#'account. A weighted score statistic of the FT components is calculated, which
#'follows a standard normal distribution under the null hypothesis
#'
#'This function does not allow missing genotypes
#'
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2. NO
#'missing data is allowed
#'@param perm positive integer indicating the number of permutations (100 by
#'default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem wst.stat wst statistic
#'@returnItem asym.pval asymptotic p-value
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{SCORE}}, \code{\link{SUM}}
#'@references Wang T, Elston C (2007) Improved Power by Use of a Weighted Score
#'Test for Linkage Disequilibrium Mapping. \emph{The American Journal of Human
#'Genetics}, \bold{80}: 353-360
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
#'  # apply WST with 500 permutations
#'  mywst = WST(phenotype, genotype, perm=500)
#'  mywst
#'  }
#'
WST <-
function(y, X, perm=100)
{
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  
  ## fast fourier transform
  X.new = t(mvfft(t(X)))
  ## take real part
  X.new = Re(X.new)
  ## getUV
  getuv = my_getUV(y, X.new)
  U = getuv$U
  V = getuv$V
  ## run score method
  wst.stat = my_wst_method(U, V)
  asym.pval = 1 - pnorm(wst.stat)
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    ymean = mean(y)
    for (i in 1:perm)
    {
      perm.sample <- sample(1:length(y))
      # center phenotype y
      y.perm = y[perm.sample] - ymean
      # get score vector
      U.perm = colSums(y.perm * X.new, na.rm=TRUE)
      perm.wst = my_wst_method(U.perm, V)
      x.perm[i] = perm.wst
    }
    # permuted p-value 
    perm.pval = sum(x.perm > wst.stat) / perm	
  }
  
  ## results
  name = "WST: Weighted Score Test"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(wst.stat = wst.stat, 
             asym.pval = asym.pval, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

