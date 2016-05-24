#'CAST: Cohort Allelic Sums Test
#'
#'CAST is a pooled association test applied to discover if the difference in
#'the sums of allelic mutation frequencies in case and control cohorts is
#'greater than would be expected by chance. CAST works by first collapsing the
#'genotypes across rare variants to generate a super-variant. It then tests the
#'association between the trait and this new super-variant.
#'
#'If no variants are below the specified \code{maf} threshold, the function
#'will stop and return an error message \cr
#'
#'The argument \code{test="fisher"} involves a fisher exact test. Conversely,
#'\code{test="chisq"} indicates a chi-square test. \cr
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
#'@param test character string indicating the type of test to be applied.
#'Possible values are \code{"fisher"} and \code{"chisq"} (\code{test="fisher"}
#'by default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem cast.stat cast statistic
#'@returnItem asym.pval asymptotic p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, maf, and applied test
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{link{TTEST}}
#'@references Morgenthaler S, Thilly WG (2007) A strategy to discover genes
#'that carry multi-allelic or mono-allelic risk for common diseases: A cohort
#'allelic sums test (CAST). \emph{Mutation Research}, \bold{615}: 28-56
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
#'  phenotype = c(rep(1,cases), rep(0,controls))
#'
#'  # genotype matrix with 10 variants (random data)
#'  set.seed(1234)
#'  genotype = matrix(rbinom(total*10, 2, 0.051), nrow=total, ncol=10)
#' 
#'  # apply CAST with fisher exact test 
#'  mycast1 = CAST(phenotype, genotype, maf=0.05, test = "fisher")
#'  mycast1
#'
#'  # apply CAST with chi-square test 
#'  mycast2 = CAST(phenotype, genotype, maf=0.05, test = "chisq")
#'  mycast2
#'  }
#'
CAST <- 
function(y, X, maf=0.05, test="fisher")
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
  if (mode(maf) != "numeric" || length(maf) != 1
      || maf <= 0  || maf >= 1)
    stop("argument 'maf' must be a value between 0 and 1")
  if (!(test %in% c("fisher", "chisq")))
    stop("argument 'test' must be 'fisher' or 'chisq'")
  
  # number of individuals 'N' and variants 'M'
  N = nrow(X)
  ## get minor allele frequencies
  MAFs = colSums(X, na.rm=TRUE) / (2 * N)   
  ## are there any rare variants?
  rare = sum(MAFs < maf)
  if (rare == 0)
    stop(paste("\n", "Oops: No rare variants below maf=", 
               maf, " were detected. Try a larger maf", sep=""))
  # collapsing genotypes across variants
  X.collapsed = rowSums(X[ , MAFs < maf], na.rm=TRUE)	
  ## convert to 0 and 1 (0: no rare copies; 1: one or more rare copies)
  X.collapsed[X.collapsed != 0] = 1 
  
  # testing proportions between cases and controls
  props = table(y, X.collapsed)
  stat = switch(test, 
                "fisher" = unlist(fisher.test(props)[c(3,1)]),
                "chisq" = unlist(chisq.test(props)[c(1,3)])
                )
  names(stat) = NULL
  cast.stat = stat[1]
  pval = stat[2]
  
  ## results
  name = "CAST: Cohort Allelic Sums Test"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), rare, maf, test)
  names(arg.spec) = c("cases", "controls", "variants", "rarevar", "maf", "test")	
  res = list(cast.stat = cast.stat, 
             asym.pval = pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

