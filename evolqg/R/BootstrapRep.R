#' Bootstrap analysis via resampling
#'
#'
#'   Calculates the repeatability of the covariance matrix of the suplied data
#'   via bootstrap resampling
#'
#'   Samples with replacement are taken from the full population, a statistic calculated
#'   and compared to the full population statistic. 
#'
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param ComparisonFunc comparison function
#' @param iterations Number of resamples to take
#' @param sample.size Size of ressamples, default is the same size as ind.data 
#' @param correlation If TRUE, correlation matrix is used, else covariance matrix. 
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @return returns the mean repeatability, that is, the mean value of comparisons from samples to original statistic.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{MonteCarloStat}}, \code{\link{AlphaRep}}
#' @export
#' @examples
#' BootstrapRep(iris[,1:4], MantelCor, iterations = 5, correlation = TRUE)
#'              
#' BootstrapRep(iris[,1:4], RandomSkewers, iterations = 50)
#' 
#' BootstrapRep(iris[,1:4], KrzCor, iterations = 50, correlation = TRUE)
#' 
#' BootstrapRep(iris[,1:4], PCAsimilarity, iterations = 50)
#' 
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' #library(doParallel)
#' ##Windows:
#' #cl <- makeCluster(2)
#' #registerDoParallel(cl)
#' ##Mac and Linux:
#' #registerDoParallel(cores = 2)
#' #BootstrapRep(iris[,1:4], PCAsimilarity,
#' #             iterations = 5,
#' #             parallel = TRUE)
#' @keywords bootstrap
#' @keywords repetabilities


BootstrapRep <- function(ind.data,
                         ComparisonFunc,
                         iterations = 1000, 
                         sample.size = dim (ind.data)[1],
                         correlation = FALSE, 
                         parallel = FALSE){
  if(correlation)  {StatFunc <- cor; c2v <- cov2cor
  } else {StatFunc <- cov; c2v <- function(x) x}
  repeatability <- BootstrapStat(ind.data, iterations,
                                 ComparisonFunc = function(x, y) ComparisonFunc(c2v(x), 
                                                                                c2v(y)),
                                 StatFunc = StatFunc,
                                 sample.size = sample.size,
                                 parallel = parallel)
  return(mean(repeatability[,2]))
}


#' Non-Parametric population samples and statistic comparison
#'
#' Random populations are generated via ressampling 
#' using the suplied population. A statistic is calculated on the
#' random population and compared to the statistic calculated on the
#' original population. 
#'
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param iterations Number of resamples to take
#' @param ComparisonFunc comparison function
#' @param StatFunc Function for calculating the statistic
#' @param sample.size Size of ressamples, default is the same size as ind.data 
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @return returns the mean repeatability, that is, the mean value of comparisons from samples to original statistic.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}, \code{\link{AlphaRep}}
#' @rdname BootstrapStat
#' @export
#' @import plyr
#' @examples
#' cov.matrix <- RandomMatrix(5, 1, 1, 10)
#'
#' BootstrapStat(iris[,1:4], iterations = 50,
#'                ComparisonFunc = function(x, y) PCAsimilarity(x, y)[1],
#'                StatFunc = cov)
#'
#' #Calculating R2 confidence intervals
#' r2.dist <- BootstrapR2(iris[,1:4], 30)
#' quantile(r2.dist)
#' 
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' #library(doParallel)
#' ##Windows:
#' #cl <- makeCluster(2)
#' #registerDoParallel(cl)
#' ##Mac and Linux:
#' #registerDoParallel(cores = 2)
#' #BootstrapStat(iris[,1:4], iterations = 100,
#' #               ComparisonFunc = function(x, y) KrzCor(x, y)[1],
#' #               StatFunc = cov,
#' #               parallel = TRUE)
#' @keywords parametricsampling
#' @keywords montecarlo
#' @keywords repeatability
BootstrapStat <- function (ind.data, 
                           iterations,
                           ComparisonFunc, 
                           StatFunc, 
                           sample.size = dim (ind.data)[1],
                           parallel = FALSE){
  if(isSymmetric(as.matrix(ind.data))) stop("input appears to be a matrix, use residuals.")
  c.matrix <- StatFunc(ind.data)
  n_ind <- dim (ind.data)[1]
  populations  <- alply(1:iterations, 1,
                        function(x) ind.data[sample (1:n_ind, sample.size, TRUE),],
                        .parallel = parallel)
  comparisons <- ldply (populations, 
                        doComparisonBS, c.matrix, ComparisonFunc, StatFunc, 
                        ind.data, sample.size, n_ind,
                        .parallel = parallel)
  return (comparisons)
}

doComparisonBS <- function (x, c.matrix, ComparisonFunc, StatFunc, ind.data, sample.size, n_ind){
  while(TRUE){
    out = tryCatch(ComparisonFunc (c.matrix, StatFunc(x)), 
                   error = function(c) {
                     "Error in BootStrap sample, trying again"
                     return(NA)
                   })
    if(is.na(out))
      x = ind.data[sample (1:n_ind, sample.size, TRUE),]
    else
      break
  }
  return(out)
}

#' R2 confidence intervals by bootstrap resampling
#'
#' Random populations are generated by  resampling 
#' the suplied data or residuals. R2 is calculated on all the
#' random population's correlation matrices, provinding a distribution based on the original data.
#'
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param iterations Number of resamples to take
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @return returns a vector with the R2 for all populations
#' @author Diogo Melo Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}, \code{\link{AlphaRep}}
#' @export
#' @import plyr
#' @importFrom mvtnorm rmvnorm
#' @examples
#' r2.dist <- BootstrapR2(iris[,1:4], 30)
#' quantile(r2.dist)
#' @keywords bootstrap
#' @keywords integration
#' @keywords repeatability
BootstrapR2 <- function (ind.data, iterations = 1000, parallel = FALSE) {
  it.r2 <- BootstrapStat(ind.data, iterations,
                          ComparisonFunc = function(x, y) y,
                          StatFunc = function(x) CalcR2(cor(x)),
                          parallel = parallel)
  return (it.r2[,2])
}
