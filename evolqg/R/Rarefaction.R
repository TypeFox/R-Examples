#' Rarefaction analysis via ressampling
#' 
#' Calculates the repeatability of a statistic of the data, such as
#' correlation or covariance matrix, via bootstrap resampling with
#' varying sample sizes, from 2 to the size of the original data.
#'
#' Samples of various sizes, with replacement, are taken from the full population, a statistic calculated
#' and compared to the full population statistic. 
#'
#' A specialized ploting function displays the results in publication quality.
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param ComparisonFunc comparison function
#' @param ... Aditional arguments passed to ComparisonFunc
#' @param num.reps number of populations sampled per sample size
#' @param correlation If TRUE, correlation matrix is used, else covariance matrix. MantelCor always uses correlation matrix.
#' @param replace If true, samples are taken with replacement
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @details Bootstraping may be misleading with very small sample sizes. Use with caution if original sample sizes are small.
#' @return returns the mean value of comparisons from samples to original statistic, for all sample sizes.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}
#' @rdname Rarefaction
#' @export
#' @import plyr
#' @importFrom ggplot2 ggplot aes_string geom_boxplot scale_x_continuous scale_y_continuous theme_bw
#' @importFrom reshape2 melt
#' @examples
#' ind.data <- iris[1:50,1:4]
#' 
#' results.RS <- Rarefaction(ind.data, PCAsimilarity, num.reps = 5)
#' results.Mantel <- Rarefaction(ind.data, MatrixCor, correlation = TRUE, num.reps = 5)
#' results.KrzCov <- Rarefaction(ind.data, KrzCor, num.reps = 5)
#' results.PCA <- Rarefaction(ind.data, PCAsimilarity, num.reps = 5)
#' 
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' #library(doParallel)
#' ##Windows:
#' #cl <- makeCluster(2)
#' #registerDoParallel(cl)
#' ##Mac and Linux:
#' #registerDoParallel(cores = 2)
#' #results.KrzCov <- Rarefaction(ind.data, KrzCor, num.reps = 5, parallel = TRUE)
#' 
#' #Easy access
#' library(reshape2)
#' melt(results.RS)
#'
#' @keywords rarefaction
#' @keywords bootstap
#' @keywords repeatability
Rarefaction <- function(ind.data,
                        ComparisonFunc,
                        ...,
                        num.reps = 10,
                        correlation = FALSE, 
                        replace = FALSE,
                        parallel = FALSE){
  if(correlation)  {StatFunc <- cor; c2v <- cov2cor
  } else {StatFunc <- cov; c2v <- function(x) x}
  rarefaction.list <- RarefactionStat(ind.data,
                                      StatFunc = StatFunc,
                                      ComparisonFunc = function(x, y) ComparisonFunc(c2v(x), 
                                                                                     c2v(y), ...),
                                      num.reps = num.reps,
                                      parallel = parallel)
  return(rarefaction.list)
}

#' Non-Parametric rarefacted population samples and statistic comparison
#' 
#' Calculates the repeatability of a statistic of the data, such as
#' correlation or covariance matrix, via resampling with
#' varying sample sizes, from 2 to the size of the original data.
#'
#' Samples of various sizes, without replacement, are taken from the full population, a statistic calculated and compared to the full population statistic. 
#'
#' A specialized ploting function displays the results in publication quality.
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param StatFunc Function for calculating the statistic
#' @param ComparisonFunc comparison function
#' @param ... Aditional arguments passed to ComparisonFunc
#' @param num.reps number of populations sampled per sample size
#' @param replace If true, samples are taken with replacement
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @details Bootstraping may be misleading with very small sample sizes. Use with caution.
#' @return returns the mean value of comparisons from samples to original statistic, for all sample sizes.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}
#' @export
#' @import plyr
#' @importFrom ggplot2 ggplot aes_string layer scale_x_continuous scale_y_continuous theme_bw
#' @importFrom reshape2 melt
#' @examples
#' ind.data <- iris[1:50,1:4]
#' 
#' #Can be used to calculate any statistic via Rarefaction, not just comparisons
#' #Integration, for instanse:
#' results.R2 <- RarefactionStat(ind.data, cor, function(x, y) CalcR2(y), num.reps = 5)
#' 
#' #Easy access
#' library(reshape2)
#' melt(results.R2)
#' 
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' #library(doParallel)
#' ##Windows:
#' #cl <- makeCluster(2)
#' #registerDoParallel(cl)
#' ##Mac and Linux:
#' #registerDoParallel(cores = 2)
#' #results.R2 <- RarefactionStat(ind.data, cor, function(x, y) CalcR2(y), parallel = TRUE)
#'
#' @keywords rarefaction
#' @keywords bootstap
#' @keywords repeatability
RarefactionStat <- function(ind.data,
                            StatFunc,
                            ComparisonFunc,
                            ...,
                            num.reps = 10,
                            replace = FALSE,
                            parallel = FALSE)
{
  if(isSymmetric(as.matrix(ind.data))) stop("input appears to be a matrix, use residuals.")
  observed.stat = StatFunc(ind.data)
  num.ind = dim(ind.data)[1]
  MapStatFunc <- function(n){
    SampleFunction <- function(x){
      while(TRUE){
        local.sample = sample(1:num.ind, n, replace=replace)
        output <- tryCatch(StatFunc(ind.data[local.sample,]), warning=function(w) w)
        if(!is(output, "warning"))
          return(output)
      }
    }
    return(alply(1:num.reps, 1, SampleFunction))
  }
  sample.stats = alply(3:num.ind, 1, MapStatFunc, .parallel = parallel)
  names(sample.stats) <- 3:num.ind
  MapComparisonFunc <- function(stat.list, ...){
    ldply(stat.list, function(x, ...) ComparisonFunc(observed.stat, x, ...))[,2]
  }
  comparison.list = llply(sample.stats, MapComparisonFunc, .parallel = parallel)
  return(comparison.list)
}

#' Plot Rarefaction analysis
#' 
#' A specialized ploting function displays the results from Rarefaction functions in publication quality.
#' 
#' @param comparison.list output from rarefaction functions can be used in ploting
#' @param y.axis Y axis lable in plot
#' @param x.axis Y axis lable in plot
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}
#' @export
#' @import plyr
#' @importFrom ggplot2 ggplot aes_string layer scale_x_continuous scale_y_continuous theme_bw
#' @importFrom reshape2 melt
#' @examples
#' ind.data <- iris[1:50,1:4]
#' 
#' results.RS <- Rarefaction(ind.data, PCAsimilarity, num.reps = 5)
#' results.Mantel <- Rarefaction(ind.data, MatrixCor, correlation = TRUE, num.reps = 5)
#' results.KrzCov <- Rarefaction(ind.data, KrzCor, num.reps = 5)
#' results.PCA <- Rarefaction(ind.data, PCAsimilarity, num.reps = 5)
#'
#' #Plotting using ggplot2
#' a <- PlotRarefaction(results.RS, "Random Skewers")
#' b <- PlotRarefaction(results.Mantel, "Mantel")
#' c <- PlotRarefaction(results.KrzCov, "KrzCor")
#' d <- PlotRarefaction(results.PCA, "PCAsimilarity")
#'
#' library(grid)
#' grid.newpage()
#' pushViewport(viewport(layout = grid.layout(2, 2)))
#' vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#' print(a, vp = vplayout(1, 1))
#' print(b, vp = vplayout(1, 2))
#' print(c, vp = vplayout(2, 1))
#' print(d, vp = vplayout(2, 2))
#'
#' @keywords rarefaction
#' @keywords bootstap
#' @keywords repeatability
#' @export
PlotRarefaction <- function(comparison.list, y.axis = "Statistic", 
                            x.axis = "Number of sampled specimens"){
  plot.df = melt(comparison.list, value.name = 'avg.corr')
  plot.df = as.data.frame(lapply(plot.df, as.numeric))
  rarefaction.plot = ggplot(plot.df, aes_string('L1', 'avg.corr', group = 'L1')) +
  geom_boxplot() +
  scale_x_continuous(x.axis) +
  scale_y_continuous(y.axis) +
  theme_bw()
  return(rarefaction.plot)
}
