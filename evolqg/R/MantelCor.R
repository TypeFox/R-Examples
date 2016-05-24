#' Compare matrices via Mantel Correlation
#'
#' Calculates correlation matrix correlation and significance via Mantel test.
#'
#' @param cor.x Single correlation matrix or list of correlation matrices.
#'
#' If single matrix is suplied, it is compared to cor.y.
#'
#' If list is suplied and no cor.y is suplied, all matrices
#' are compared.
#'
#' If cor.y is suplied, all matrices in list are compared to it.
#' @param cor.y First argument is compared to cor.y.
#' Optional if cor.x is a list.
#' @param permutations Number of permutations used in significance calculation.
#' @param mod Set TRUE to use mantel in testing modularity hypothesis. Will return
#' AVG+, AVG- and AVG Ratio based on binary hipotesis matrix.
#' @param MHI Indicates if modularity hypothesis test should use Modularity Hypothesis Index instead of AVG Ratio. Ignored if mod is FALSE.
#' @param repeat.vector Vector of repeatabilities for correlation correction.
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @param ... aditional arguments passed to other methods
#' @return If cor.x and cor.y are passed, returns matrix pearson
#' correlation and significance via Mantel permutations.
#'
#' If cor.x is a list of matrices and cor.y is passed, same as above, but for all matrices in cor.x.
#'
#' If only cor.x is passed, a matrix of MantelCor average
#' values and probabilities of all comparisons.
#' If repeat.vector is passed, comparison matrix is corrected above
#' diagonal and repeatabilities returned in diagonal.
#' @note If the significance is not needed, MatrixCor provides the 
#' correlation and skips the permutations, so it is much faster.
#' @export
#' @importFrom vegan mantel
#' @rdname MantelCor
#' @references http://en.wikipedia.org/wiki/Mantel_test
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{KrzCor}},\code{\link{RandomSkewers}},\code{\link{mantel}},\code{\link{RandomSkewers}},\code{\link{TestModularity}}
#' @examples
#' c1 <- RandomMatrix(10, 1, 1, 10)
#' c2 <- RandomMatrix(10, 1, 1, 10)
#' c3 <- RandomMatrix(10, 1, 1, 10)
#' MantelCor(cov2cor(c1), cov2cor(c2))
#' 
#' cov.list <- list(c1, c2, c3)
#' cor.list <- llply(list(c1, c2, c3), cov2cor)
#'
#' MantelCor(cor.list)
#'
#'# For repeatabilities we can use MatrixCor, which skips the significance calculation
#' reps <- unlist(lapply(cov.list, MonteCarloRep, 10, MatrixCor, correlation = TRUE))
#' MantelCor(llply(cor.list, repeat.vector = reps))
#'
#' c4 <- RandomMatrix(10)
#' MantelCor(cor.list, c4)
#' 
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' #library(doParallel)
#' ##Windows:
#' #cl <- makeCluster(2)
#' #registerDoParallel(cl)
#' ##Mac and Linux:
#' #registerDoParallel(cores = 2)
#' #MantelCor(cor.list, parallel = TRUE) 
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords randomskewers
MantelCor <- function (cor.x, cor.y, ...) UseMethod("MantelCor")

#' @rdname MantelCor
#' @method MantelCor default
#' @export
MantelCor.default <- function (cor.x, cor.y, permutations = 1000, mod = FALSE, MHI = FALSE, ...) {
  mantel.output <- mantel(cor.x, cor.y, permutations = permutations)
  correlation <- mantel.output$statistic
  prob <- mantel.output$signif
  if (mod){
    index <- cor.y[lower.tri(cor.y)]
    avg.plus <- mean (cor.x [lower.tri(cor.x)] [index != 0])
    avg.minus <- mean (cor.x [lower.tri(cor.x)] [index == 0])
    if(MHI){
      avg.index <- (avg.plus - avg.minus)/CalcICV(cor.x)
      output <- c(correlation, prob, avg.plus, avg.minus, avg.index)
      names(output) <- c("Rsquared", "Probability", "AVG+", "AVG-", "MHI")
    } else{
      avg.ratio <- avg.plus / avg.minus
      output <- c(correlation, prob, avg.plus, avg.minus, avg.ratio)
      names(output) <- c("Rsquared", "Probability", "AVG+", "AVG-", "AVG Ratio")
    }
  } else{
    if(sum(diag(cor.x)) != dim(cor.x)[1] | sum(diag(cor.y))!= dim(cor.y)[1])
      warning("Matrices do not appear to be correlation matrices. Use with caution.")
    output <- c(correlation, prob)
    names(output) <- c("Rsquared", "Probability")
  }
  return (output)
}

#' @rdname MantelCor
#' @method MantelCor list
#' @export
MantelCor.list <- function (cor.x, cor.y = NULL,
                            permutations = 1000, repeat.vector = NULL,
                            mod = FALSE, MHI = FALSE, parallel = FALSE, ...)
{
  if (is.null (cor.y)) {
    output <- ComparisonMap(cor.x,
                         function(x, cor.y) MantelCor(x, cor.y, permutations),
                         repeat.vector = repeat.vector,
                         parallel = parallel)
  } else{
    output <- SingleComparisonMap(cor.x, cor.y,
                               function(x, y) MantelCor(y,
                                                        x,
                                                        permutations, mod = mod, MHI = MHI),
                               parallel = parallel)
  }
  return(output)
}

#' @export
#' @rdname MantelCor
MatrixCor <- function (cor.x, cor.y, ...) UseMethod("MatrixCor")

#' @rdname MantelCor
#' @method MatrixCor default
#' @export
MatrixCor.default <- function (cor.x, cor.y, ...)                           
{
  if(sum(diag(cor.x)) != dim(cor.x)[1] | sum(diag(cor.y))!= dim(cor.y)[1])
    warning("Matrices do not appear to be correlation matrices. Use with caution.")
  cor(cor.x[lower.tri(cor.x)], cor.y[lower.tri(cor.y)], ...)
}

#' @rdname MantelCor
#' @method MatrixCor list
#' @export
MatrixCor.list <- function (cor.x, cor.y = NULL,
                            permutations = 1000, repeat.vector = NULL,
                            mod = FALSE, parallel = FALSE, ...)
{
  if (is.null (cor.y)) {
    output <- ComparisonMap(cor.x,
                            function(x, cor.y) c(MatrixCor(x, cor.y), NA),
                            repeat.vector = repeat.vector,
                            parallel = parallel)[[1]]
  } else{
    output <- SingleComparisonMap(cor.x, cor.y,
                                  function(x, y) MatrixCor(x, y),                                                    
                                  parallel = parallel)
  }
  return(output)
}
