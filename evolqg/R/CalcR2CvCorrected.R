#' Corrected integration value
#'
#' Calculates the Young correction for integration, using bootstrap resampling
#'
#' @param ind.data Matrix of indiviual measurments, or ajusted linear model
#' @param cv.level Coeficient of variation level choosen for integration index ajustment in linear model. Defaults to 0.06.
#' @param iterations Number of resamples to take
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @param ... aditional arguments passed to other methods
#' @return List with adjusted integration indexes, fitted models and simulated distributions of integration indexes and mean coeficient of variation.
#' @references Young, N. M., Wagner, G. P., and Hallgrimsson, B. (2010).
#' Development and the evolvability of human limbs. Proceedings of the
#' National Academy of Sciences of the United States of America, 107(8),
#' 3400-5. doi:10.1073/pnas.0911856107
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{MeanMatrixStatistics}}, \code{\link{CalcR2}}
#' @rdname CalcR2CvCorrected
#' @export
#' @examples
#' integration.dist = CalcR2CvCorrected(iris[,1:4])
#'
#' #adjusted values
#' integration.dist[[1]]
#'
#' #ploting models
#' library(ggplot2)
#' ggplot(integration.dist$dist, aes(r2, mean_cv)) + geom_point() + 
#'        geom_smooth(method = 'lm', color= 'black') + theme_bw()
#'        
#' ggplot(integration.dist$dist, aes(eVals_cv, mean_cv)) + geom_point() + 
#'        geom_smooth(method = 'lm', color= 'black') + theme_bw()
#' @keywords correlation
#' @keywords integration

CalcR2CvCorrected  <- function(ind.data, ...) UseMethod("CalcR2CvCorrected")

#' @rdname CalcR2CvCorrected
#' @method CalcR2CvCorrected default
#' @export
CalcR2CvCorrected.default <- function (ind.data, cv.level = 0.06, iterations = 1000, parallel = FALSE, ...) {
  cv <- function (x) return (sd(x)/mean(x))
  Stats = function(x) {
    cov.matrix = var(x)
    cor.matrix = cov2cor(cov.matrix)
    return(c(CalcR2(cor.matrix), cv(eigen(cov.matrix)$values), mean (apply (x, 2, cv))))
  }
  it.stats <- BootstrapStat(ind.data, iterations,
                           ComparisonFunc = function(x, y) y,
                           StatFunc = Stats,
                           parallel = parallel)[,-1]
  colnames(it.stats) <- c("r2", "eVals_cv", "mean_cv")
  lm.r2 <- lm(it.stats[,1]~it.stats[,3])
  lm.eVals.cv <- lm(it.stats[,2]~it.stats[,3])
  adjusted.r2 <- lm.r2$coefficients %*% c(1, cv.level)
  adjusted.eVals.cv <- lm.eVals.cv$coefficients %*% c(1, cv.level)
  adjusted.integration  <-  c(adjusted.r2, adjusted.eVals.cv)
  names(adjusted.integration) <- c("r2", "eVals_cv")
  models <- list("r2" = lm.r2, "eVals_cv" = lm.eVals.cv)
  output <- list("adjusted.integration.index" = adjusted.integration, "models" = models, "dist" = it.stats)
  return (output)
}

#' @rdname CalcR2CvCorrected
#' @method CalcR2CvCorrected lm
#' @export
CalcR2CvCorrected.lm <- function (ind.data, cv.level = 0.06, iterations = 1000, ...) {
    cv <- function (x) return (sd(x)/mean(x))
    model <- ind.data
    ind.data <- model$model[[1]]
    fac <- model$model[-1]
    df <- model$df.residual
    res <- residuals (model)
    size <- dim (res)
    r2 <- mcv <- evar <- c()
    i <- 1
    while (i <= iterations) {
        current.sample <- sample ((size[1]), replace = TRUE)
        corr <- cov2cor (var (res[current.sample,]) * (size[1] - 1)/df)
        r2[i] <- mean (corr[lower.tri(corr)]^2)
        evar[i] <- cv (eigen(var (res[current.sample,] * (size[1] - 1)/df))$values)
        data <- data.frame(ind.data.s = ind.data[current.sample,], fac.s = fac[current.sample,])
        tmp1 <- table (data$fac.s)
        tmp2 <- ddply(data, ~ fac.s, function(x) apply(x[,1:size[2]], 2, cv))[-1]
        tmp3 <- rowMeans(tmp2)
        mcv[i] <- sum ((tmp3 * tmp1)/sum (tmp1))
        if (!is.na(mcv[i]))
            i <- i + 1
      }
    it.stats <- cbind(r2, evar, mcv)
    colnames(it.stats) <- c("r2", "eVals_cv", "mean_cv")
    lm.r2 <- lm(it.stats[,1]~it.stats[,3])
    lm.eVals.cv <- lm(it.stats[,2]~it.stats[,3])
    adjusted.r2 <- lm.r2$coefficients %*% c(1, cv.level)
    adjusted.eVals.cv <- lm.eVals.cv$coefficients %*% c(1, cv.level)
    adjusted.integration  <-  c(adjusted.r2, adjusted.eVals.cv)
    names(adjusted.integration) <- c("r2", "eVals_cv")
    models <- list("r2" = lm.r2, "eVals_cv" = lm.eVals.cv)
    output <- list("adjusted.integration.index" = adjusted.integration, "models" = models, "dist" = it.stats)
    return (output)
}
