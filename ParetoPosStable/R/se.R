#' @name se
#' @aliases se se.default print.se
#' @title Approximated standard errors of Pareto Positive Stable (PPS) parameter estimates
#' @description It approximates the stantard errors of PPS parameter estimates by bootstrapping.
#' @param PPSfit a \code{PPSfit} Object, typically from \code{PPS.fit()}.
#' @param k the number of steps in the bootstrapping procedure.
#' @param parallel A logical argument specifying if parallelization is allowed in the bootstrapping procedure.
#' @param ncores is the number of cores that we use if parallel is TRUE
#' @param \dots other arguments.
#' 
#' @details
#' The function simulates \code{k} samples from the model given in the \code{PPSfit} argument, fits them with the same method of estimation and uses the parameter estimates to approximate the standard errors.
#'
#' @return A list with the standard errors.
#' 
#' @references Sarabia, J.M and Prieto, F. (2009). The Pareto-positive stable distribution: A new descriptive model for city size data, \emph{Physica A: Statistical Mechanics and its Applications}, \bold{388}(19), 4179-4191.
#' 
#' @seealso \code{\link{PPS.fit}}
#'
#' @examples
#' x <- rPPS(50, 1.2, 100, 2.3)
#' fit <- PPS.fit(x)
#' coef(fit)
#' se(fit, k = 50)

#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster clusterExport detectCores
#' @importFrom foreach foreach %dopar% %do%

#' @export
se <- function(PPSfit, k = 2000, parallel = TRUE, ncores = 2, ...)
  UseMethod("se")

#' @export
se.default <- 
  function(PPSfit, k = 2000, parallel = TRUE, ncores = 2, ...){
    if (is.null(PPSfit$Pareto)) PPSfit$Pareto <- FALSE
    if (PPSfit$Pareto == FALSE){
      if (is.null(PPSfit$sigma)){
        if (parallel == TRUE){
          cl <- makeCluster(ncores)
          registerDoParallel(cl)
          clusterExport(cl, list("PPSfit"), envir = environment())
          muestra <- foreach (i = 1 : k, .combine = cbind, .multicombine = TRUE, .inorder = FALSE, .packages = c("ParetoPosStable")) %dopar% {
            converged <- FALSE
            datos.sim <- sample(x = PPSfit$obs, size = PPSfit$n, replace = TRUE)
            ajuste <- PPS.fit(datos.sim, estim.method = PPSfit$estim.method, start = PPSfit$estimate)
            c(ajuste$estimate[[1]], ajuste$estimate[[2]], ajuste$estimate[[3]])
          }
          stopCluster(cl)
        }
        else{#if no parallel
          muestra <- foreach (i = 1 : k, .combine = cbind, .multicombine = TRUE, .inorder = FALSE, .packages = c("ParetoPosStable")) %do% {
            converged <- FALSE
            datos.sim <- sample(x = PPSfit$obs, size = PPSfit$n, replace = TRUE)
            ajuste <- PPS.fit(datos.sim, estim.method = PPSfit$estim.method, start = PPSfit$estimate)
            c(ajuste$estimate[[1]], ajuste$estimate[[2]], ajuste$estimate[[3]])
          }
        }
        se.lambda <- sqrt(sum((muestra[1,] - PPSfit$estimate[[1]]) ^ 2) / k)
        se.sigma <- sqrt(sum((muestra[2,] - PPSfit$estimate[[2]]) ^ 2) / k)
        se.nu <- sqrt(sum((muestra[3,] - PPSfit$estimate[[3]]) ^ 2) / k)
        se.value <- data.frame(lambda = se.lambda, sigma = se.sigma, nu = se.nu)
        
      }
      else {
        if (parallel == TRUE){
          cl <- makeCluster(ncores)
          registerDoParallel(cl)
          clusterExport(cl, list("PPSfit"), envir = environment())
          muestra <- foreach (i = 1 : k, .combine = cbind, .multicombine = TRUE, .inorder = FALSE, .packages = c("ParetoPosStable")) %dopar% {
            converged <- FALSE
            datos.sim <- sample(x = PPSfit$obs, size = PPSfit$n, replace = TRUE)
            ajuste <- PPS.fit(datos.sim, estim.method = PPSfit$estim.method, sigma = PPSfit$sigma, start = PPSfit$estimate)
            c(ajuste$estimate[[1]], ajuste$estimate[[2]])
          }
          stopCluster(cl)
        }
        else{
          muestra <- foreach (i = 1 : k, .combine = cbind, .multicombine = TRUE, .inorder = FALSE, .packages = c("ParetoPosStable")) %do% {
            converged <- FALSE
            datos.sim <- sample(x = PPSfit$obs, size = PPSfit$n, replace = TRUE)
            ajuste <- PPS.fit(datos.sim, estim.method = PPSfit$estim.method, sigma = PPSfit$sigma, start = PPSfit$estimate)
            c(ajuste$estimate[[1]], ajuste$estimate[[2]])
          }
        }
        se.lambda <- sqrt(sum((muestra[1,] - PPSfit$estimate[[1]]) ^ 2) / k)
        se.nu <- sqrt(sum((muestra[2,] - PPSfit$estimate[[2]]) ^ 2) / k)
        se.value <- data.frame(lambda = se.lambda, nu = se.nu)
      }
    }
    if (PPSfit$Pareto == TRUE){
      if (is.null(PPSfit$sigma)){
        if (parallel == TRUE){
          cl <- makeCluster(ncores)
          registerDoParallel(cl)
          clusterExport(cl, list("PPSfit"), envir = environment())
          muestra <- foreach (i = 1 : k, .combine = cbind, .multicombine = TRUE, .inorder = FALSE, .packages = c("ParetoPosStable")) %dopar% {
            converged <- FALSE
            datos.sim <- sample(x = PPSfit$obs, size = PPSfit$n, replace = TRUE)
            ajuste <- PPS.fit(datos.sim, estim.method = PPSfit$estim.method, start = PPSfit$estimate, Pareto = TRUE)
            c(ajuste$estimate[[1]], ajuste$estimate[[2]])
          }
          stopCluster(cl)
        }
        else{
          muestra <- foreach (i = 1 : k, .combine = cbind, .multicombine = TRUE, .inorder = FALSE, .packages = c("ParetoPosStable")) %do% {
            converged <- FALSE
            datos.sim <- sample(x = PPSfit$obs, size = PPSfit$n, replace = TRUE)
            ajuste <- PPS.fit(datos.sim, estim.method = PPSfit$estim.method, start = PPSfit$estimate, Pareto = TRUE)
            c(ajuste$estimate[[1]], ajuste$estimate[[2]])
          }
        }
        se.lambda <- sqrt(sum((muestra[1,] - PPSfit$estimate[[1]]) ^ 2) / k)
        se.sigma <- sqrt(sum((muestra[2,] - PPSfit$estimate[[2]]) ^ 2) / k)
        se.value <- data.frame(lambda = se.lambda, sigma = se.sigma)
      }
      else {
        if (parallel == TRUE){
          cl <- makeCluster(ncores)
          registerDoParallel(cl)
          clusterExport(cl, list("PPSfit"), envir = environment())
          muestra <- foreach (i = 1 : k, .combine = cbind, .multicombine = TRUE, .inorder = FALSE, .packages = c("ParetoPosStable")) %dopar% {
            converged <- FALSE
            datos.sim <- sample(x = PPSfit$obs, size = PPSfit$n, replace = TRUE)
            ajuste <- PPS.fit(datos.sim, estim.method = PPSfit$estim.method, start = PPSfit$estimate, sigma = PPSfit$sigma, Pareto = TRUE)
            c(ajuste$estimate[[1]])
          }
          stopCluster(cl)
        }
        else{
          muestra <- foreach (i = 1 : k, .combine = cbind, .multicombine = TRUE, .inorder = FALSE, .packages = c("ParetoPosStable")) %do% {
            converged <- FALSE
            datos.sim <- sample(x = PPSfit$obs, size = PPSfit$n, replace = TRUE)
            ajuste <- PPS.fit(datos.sim, estim.method = PPSfit$estim.method, start = PPSfit$estimate, sigma = PPSfit$sigma, Pareto = TRUE)
            c(ajuste$estimate[[1]])
          }
        }
        se.lambda <- sqrt(sum((muestra[1,] - PPSfit$estimate[[1]]) ^ 2) / k)
        se.value <- data.frame(lambda = se.lambda)
      }
    }
    ans <- list(se = se.value)
    ans$call <- match.call()
    class(ans) <- "se"
    ans
  }

#' @export
print.se <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  if(!inherits(x, "se")) stop("first argument should be of class 'se'")
  
  cat("\nStandard errors:\n")
  print(x$se)
  cat("\n")
}