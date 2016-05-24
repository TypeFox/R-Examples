#' @name GoF
#' @aliases GoF GoF.default
#' @title Goodness of fit tests for the Pareto Positive Stable (PPS) distribution
#' @description Kolmogorov-Smirnov, Anderson-Darling and PPS goodness of fit tests to validate a PPS fit (typically from \code{PPS.fit()}).
#' 
#' @param PPSfit A \code{PPSfit} Object.
#' @param k The number of iterations in the bootstrap procedure to approximate the p-values.
#' @param parallel A logical argument specifying if parallelization is allowed in the bootstrap iteration procedure.
#' @param ncores is the number of cores that we use if parallel is TRUE.
#' @param \dots Other arguments.
#'
#' @details
#'   It returns the Kolmogorov-Smirnov, the Anderson-Darling tests and a specific test for PPS distributions. p-values are approximated by a bootstrap procedure.   
#' The specific goodness of fit test for PPS distributions is based on the linearity of the survival function vs. the scaled observations in a double log-log scale (see Sarabia and Prieto, 2009).
#' @return A list with the values of the tests statistics and the approximated p-values.
#' @references Sarabia, J.M and Prieto, F. (2009). The Pareto-positive stable distribution: A new descriptive model for city size data, \emph{Physica A: Statistical Mechanics and its Applications}, \bold{388}(19), 4179-4191.
#' @seealso \code{\link{PPS.fit}}, \code{\link{plot.PPSfit}}
#' @examples
#' x <- rPPS(50, 1.2, 100, 2.3)
#' fit <- PPS.fit(x)
#' GoF(fit, k = 50)

#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster clusterExport detectCores
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom stats ecdf ks.test
#' @importFrom ADGofTest ad.test

#' @export
GoF <- function(PPSfit, k = 2000, parallel = TRUE, ncores = 2, ...) 
  UseMethod("GoF")

#' @export
GoF.default <-
  function (PPSfit, k = 2000, parallel = TRUE, ncores = 2, ...) {  
    if (is.null(PPSfit$Pareto)) PPSfit$Pareto <- FALSE
    if (PPSfit$Pareto == TRUE & !is.null(PPSfit$sigma)) pars <- c(as.numeric(PPSfit$estimate), PPSfit$sigma, 1)
    if (PPSfit$Pareto == TRUE & is.null(PPSfit$sigma)) pars <- c(as.numeric(PPSfit$estimate), 1)
    if (PPSfit$Pareto == FALSE & is.null(PPSfit$sigma)) pars <- as.numeric(PPSfit$estimate)
    if (PPSfit$Pareto == FALSE & !is.null(PPSfit$sigma)) pars <- c(PPSfit$estimate[[1]], PPSfit$sigma, PPSfit$estimate[[2]])
    datos <- PPSfit$obs
    estadistico.ks <- ks.test(datos, pPPS, lam = pars[1], sc = pars[2], v = pars[3])$statistic
    estadistico.ad <- ad.test(datos, pPPS, lam = pars[1], sc = pars[2], v = pars[3])$statistic
    rango <- function(y) - (length(y) * ecdf(y)(y) - (length(y) + 1))
    z <- log(log(datos[datos != pars[2]] / pars[2]))
    y <- log(-log(rango(datos[datos != pars[2]]) / (PPSfit$n + 1)))
    estadistico.pps <- sum((y - (log(pars[1]) + pars[3] * z)) ^ 2)
    metodo <- PPSfit$estim.method
    n <- PPSfit$n
    
    if (parallel == TRUE){
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      clusterExport(cl, list("PPSfit", "pars", "rPPS", "pPPS"), envir = environment())
      muestra <- foreach (i = 1 : k, .combine = cbind, .multicombine = TRUE, .inorder = FALSE, .packages = c("ParetoPosStable", "ADGofTest")) %dopar% {
        converged <- FALSE
        datos.sim <- rPPS(n, pars[1], pars[2], pars[3])
        if (PPSfit$Pareto == FALSE & is.null(PPSfit$sigma)){
          ajuste.sim <- PPS.fit(datos.sim, estim.method = metodo, Pareto = PPSfit$Pareto)
          pars.sim <- as.numeric(ajuste.sim$estimate)
        }
        if (PPSfit$Pareto == FALSE & !is.null(PPSfit$sigma)){
          ajuste.sim <- PPS.fit(datos.sim, estim.method = metodo, sigma = PPSfit$sigma, Pareto = PPSfit$Pareto)
          pars.sim <- c(ajuste.sim$estimate[[1]], PPSfit$sigma, ajuste.sim$estimate[[2]])   
        }
        if (PPSfit$Pareto == TRUE & is.null(PPSfit$sigma)){
          ajuste.sim <- PPS.fit(datos.sim, estim.method = metodo, Pareto = PPSfit$Pareto)
          pars.sim <- c(as.numeric(ajuste.sim$estimate), 1)
        }
        if (PPSfit$Pareto == TRUE & !is.null(PPSfit$sigma)){
          ajuste.sim <- PPS.fit(datos.sim, estim.method = metodo, sigma = PPSfit$sigma, Pareto = PPSfit$Pareto)
          pars.sim <- c(as.numeric(ajuste.sim$estimate), PPSfit$sigma, 1) 
        }
    test.ks<-ks.test(datos.sim, pPPS, lam = pars.sim[1], sc = pars.sim[2], v = pars.sim[3])
    ks.sample <- test.ks$statistic
    test.ad <- ad.test(datos.sim, pPPS, lam = pars.sim[1], sc = pars.sim[2], v = pars.sim[3])
    ad.sample <- test.ad$statistic
    z <- log(log(datos.sim[datos.sim != pars.sim[2]] / pars.sim[2]))
    y <- log(-log(rango(datos.sim[datos.sim != pars.sim[2]]) / (PPSfit$n + 1)))
    pps.sample <- sum((y - (log(pars.sim[1]) + pars.sim[3] * z)) ^ 2)
    c(ks.sample, ad.sample, pps.sample)
  }
  stopCluster(cl)
    } 
  else{
    muestra <- foreach (i = 1 : k, .combine = cbind, .multicombine = TRUE, .inorder = FALSE, .packages = c("ParetoPosStable", "ADGofTest")) %do% {
      converged <- FALSE
      datos.sim <- rPPS(n, pars[1], pars[2], pars[3])
      if (PPSfit$Pareto == FALSE & is.null(PPSfit$sigma)){
        ajuste.sim <- PPS.fit(datos.sim, estim.method = metodo, Pareto = PPSfit$Pareto)
        pars.sim <- as.numeric(ajuste.sim$estimate)
      }
      if (PPSfit$Pareto == FALSE & !is.null(PPSfit$sigma)){
        ajuste.sim <- PPS.fit(datos.sim, estim.method = metodo, sigma = PPSfit$sigma, Pareto = PPSfit$Pareto)
        pars.sim <- c(ajuste.sim$estimate[[1]], PPSfit$sigma, ajuste.sim$estimate[[2]])   
      }
      if (PPSfit$Pareto == TRUE & is.null(PPSfit$sigma)){
        ajuste.sim <- PPS.fit(datos.sim, estim.method = metodo, Pareto = PPSfit$Pareto)
        pars.sim <- c(as.numeric(ajuste.sim$estimate), 1)
      }
      if (PPSfit$Pareto == TRUE & !is.null(PPSfit$sigma)){
        ajuste.sim <- PPS.fit(datos.sim, estim.method = metodo, sigma = PPSfit$sigma, Pareto = PPSfit$Pareto)
        pars.sim <- c(as.numeric(ajuste.sim$estimate), PPSfit$sigma, 1) 
      }
      test.ks<-ks.test(datos.sim, pPPS, lam = pars.sim[1], sc = pars.sim[2], v = pars.sim[3])
      ks.sample <- test.ks$statistic
      test.ad <- ad.test(datos.sim, pPPS, lam = pars.sim[1], sc = pars.sim[2], v = pars.sim[3])
      ad.sample <- test.ad$statistic
      z <- log(log(datos.sim[datos.sim != pars.sim[2]] / pars.sim[2]))
      y <- log(-log(rango(datos.sim[datos.sim != pars.sim[2]]) / (PPSfit$n + 1)))
      pps.sample <- sum((y - (log(pars.sim[1]) + pars.sim[3] * z)) ^ 2)
      c(ks.sample, ad.sample, pps.sample)
    }
  }
  
    ks.sample <- muestra[1,]
    ad.sample <- muestra[2,]
    pps.sample <- muestra[3,]
    
    tests <- list(ks.statistic = estadistico.ks, ad.statistic = estadistico.ad, pps.statistic = estadistico.pps, ks.p.value = sum(estadistico.ks < ks.sample) / (k + 1), ad.p.value = sum(estadistico.ad < ad.sample) / (k + 1), pps.p.value = sum(estadistico.pps < pps.sample) / (k + 1), PPSfit)
    ans <- list(tests.results = tests)
    ans$call <- match.call()
    class(ans) <- "GoF"
    ans
  }

#' @export
print.GoF <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  if (!class(x) == "GoF") stop("Object must be returned by 'GoF' function") 
  cat("\nGoodness of fit tests:\n")
  tabla <- data.frame(Statistic.value = c(x$tests.results$ks.statistic, x$tests.results$ad.statistic, x$tests.results$pps.statistic), p.value = c(x$tests.results$ks.p.value, x$tests.results$ad.p.value, x$tests.results$pps.p.value))
  row.names(tabla) <- c("KS test", "AD test", "PPS test")
  print(tabla)
  cat("\n")
  invisible(x)
}
