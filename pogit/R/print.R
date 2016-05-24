#' Print an object of class \code{pogit}
#' 
#' The default print method for a \code{pogit} object.
#' 
#' Returns basic information about the model, the number of observations and 
#' covariates used, the number of regression effects subject to selection, 
#' MCMC options and the runtime used for the sampling algorithm. See
#' \code{\link{summary.pogit}} for more details. 
#'  
#' @param x an object of class \code{pogit}
#' @param ... further arguments passed to or from other methods (not used)
#'
#' @author Michaela Dvorzak <m.dvorzak@@gmx.at>
#' @export

print.pogit <- function(x, ...){
  
  stopifnot(class(x) == "pogit")
  if (x$family == "logit"){
    bin <- "binomial "
    if (all(x$data$N == 1)) bin <- ""
  }
  cat(paste(
    switch(as.character(x$BVS), 
           "TRUE"  = "Bayesian variable selection",
           "FALSE" = "MCMC"), 
    "for the", switch(as.character(x$family),
                      "logit"   = paste(bin, "logit", sep=""),
                      "pogit"   = switch(as.character(x$fun), 
                                         "select_poisson" = "Pogit",
                                         "select_poissonOD" = "overdispersed Pogit"),
                      "poisson" = switch(as.character(x$fun), 
                                         "select_poisson" = "Poisson",
                                         "select_poissonOD" = "overdispersed Poisson"),
                      "negbin"  = "negative binomial"), "model:\n"))
  
  cat("\nCall:\n")
  print(x$call)
  
  cat("\nModel:", length(x$data$y), "observations") 
  
  if (x$family %in% c("pogit", "poisson")){
    if (x$family == "pogit") cat("\n\n-Poisson:")
    cat("\n Covariates:", k <- x$model.pois$d)
    cat("\n --- subject to selection:", k - sum(x$model.pois$deltafix))
    rid <- as.logical(x$model.pois$ri)
    if (rid) rids <- as.logical(1 - x$model.pois$gammafix)
    cat("\n Random intercept included:", switch(as.character(rid),
                                                "FALSE" = "no",
                                                "TRUE" = "yes"))      
    if (rid) cat("\n --- subject to selection:",switch(as.character(rids),
                                                       "FALSE" = "no",
                                                       "TRUE" = "yes")) 
  }
  
  if (x$family %in% c("logit", "pogit")){
    if (x$family == "pogit") cat("\n\n-Logit:")
    cat("\n Covariates:", k <- x$model.logit$d)
    cat("\n --- subject to selection:", k - sum(x$model.logit$deltafix))
    rid <- as.logical(x$model.logit$ri)
    if (rid) rids <- as.logical(1 - x$model.logit$gammafix)
    cat("\n Random intercept included:", switch(as.character(rid),
                                                "FALSE" = "no",
                                                "TRUE" = "yes"))      
    if (rid) cat("\n --- subject to selection:", switch(as.character(rids),
                                                        "FALSE" = "no",
                                                        "TRUE" = "yes")) 
  }
  
  if (x$family == "pogit"){
    cat(paste("\n\nMethod:", 
              switch(x$method, 
                     "val"  = "validation data",
                     "infprior" = "informative prior\n")))
    if (x$method == "val"){
      cat("\nSample of validation data for", nrow(x$data$val), "categories\n")
    }
    if (x$method == "infprior"){
      pM <- as.character(paste("a0[", 0:(x$model.logit$d + x$model.logit$ri),"]", sep=""))
      priorSet2 <- round(with(x$prior.logit, c(priorMean = c(unname(m0), unname(aj0)))), 3)
      names(priorSet2) <- pM
      print(priorSet2)
    }
  }
  
  if (x$family == "negbin"){
    cat("\nCovariates:", k <- x$model.nb$d)
    cat("\n--- subject to selection:", k - sum(x$model.nb$deltafix))
  }
  
  cat("\n\nMCMC:")
  cat("\nM =", x$mcmc$M, "draws after a burn-in of", x$mcmc$burnin)
  
  if (x$BVS) cat("\nBVS started after", x$mcmc$startsel, "iterations")
  cat("\nThinning parameter:", x$mcmc$thin)
  
  if(x$family == "negbin"){
    cat(paste0("\n\nAcceptance rate for rho: ", x$acc.rho, "%"))
  }
  
  cat("\n\nRuntime:")
  cat("\ntotal:", round(x$dur$total, 1), "sec.")
  cat("\nsince burn-in:", round(x$dur$durM, 1), "sec.")
}
