#' Simulation to Estimate Statistical Power of a Rasch Model Test
#' 
#' This function conducts a simulation to estimate statistical power 
#' of a Rasch model test for user-specified item and person parameters.
#' 
#' The F-test in a three-way analysis of variance design   \eqn{(A \succ \mathbf{B}) x C}{}(A > \strong{B}) x C
#' with mixed classification (fixed factor A = subgroup, random factor B = testee, 
#' and fixed factor C = items) is used to simulate statistical power of a
#' Rasch model test. This approach using a F-distributed statistic, where 
#' the sample size directly affects the degree of freedom enables determination
#' of the sample size according to a given type I  and type II risk, and according 
#' to a certain effect of model misfit which is of practical relevance. 
#' Note, that this approach works as long as there exists no main effect of
#' A (subgroup). Otherwise an artificially high type I risk of the A x C interaction
#' F-test results - that is, the approach works as long as no statistically significant
#' main effect of A occurs.
#' 
#' @param b           Either a vector or an integer indicating the number of observations in each group.
#' @param ipar        Item parameters in both groups specified in a list.
#' @param ppar        Person parameters specified by a distribution for each group.
#' @param runs        Number of simulation runs. 
#' @param H0          If \code{TRUE}, null hypothesis condition is simulated.
#' @param sig.level   Nominal significance level.
#' @param method      Simulation method: for-loop or vectorized.
#' @param output      If \code{TRUE}, output is shown.
#' 
#' @author 
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at},
#' Jan Steinfeld \email{jan.steinfeld@@univie.ac.at}
#' 
#' @references
#' Kubinger, K. D., Rasch, D., & Yanagida, T. (2009). On designing data-sampling for Rasch model 
#' calibrating an achievement test. \emph{Psychology Science Quarterly, 51}, 370-384.
#'
#' Kubinger, K. D., Rasch, D., & Yanagida, T. (2011). A new approach for testing the Rasch model.
#' \emph{Educational Research and Evaluation, 17}, 321-333.
#' 
#' @return Returns a list with following entries:
#' 
#' \tabular{ll}{ 
#'   \code{b}         \tab number of observations in each group \cr
#'   \code{ipar}      \tab item parameters in both subgroups \cr
#'   \code{c}         \tab number of items \cr
#'   \code{ppar}      \tab distribution of person parameters \cr
#'   \code{runs}      \tab number of simulation runs \cr
#'   \code{sig.level} \tab nominal significance level \cr
#'   \code{H0.AC.p}   \tab \emph{p}-values of the interaction A x C in the null hypothesis condition (if \code{H0 = TRUE}) \cr
#'   \code{H1.AC.p}   \tab \emph{p}-values of the interaction A x C in the alternative hypothesis condition \cr
#'   \code{power}     \tab estimated statistical power \cr
#'   \code{type1}     \tab estimated significance level \cr
#' }
#' 
#' @seealso 
#' \code{\link{aov.rasch}}
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' # item parameters
#' ipar2 <- ipar1 <- seq(-3, 3, length.out = 20)
#' # model differential item function (DIF)
#' ipar2[10] <- ipar1[11]
#' ipar2[11] <- ipar1[10]
#' # simulation for b = 200 
#' pwr.rasch(200, ipar = list(ipar1, ipar2))
#' 
#' # simulation for b = 100, 200, 300, 400, 500 
#' pwr.rasch(seq(100, 500, by = 100), ipar = list(ipar1, ipar2))
#' 
#' # simulation for b = 100, 200, 300, 400, 500 
#' # uniform distribution [-3, 3] of person parameters
#' pwr.rasch(200, ipar = list(ipar1, ipar2), ppar = list("runif(b, -3, 3)", "runif(b, -3, 3)"))
#' }
pwr.rasch <- function(b, ipar = list(),
                      ppar = list("rnorm(b, mean = 0, sd = 1.5)", "rnorm(b, mean = 0, sd = 1.5)"),
                      runs = 1000, H0 = TRUE, sig.level = 0.05, 
                      method = c("loop", "vectorized"), output = TRUE) {
    
  #--------------------------------------------------------------------------------------------------------#
  # Input Check
  l.call <- match.call()
  
  # sum of item parameters = 0
  if (all(round(unlist(lapply(ipar, sum)), 3) != 0)) {
    
     stop("Item pararameters are not normalized to sum-0")
    
  } 
  
  # number of items
  if (length(unique(unlist(lapply(ipar, length)))) != 1) {
    
     stop("Different number of items specified in both groups")
    
  }   
  
  #--------------------------------------------------------------------------------------------------------#

  cat("--------------------------------------------------------\n")
  cat("  Call:    "); print(l.call)
  cat("  Time:   ", time <- paste(Sys.time()), "\n")
  cat("  R:      ", R.version$version.string, "\n")
  cat("  Package:", pkg.version<- paste0("pwrRasch version ", packageDescription("pwrRasch")$Version,
                                         " (", packageDescription("pwrRasch")$Date, ")"), "\n")
  cat("--------------------------------------------------------\n") 
 
  if (length(b) == 1) {

      simres <- pwr.rasch.internal(b = b, ipar = ipar, ppar = ppar, 
                                  runs = runs, H0 = H0, sig.level = sig.level, 
                                  method = method, output = output)

  } else { # length(b) > 1
      
    j <- 1 
    simres <- NULL 
    for (i in b) {
              
      if (i == b[1]) {
          
        simres[[j]] <- pwr.rasch.internal(b = i, ipar = ipar, ppar = ppar, 
                                          runs = runs, H0 = H0, sig.level = sig.level, 
                                          method = method, output = FALSE)
        j <- j + 1 
          
      } else {
          
        simres[[j]] <- pwr.rasch.internal(b = i, ipar = ipar, ppar = ppar, 
                                          runs = runs, H0 = H0, sig.level = sig.level, 
                                          method = method, output = FALSE)     
        j <- j + 1 
      
      }
        
    }
      
  }  
    
  #------------------------------------------------------------------------------------------------------#
  # Output
           
  if (length(b) == 1) {
      
    if (output == TRUE) {

      c <- unique(unlist(lapply(ipar, length)))
        
      if (H0 == TRUE) {
          
        cat("\n  Statistical Power Simulation for the Rasch model \n\n",
              
            "    b (number of persons in each group): ", b, "\n",
            "    c (number of items):                 ", c, "\n",        
            "    simulation runs:                     ", runs, "\n\n",
              
            "    Estimated statistical power: ", formatC(simres$power, format = "f", digits = 3), "\n",
            "    Nominal significance level:  ", sig.level, "\n",
            "    Empirical significance level:", formatC(simres$type1, format = "f", digits = 3), "\n")
                        
      } else {
          
        cat("\n  Statistical Power Simulation for the Rasch model \n\n", 
              
            "    b (numer of persons in each group): ", b, "\n",
            "    c (numer of items):                 ", c, "\n",   
            "    simulation runs:                    ", runs, "\n\n",
              
            "    Estimated statistical power: ", formatC(simres$power, format = "f", digits = 3), "\n",
            "    Nominal significance level:  ", sig.level, "\n")  
                
      }
        
      cat("--------------------------------------------------------\n")
        
    }

    #------------------------------------------------------------------------------------------------------#
    # Warnings
      
    if (H0 == TRUE) {
        
      # Warning: Nominal vs. empirical significane level
      if (round(sum(simres$H0.AC.p < sig.level) / runs, 3) < (sig.level - sig.level / 100 * 20) |
          round(sum(simres$H0.AC.p < sig.level) / runs, 3) > (sig.level + sig.level / 100 * 20)) {
          
        warning("F-test does not hold its type-I-risk with 20%-robusteness, i.e, results may not be trustworthy.")
        
      }
        
    } 
          
  } else {  # length(b) > 1  
    
    if (output == TRUE) {
          
      cat("\n  Statistical Power Simulation for the Rasch model \n\n",
              
          "    b (numer of persons in each group): ", paste(b, collapse = ", "), "\n",
          "    c (numer of items):                 ", unique(unlist(lapply(ipar, length))), "\n",
          "    simulation runs:                    ", runs, "\n\n",
              
          "    Estimated statistical power: \n")
          
      for (i in 1:length(simres)) {
            
        cat(paste0("      b = ", formatC(b[i], digits = max(nchar(b)) - 1, format = "d"), ":"),
            formatC(simres[[i]][["power"]], format = "f", digits = 3), "\n")
            
      }  
          
          
      cat("\n     Nominal significance level:  ", sig.level, "\n")
          
      if (H0 == TRUE) {
             
        cat("\n     Empirical significance level: \n")    
            
        for (i in 1:length(simres)) {
              
          cat(paste0("      b = ", formatC(b[i], digits = max(nchar(b)) - 1, format = "d"), ":"),
              formatC(simres[[i]][["type1"]], format = "f", digits = 3), "\n")
              
        }  
            
      }
          
      cat("------------------------------------------------------\n")
          
    }

    #------------------------------------------------------------------------------------------------------#
    # Warnings
    
    if (H0 == TRUE) {
      
      # Warning: Nominal vs. empirical significane level
      sig.level.act <- unlist(lapply(simres, function(x) sum(x$H0.AC.p < sig.level) / runs))
      if (any(sig.level.act < (sig.level - sig.level / 100 * 20)) |
          any(sig.level.act > (sig.level + sig.level / 100 * 20))) {
            
        warning("F-test does not hold its type-I-risk with 20%-robusteness, i.e, results may not be trustworthy.")
            
      }
      
    }         
            
  } 
    
  class(simres) <- "pwrrasch"
  return(invisible(simres))
    
} 