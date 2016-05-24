#' Convert a virtual species distribution (or a suitability raster) into presence-absence
#' 
#' This functions converts the probabilities of presence from the output of
#'  \code{\link{generateSpFromFun}}, \code{\link{generateSpFromPCA}}, \code{\link{generateRandomSp}}
#' or a suitability raster into
#' a presence-absence raster. The conversion can be threshold-based, or based
#' on a probability of conversion (see details).
#' 
#' @param x a suitability raster, or the output from functions 
#' \code{\link{generateSpFromFun}}, \code{\link{generateSpFromPCA}} 
#' or \code{\link{generateRandomSp}}
#' @param PA.method \code{"threshold"} or \code{"probability"}. If 
#' \code{"threshold"}, then occurrence probabilities are simply converted into
#' presence-absence according to the threshold \code{beta}. If \code{"probability"}, then
#' probabilities are converted according to a logistic function of threshold 
#' \code{beta} and slope \code{alpha}.
#' @param beta \code{"random"}, a numeric value in the range of your 
#' probabilities or \code{NULL}. This is the threshold of conversion into
#' presence-absence (= the inflexion point if \code{PA.method = "probability"}).
#' If \code{"random"}, a numeric value will be randomly generated within the range
#' of \code{x}.
#' @param alpha \code{NULL} or a negative numeric value. Only useful if 
#' \code{PA.method = "probability"}. The value of \code{alpha} will
#' shape the logistic function transforming occurrences into presence-absences.
#' See \code{\link{logisticFun}} and examples therein for the choice of \code{alpha}
#' @param species.prevalence \code{NULL} or a numeric value between 0 and 1.
#' The species prevalence is the proportion of sites actually occupied by the
#' species.
#' @param plot \code{TRUE} or \code{FALSE}. If \code{TRUE}, maps of probabilities
#' of occurrence and presence-absence will be plotted.
#' 
#' @export
#' @import raster
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, C. Bellard & F. Courchamp
#' @references
#' Meynard C.N. & Kaplan D.M. 2013. Using virtual species to study species 
#' distributions and model performance. 
#' \emph{Journal of Biogeography} \bold{40}:1-8
#' 
#' Meynard C.N. & Kaplan D.M. 2011. The effect of a gradual response to the 
#' environment on species distribution model performance.
#' \emph{Ecography} \bold{35}:499-509
#' 
#' @details 
#' The conversion of probabilities of occurrence into presence - absence is
#' usually performed by selecting a threshold above which presence always occurs,
#' and never below. However, this approach may be too much unrealistic because
#' species may sometime be present in areas with a low probability of occurrence,
#' or be absent from areas with a high probability of occurrence. In addition,
#' when using a threshold you erase the previously generated response shapes: 
#' it all becomes threshold. Thus, this threshold approach should be avoided.
#' 
#'  
#' Hence, a more
#' realistic conversion consists in converting probabilities into presence -
#' absence with a probability function (see references). Such a probability 
#' conversion can be performed here with a logit function 
#' (see \code{\link{logisticFun}}).
#' 
#' To perform the probability conversion you have to define two of the
#' three following parameters:
#' \itemize{
#' \item{\code{beta}: the 'threshold' of the logistic function (i.e. the 
#' inflexion point)}
#' \item{\code{alpha}: the slope of the logistic function}
#' \item{\code{species.prevalence}: the proportion of sites in which the species
#' occur}
#' }
#' 
#' If you provide \code{beta} and \code{alpha}, the \code{species.prevalence}
#' is calculated immediately calculated after conversion into presence-absence.
#' 
#' On the other hand, if you provide \code{species.prevalence} and either
#' \code{beta} or \code{alpha}, the function will try to determine \code{alpha}
#' (if you provided \code{beta}) or \code{beta} (if you provided \code{alpha}).
#' 
#' The relationship between species prevalence, alpha and beta is dependent
#' on the available range of environmental conditions (see Meynard and Kaplan,
#' 2011 and especially the Supporting Information). As a consequence, the 
#' desired species prevalence may not be available for the defined \code{alpha} 
#' or \code{beta}. In these conditions, the function will retain the \code{alpha} or
#' \code{beta} which provides the closest prevalence to your \code{species.prevalence},
#' but you may also provide another value of \code{alpha} or \code{beta} to obtain
#' a closer prevalence. 
#'  
#' In all cases, the \code{species.prevalence} indicated in the output is the
#' prevalence measured on the output presence-absence map.
#' @note
#' The approximation of \code{alpha} or \code{beta} to the chosen 
#' \code{species.prevalence} may take time if you work on very large rasters.
#' @return
#' a \code{list} containing 5 elements:
#' \itemize{
#' \item{\code{approach}: the approach used to generate the species, \emph{i.e.}, \code{"response"}}
#' \item{\code{details}: the details and parameters used to generate the species}
#' \item{\code{suitab.raster}: the virtual species distribution, as a Raster object containing the
#' environmental suitability)}
#' \item{\code{PA.conversion}: the parameters used to convert the suitability into presence-absence}
#' \item{\code{pa.raster}: the presence-absence map, as a Raster object containing 0 (absence) / 1 (presence) / NA}
#' }
#' The structure of the virtualspecies object can be seen using str()
#' @examples
#' # Create an example stack with two environmental variables
#' a <- matrix(rep(dnorm(1:100, 50, sd = 25)), 
#'             nrow = 100, ncol = 100, byrow = TRUE)
#' env <- stack(raster(a * dnorm(1:100, 50, sd = 25)),
#'              raster(a * 1:100))
#' names(env) <- c("variable1", "variable2")
#' 
#' # Creation of the parameter list
#' parameters <- formatFunctions(variable1 = c(fun = 'dnorm', mean = 0.00012,
#'                                             sd = 0.0001),
#'                               variable2 = c(fun = 'linearFun', a = 1, b = 0))
#' sp1 <- generateSpFromFun(env, parameters, plot = FALSE)
#' 
#' # Conversion into presence-absence with a threshold-based approach
#' convertToPA(sp1, PA.method = "threshold", beta = 0.2,  plot = TRUE)
#' convertToPA(sp1, PA.method = "threshold", beta = "random", plot = TRUE)
#' 
#' # Conversion into presence-absence with a probability-based approach
#' convertToPA(sp1, PA.method = "probability", beta = 0.4, 
#'               alpha = -0.05, plot = TRUE)
#' convertToPA(sp1, PA.method = "probability", beta = "random", 
#'               alpha = -0.1, plot = TRUE)
#'               
#' # Conversion into presence-absence by choosing the prevalence
#' # Threshold method
#' convertToPA(sp1, PA.method = "threshold",
#'               species.prevalence = 0.3, plot = TRUE)
#' # Probability method, with alpha provided              
#' convertToPA(sp1, PA.method = "probability", alpha = -0.1, 
#'               species.prevalence = 0.2, plot = TRUE)        
#' # Probability method, with beta provided              
#' convertToPA(sp1, PA.method = "probability", beta = 0.5, 
#'               species.prevalence = 0.2, alpha = NULL, 
#'               plot = TRUE)    
#'  
#' # Plot the output Presence-Absence raster only             
#' sp1 <- convertToPA(sp1, plot = FALSE)
#' plot(sp1$pa.raster)
                    

convertToPA <- function(x, 
                        PA.method = "probability",
                        beta = "random",
                        alpha = -.05,
                        species.prevalence = NULL,
                        plot = TRUE)

{
  if("virtualspecies" %in% class(x))
  {
    if(class(x$suitab.raster) == "RasterLayer")
    {
      sp.raster <- x$suitab.raster
    } else stop("x must be:\n- a raster layer object\nor\n- the output list from functions
               generateSpFromFun(), generateSpFromPCA() or generateRandomSp()")
  } else if (is.raster(x))
  {
    sp.raster <- x
  } else stop("x must be:\n- a raster layer object\nor\n- the output list from functions
               generateSpFromFun(), generateSpFromPCA() or generateRandomSp()")
  
  if(any(is.na(maxValue(sp.raster))))
  {
    sp.raster <- setMinMax(sp.raster)
  }
  
  if(PA.method == "threshold")
  {
    if(is.numeric(beta))
    {
      if(is.numeric(species.prevalence))
      {
        warning("Both beta and species.prevalence were provided. beta will be
                ignored.")
        beta <- NULL
      } else if(beta < sp.raster@data@min | beta > sp.raster@data@max)
      {
        stop("beta must be in the range of your suitability raster")
      }
    } else if(beta == "random")
    {
      if(is.numeric(species.prevalence))
      {
        beta <- NULL 
      } else
      {
        beta <- sample(seq(sp.raster@data@min, 
                           sp.raster@data@max, length = 1000), 1)
      }
    } else if(is.null(beta))
    {
      if(is.null(species.prevalence))
      {
        stop("Either provide beta or species.prevalence when choosing
             PA.method = 'threshold'")
      }
      } else
      {
        stop("beta must either be 'random', a numeric value within the range of
             your data or NULL")
      }
      } else if(PA.method == "probability")
      {
        if(length(c(alpha, beta, species.prevalence)) < 1)
        {
          if(!is.null(species.prevalence))
          {
            warning("Neither alpha nor beta were provided. As a consequence, alpha
                    will be determined to a random value, and beta will be adjusted 
                    automatically to the desired species prevalence.")
            alpha <- -sample(c(seq((sp.raster@data@max - sp.raster@data@min)/1000,
                                   (sp.raster@data@max - sp.raster@data@min)/100, length = 10),
                               seq((sp.raster@data@max - sp.raster@data@min)/100,
                                   (sp.raster@data@max - sp.raster@data@min)/10, length = 100),
                               seq((sp.raster@data@max - sp.raster@data@min)/10,
                                   (sp.raster@data@max - sp.raster@data@min)*10, length = 10)), size = 1)
          } else
          {
            stop("If you choose PA.method = 'probability', you must provide two of the
                 three following parameters: beta, alpha and species.prevalence.")
          }
          } else if(length(c(alpha, beta, species.prevalence)) > 2)
          {
            if(beta != "random")
            {
              warning("You should not provide the three parameters beta, alpha and 
                      species.prevalence. beta will be ignored and automatically defined by the
                      function")
            }
            beta <- NULL
            } 
        # Checking the arguments. At this stage only two of them should be not NULL
        if(!is.null(beta))
        {
          if(is.numeric(beta))
          {
            if(beta < sp.raster@data@min | beta > sp.raster@data@max)
            {
              stop("beta must be in the range of your suitability raster")
            }
          } else if(beta == "random")
          {
            beta <- sample(seq(sp.raster@data@min, 
                               sp.raster@data@max, length = 1000), 1)
          } else
          {
            stop("beta must either be 'random', a numeric value within the range of
                 your data or NULL")
          }
          }
        
        if(!is.null(species.prevalence))
        {
          if(is.numeric(species.prevalence))
          {
            if(!(species.prevalence >= 0 & 
                   species.prevalence <= 1))
            {
              stop("species.prevalence must be a numeric value between 0 and 1.")
            }
          } else 
          {
            stop("species.prevalence must either be a numeric value between 0 and 1
                 or NULL")
          }
          }
        
        if(!is.null(alpha))
        {
          if(!is.numeric(alpha))
          {
            stop("Please provide a numeric value to alpha")
          } else if(alpha > 0)
          {
            warning("alpha was provided > 0. 
                    This means that low probabilities will be converted to presences, 
                    and high probabilities will be converted to absences.
                    If this is not what was intended, provide a negative alpha.")
          }
          }
        
        if(!is.null(species.prevalence))
        {
          if(!is.null(beta))
          {
            message("   --- Determing alpha automatically according to beta and species.prevalence\n\n")
          } else
          {
            message("   --- Determing beta automatically according to alpha and species.prevalence\n\n")
          }
        } else
        {
          message("   --- Determing species.prevalence automatically according to alpha and beta\n\n")
        }
          }
  
  if (PA.method == "probability")
  {
    if(!is.null(species.prevalence))
    {
      if(!is.null(beta))
      {
        alpha.test <- NULL
        for (alpha in c((sp.raster@data@max - sp.raster@data@min)/1000, (sp.raster@data@max - sp.raster@data@min) * 10))
        {
          if(alpha > 0) alpha <- -alpha
          PA.raster <- calc(sp.raster, fun = function(x)
          {
            logisticFun(x, beta = beta, alpha = alpha)
          })
          
          PA.raster <- calc(PA.raster, fun = function(x)
          {
            sapply(x, FUN = function(y)
            {
              if(is.na(y))
              { NA } else
              {
                sample(x = c(0, 1),  size = 1, prob = c(1 - y, y))
              }
            }
            )
          }
          )
          freqs <- freq(PA.raster)
          if(any(is.na(freqs[, 1])))
          {
            freqs <- freqs[-which(is.na(freqs[, 1])), ]
          }
          alpha.test <- rbind(alpha.test, c(alpha, 
                                            ifelse(nrow(freqs) == 2,
                                                   freqs[freqs[, "value"] == 1, "count"] / sum(freqs[, "count"]),
                                                   ifelse(freqs[, "value"] == 0,
                                                          0, 1))))
        }
        epsilon <- species.prevalence - alpha.test[, 2]
        if(all(epsilon > 0))
        {
          warning(paste("Warning, the desired species prevalence cannot be obtained, because of the chosen beta and available environmental conditions (see details).
                        The closest possible estimate of prevalence was", round(alpha.test[2, 2], 2),
                        "\nPerhaps you can try a lower beta value."))
          alpha <- alpha.test[2, 1]
        } else if (all(epsilon < 0))
        {
          warning(paste("Warning, the desired species prevalence cannot be obtained, because of the chosen beta and available environmental conditions (see details).
                        The closest possible estimate of prevalence was", round(alpha.test[1, 2], 2),
                        "\nPerhaps you can try a higher beta value."))
          alpha <- alpha.test[1, 1]
        } else 
        {
          while (all(abs(epsilon) > 0.01))
          {
            alpha <- (alpha.test[which(epsilon == max(epsilon[epsilon < 0])), 1] + 
                        alpha.test[which(epsilon == min(epsilon[epsilon > 0])), 1]) / 2
            PA.raster <- calc(sp.raster, fun = function(x)
            {
              logisticFun(x, beta = beta, alpha = alpha)
            })
            
            PA.raster <- calc(PA.raster, fun = function(x)
            {
              sapply(x, FUN = function(y)
              {
                if(is.na(y))
                { NA } else
                {
                  sample(x = c(0, 1),  size = 1, prob = c(1 - y, y))
                }
              }
              )
            }
            )
            freqs <- freq(PA.raster)
            if(any(is.na(freqs[, 1])))
            {
              freqs <- freqs[-which(is.na(freqs[, 1])), ]
            }
            alpha.test <- rbind(alpha.test, c(alpha, 
                                              ifelse(nrow(freqs) == 2,
                                                     freqs[freqs[, "value"] == 1, "count"] / sum(freqs[, "count"]),
                                                     ifelse(freqs[, "value"] == 0,
                                                            0, 1))))
            
            epsilon <- species.prevalence - alpha.test[, 2]
          }
        }
      } else
      {
        beta.test <- NULL
        for (beta in c(sp.raster@data@min, sp.raster@data@max))
        {
          PA.raster <- calc(sp.raster, fun = function(x)
          {
            logisticFun(x, beta = beta, alpha = alpha)
          })
          
          PA.raster <- calc(PA.raster, fun = function(x)
          {
            sapply(x, FUN = function(y)
            {
              if(is.na(y))
              { NA } else
              {
                sample(x = c(0, 1),  size = 1, prob = c(1 - y, y))
              }
            }
            )
          }
          )
          
          freqs <- freq(PA.raster)
          if(any(is.na(freqs[, 1])))
          {
            freqs <- freqs[-which(is.na(freqs[, 1])), ]
          }
          beta.test <- rbind(beta.test, c(beta, 
                                          ifelse(nrow(freqs) == 2,
                                                 freqs[freqs[, "value"] == 1, "count"] / sum(freqs[, "count"]),
                                                 ifelse(freqs[, "value"] == 0,
                                                        0, 1))))
        }
        epsilon <- species.prevalence - beta.test[, 2]
        if(all(epsilon > 0))
        {
          warning(paste("Warning, the desired species prevalence cannot be obtained, because of the chosen alpha and available environmental conditions (see details).
                        The closest possible estimate of prevalence was", round(beta.test[1, 2], 2),
                        "\nPerhaps you can try an alpha value closer to 0."))
          beta <- beta.test[1, 1]
        } else if (all(epsilon < 0))
        {
          warning(paste("Warning, the desired species prevalence cannot be obtained, because of the chosen beta and available environmental conditions (see details).
                        The closest possible estimate of prevalence was", round(beta.test[1, 2], 2),
                        "\nPerhaps you can try an alpha value closer to 0."))
          beta <- beta.test[2, 1]
        } else 
        {
          while (all(abs(epsilon) > 0.01))
          {
            beta <- (beta.test[which(epsilon == max(epsilon[epsilon < 0])), 1] + 
                       beta.test[which(epsilon == min(epsilon[epsilon > 0])), 1]) / 2
            PA.raster <- calc(sp.raster, fun = function(x)
            {
              logisticFun(x, beta = beta, alpha = alpha)
            })
            
            PA.raster <- calc(PA.raster, fun = function(x)
            {
              sapply(x, FUN = function(y)
              {
                if(is.na(y))
                { NA } else
                {
                  sample(x = c(0, 1),  size = 1, prob = c(1 - y, y))
                }
              }
              )
            }
            )
            freqs <- freq(PA.raster)
            if(any(is.na(freqs[, 1])))
            {
              freqs <- freqs[-which(is.na(freqs[, 1])), ]
            }
            beta.test <- rbind(beta.test, c(beta, 
                                            ifelse(nrow(freqs) == 2,
                                                   freqs[freqs[, "value"] == 1, "count"] / sum(freqs[, "count"]),
                                                   ifelse(freqs[, "value"] == 0,
                                                          0, 1))))
            epsilon <- species.prevalence - beta.test[, 2]
          }
        }
      }
    } 
    
    PA.raster <- calc(sp.raster, fun = function(x)
    {
      logisticFun(x, beta = beta, alpha = alpha)
    })
    PA.raster <- calc(PA.raster, fun = function(x)
    {
      sapply(x, FUN = function(y)
      {
        if(is.na(y))
        { NA } else
        {
          sample(x = c(0, 1),  size = 1, prob = c(1 - y, y))
        }
      }
      )
    }
    )
    
      } else if (PA.method == "threshold")
      {
        if(!is.null(species.prevalence))
        {
          beta <- quantile(sp.raster, 1 - species.prevalence)
          names(beta) <- NULL
        }
        message("    - Threshold of conversion into PA:", round(beta, 3), "; method = ", PA.method, "\n")
        
        PA.raster <- reclassify(sp.raster,
                                c(-Inf, beta, 0,
                                  beta, +Inf, 1))
      } else {stop("Wrong PA.method entered (either 'probability' or 'threshold')")}
  
  species.prevalence <- round(freq(PA.raster)[2, 2] / sum(freq(PA.raster)[1:2, 2]), 3)
  names(species.prevalence) <- NULL
  
  if("virtualspecies" %in% class(x))
  {
    if(PA.method == "threshold")
    {
      x$PA.conversion = c(conversion.method = PA.method,
                          cutoff = beta,
                          species.prevalence = species.prevalence)
    } else
    {
      x$PA.conversion = c(conversion.method = PA.method,
                          alpha = alpha,
                          beta = beta,
                          species.prevalence = species.prevalence)
    }
    x$pa.raster = PA.raster
    results <- x
    if(plot) plot(stack(results$suitab.raster, results$pa.raster), main = c("Suitability", "Presence-absence"))
  } else if (is.raster(x))
  {
    if(PA.method == "threshold")
    {
      PA.conversion = c(cutoff = beta,
                        conversion.method = PA.method, 
                        species.prevalence = species.prevalence)
    } else
    {
      PA.conversion = c(conversion.method = PA.method, 
                        alpha = alpha,
                        beta = beta,
                        species.prevalence = species.prevalence)
    }
    results <- list(suitab.raster = x,
                    PA.conversion = PA.conversion,
                    pa.raster = PA.raster)
    if(plot) plot(stack(results$suitab.raster, results$pa.raster), main = c("Suitability", "Presence-absence"))
  }
  return(results)
}




