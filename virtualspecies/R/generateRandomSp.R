#' Generate a random virtual species distribution from environmental variables
#' 
#' @description
#' This function generates randomly a virtual species distribution.
#' 
#' @param raster.stack a RasterStack object, in which each layer represent an environmental 
#' variable.
#' @param approach \code{"automatic"}, \code{"random"}, \code{"response"}
#' or \code{"pca"}. This parameters defines how species will be generated. 
#' \code{"automatic"}: If less than 6 variables in \code{raster.stack}, a 
#' response approach will be used, otherwise a pca approach will be used.
#' \code{"random"}: the approach will be randomly picked. Otherwise choose
#' \code{"response"} or \code{"pca"}. See details.
#' @param rescale \code{TRUE} or \code{FALSE}. If \code{TRUE}, the final 
#' probability of presence is rescaled between 0 and 1.
#' @param convert.to.PA \code{TRUE} or \code{FALSE}. If \code{TRUE}, the 
#' virtual species distribution will also be converted into Presence-Absence.
#' @param relations [response approach] a vector containing the possible types of response function.
#' The implemented type of relations are \code{"gaussian"}, \code{"linear"},
#' \code{"logistic"} and \code{"quadratic"}.
#' @param rescale.each.response \code{TRUE} or \code{FALSE}. If \code{TRUE}, the individual responses to
#' each environmental variable are rescaled between 0 and 1
#' @param realistic.sp [response approach] \code{TRUE} or \code{FALSE}. If 
#' \code{TRUE}, the function will try to define responses that can form a viable
#' species. If \code{FALSE}, the responses will be randomly generated
#' (may result in environmental conditions that do not co-exist).
#' @param species.type [response approach] \code{"additive"} or \code{"multiplicative"}. Defines 
#' how the final probability of presence is calculated: if \code{"additive"}, responses to each
#' variable are summed; if \code{"multiplicative"}, responses are multiplicated.
#' See \code{\link{generateSpFromFun}}
#' @param niche.breadth [pca approach] \code{"any"}, \code{"narrow"} or \code{"wide"}. This parameter
#' defines how tolerant is the species regarding environmental conditions by adjusting
#' the standard deviations of the gaussian functions. See \code{\link{generateSpFromPCA}}
#' @param sample.points [pca approach] \code{TRUE} of \code{FALSE}. If you have a large
#' raster file then use this parameter to sample a number of points equal to
#' \code{nb.points}.
#' @param nb.points [pca approach] a numeric value. Only useful if \code{sample.points = TRUE}.
#' The number of sampled points from the raster, to perform the PCA. A too small
#' value may not be representative of the environmental conditions in your raster.
#' @param PA.method \code{"threshold"} or \code{"probability"}. If 
#' \code{"threshold"}, then occurrence probabilities are simply converted into
#' presence-absence according to the threshold \code{beta}. If \code{"probability"}, then
#' probabilities are converted according to a logistic function of threshold 
#' \code{beta} and slope \code{alpha}.
#' @param beta \code{"random"}, a numeric value in the range of your 
#' probabilities or \code{NULL}. This is the threshold of conversion into
#' presence-absence (= the inflexion point if \code{PA.method = "probability"}).
#' If \code{"random"}, a numeric value will be randomly generated within the range
#' of probabilities of occurrence. See \code{\link{convertToPA}}
#' @param alpha \code{NULL} or a negative numeric value. Only useful if 
#' \code{PA.method = "probability"}. The value of \code{alpha} will
#' shape the logistic function transforming occurrences into presence-absences.
#' See \code{\link{logisticFun}} and examples therein for the choice of \code{alpha}
#' @param species.prevalence \code{NULL} or a numeric value between 0 and 1.
#' The species prevalence is the proportion of sites actually occupied by the
#' species. See \code{\link{convertToPA}}
#' @param plot \code{TRUE} or \code{FALSE}. If \code{TRUE}, the generated virtual species will be plotted.
#' @details
#' This function generate random virtual species, either using a PCA approach, or using
#' a response approach. In case of a response approach, only four response functions are
#' currently used: gaussian, linear, logistic and quadratic functions.
#' 
#' Note that in case of numerouse predictor variables, the "response" approach will
#' not work well because it will often generate contradicting response functions 
#' (e.g., mean annual temperature optimum at degrees C and temperature of the coldest month at
#' 10 degrees C). In these case, it is advised to use the PCA approach (by default, a PCA approach
#' will be used if there are more than 6 predictor variables).
#' 
#' If \code{rescale.each.response = TRUE}, then the probability response to each
#' variable will be normalised between 0 and 1 according to the following formula:
#' P.rescaled = (P - min(P)) / (max(P) - min (P)). Simlarly, if \code{rescale = TRUE},
#' the final environmental suitability will be rescaled between 0 and 1 with the same formula.
#' 
#' By default, the function will perform a probabilistic conversion into presence-
#' absence, with a randomly chosen beta threshold. If you want to custmose the 
#' conversion parameters, you have to define \bold{two} of the three following parameters:
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
#' As explained in \code{\link{convertToPA}}, if you choose choose a precise
#' \code{species.prevalence}, it may not be possible to reach this particular 
#' value because of the availability of environmental conditions. Several
#' runs may be necessary to reach the desired \code{species.prevalence}.
#' @export
#' @import raster
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, C. Bellard & F. Courchamp
#' @return a \code{list} with 3 to 5 elements (depending if the conversion to presence-absence was performed):
#' \itemize{
#' \item{\code{approach}: the approach used to generate the species, \emph{i.e.}, \code{"response"}}
#' \item{\code{details}: the details and parameters used to generate the species}
#' \item{\code{suitab.raster}: the virtual species distribution, as a Raster object containing the
#' environmental suitability)}
#' \item{\code{PA.conversion}: the parameters used to convert the suitability into presence-absence}
#' \item{\code{pa.raster}: the presence-absence map, as a Raster object containing 0 (absence) / 1 (presence) / NA}
#' }
#' The structure of the virtualspecies can object be seen using str()
#' @examples
#' # Create an example stack with six environmental variables
#' a <- matrix(rep(dnorm(1:100, 50, sd = 25)), 
#'             nrow = 100, ncol = 100, byrow = TRUE)
#' env <- stack(raster(a * dnorm(1:100, 50, sd = 25)),
#'              raster(a * 1:100),
#'              raster(a * logisticFun(1:100, alpha = 10, beta = 70)),
#'              raster(t(a)),
#'              raster(exp(a)),
#'              raster(log(a)))
#' names(env) <- paste("Var", 1:6, sep = "")   
#' 
#' # More than 6 variables: by default a PCA approach will be used
#' generateRandomSp(env)
#' 
#' # Manually choosing a response approach
#' generateRandomSp(env, approach = "response")
#' 
#' # Randomly choosing the approach
#' generateRandomSp(env, approach = "random")
#' 
#' 

generateRandomSp <- function(raster.stack, 
                                  approach = "automatic",
                                  rescale = TRUE,
                                  convert.to.PA = TRUE,
                                  relations = c("gaussian", "linear", "logistic", "quadratic"),
                                  rescale.each.response = TRUE,
                                  realistic.sp = TRUE,
                                  species.type = "multiplicative",
                                  niche.breadth = "any",
                                  sample.points = FALSE, 
                                  nb.points = 10000,
                                  PA.method = "probability",
                                  alpha = -.1,
                                  beta = "random",
                                  species.prevalence = NULL,
                                  plot = TRUE)
{
  if(!(is(raster.stack, "Raster")))
  {
    stop("raster.stack must be a raster stack object")
  }
  if(any(is.na(maxValue(raster.stack))))
  {
    raster.stack <- setMinMax(raster.stack)
  }
  
  if(approach == "automatic")
  {
    if(nlayers(raster.stack) <= 5)
    {
      approach <- "response"
    } else
    {
      approach <- "pca"
    }
  } else if (approach == "random")
  {
    approach <- sample(c("response", "pca"), 1)
  } else if(!(approach %in% c("response", "pca")))
  {
    stop("Argument approach was misspecified. Either choose 'automatic', 'random', 'response' or 'pca'.")
  }

  var.names <- names(raster.stack)
  
  if(approach == "pca")
  {
    results <- generateSpFromPCA(raster.stack,
                          niche.breadth = niche.breadth,
                          sample.points = sample.points, 
                          nb.points = nb.points,
                          plot = FALSE)
  } else if (approach == "response")
  {
    parameters <- list()
    message(" - Determining species' response to predictor variables\n")
    if(any(!(relations %in% c("gaussian", "linear", "logistic", "quadratic"))))
    {
      stop(paste("Wrong relation type specified: pick among '", 
                 paste(c("gaussian", "linear", "logistic", "quadratic"), collapse = " "), "'",
                 collapse = " "))
    }
    valid.cells <- setValues(raster.stack[[1]], 1)
    var.order <- sample(var.names, length(var.names), replace = F)
    for (i in 1:length(var.order))
    {
      
      cur.var <- var.order[i]
      cur.rast <- raster.stack[[cur.var]]
      if(realistic.sp) cur.rast <- cur.rast * valid.cells # Cur.rast is here restricted to current suitable conds
      
      type <- sample(relations, 1)

      if (type == "gaussian")
      {
        parameters[[cur.var]] <- list(fun = 'dnorm',
                                      args = list(mean = sample(seq(cur.rast@data@min,
                                                                    cur.rast@data@max, 
                                                                    length = 100000), 1),
                                                  sd = sample(seq(0, 
                                                                  (raster.stack[[cur.var]]@data@max - raster.stack[[cur.var]]@data@min), 
                                                                  length = 100000), 1))
        )
      } else if (type == "linear")
      { # At the moment this is not really useful because the rescale will transforme the results in either 0:1 or 1:0, regardless of the slope
        # To be improved later
        parameters[[cur.var]] <- list(fun = 'linearFun',
                                      args = list(a = sample(seq(-1, 1, length = 100), 1),
                                                  b = sample(seq(raster.stack[[cur.var]]@data@min, 
                                                                 raster.stack[[cur.var]]@data@max, 
                                                                 length = 100000), 1))
        )
      } else if (type == "logistic")
      {
        beta.t <- sample(seq(raster.stack[[cur.var]]@data@min,
                             raster.stack[[cur.var]]@data@max,
                             length = 1000000), 1)
        alpha.t <-  sample(c(seq((raster.stack[[cur.var]]@data@max - raster.stack[[cur.var]]@data@min)/1000,
                                 (raster.stack[[cur.var]]@data@max - raster.stack[[cur.var]]@data@min)/100, length = 10),
                             seq((raster.stack[[cur.var]]@data@max - raster.stack[[cur.var]]@data@min)/100,
                                 (raster.stack[[cur.var]]@data@max - raster.stack[[cur.var]]@data@min)/10, length = 100),
                             seq((raster.stack[[cur.var]]@data@max - raster.stack[[cur.var]]@data@min)/10,
                                 (raster.stack[[cur.var]]@data@max - raster.stack[[cur.var]]@data@min)*10, length = 10)), size = 1)
        if(realistic.sp == TRUE)
        {
          if(beta.t > cur.rast@data@max)
          {
            alpha.t <- alpha.t
          } else if (beta.t < cur.rast@data@min)
          {
            alpha.t <- -alpha.t
          } else
          {
            alpha.t <- sample(c(alpha.t, -alpha.t), 1)
          }
        }
        
        parameters[[cur.var]] <- list(fun = 'logisticFun',
                                      args = list(alpha = alpha.t,
                                                  beta = beta.t)
        )
      } else if (type == "quadratic")
      {
        max.point <- sample(seq(cur.rast@data@min,
                                cur.rast@data@max, 
                                length = 1000), 1)
        a <- sample(seq(-.01, -20, length = 10000), 1)
        b <- - max.point * 2 * a
        parameters[[cur.var]] <- list(fun = 'quadraticFun',
                                      args = list(a = a,
                                                  b = b,
                                                  c = 0)
        )
        
      }
      
      # Restricting values to suitable conditions
      tmp.rast <- calc(raster.stack[[cur.var]], fun = function(x)
      {
        do.call(match.fun(parameters[[cur.var]]$fun), args = c(list(x), parameters[[cur.var]]$args))
      }
      )
      tmp.rast <- (tmp.rast - tmp.rast@data@min) / (tmp.rast@data@max - tmp.rast@data@min)
      valid.cells <- valid.cells * (tmp.rast > 0.05)
    }
    message(" - Calculating species suitability\n")
    results <- generateSpFromFun(raster.stack, parameters, species.type = species.type, plot = FALSE)
  }
  
  



  if(convert.to.PA == TRUE)
  {
    message(" - Converting into Presence - Absence\n")
    results <- convertToPA(results, 
                           PA.method = PA.method,
                           alpha = alpha,
                           beta = beta,
                           species.prevalence = species.prevalence,
                           plot = FALSE)
    
    if(plot) plot(stack(results$suitab.raster, results$pa.raster), main = c("Suitability", "Presence-absence"))
  } else
  {
    if(plot) plot(results$suitab.raster, main = "Suitability")
  }
  
  return(results)
}
