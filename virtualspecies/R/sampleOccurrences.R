#' Sample occurrences in a virtual species distribution
#' 
#' @description
#' This function samples presences within a species distribution, either
#' randomly or with a sampling bias. The sampling bias can be defined manually
#' or with a set of pre-defined biases.
#' 
#' @param x a \code{rasterLayer} object or the output list from 
#' \code{generateSpFromFun}, \code{generateSpFromPCA}, \code{generateRandomSp}, \code{convertToPA}
#' or  \code{limitDistribution}
#' The raster must contain values of 0 or 1 (or NA).
#' @param n an integer. The number of occurrence points to sample.
#' @param type \code{"presence only"} or \code{"presence-absence"}. The type of 
#' occurrence points to sample.
#' @param sampling.area a character string, a \code{polygon} or an \code{extent} object.
#' The area in which the sampling will take place. See details.
#' @param detection.probability a numeric value between 0 and 1, corresponding to the
#' probability of detection of the species. See details.
#' @param correct.by.suitability \code{TRUE} or \code{FALSE}. If \code{TRUE}, then
#' the probability of detection will be weighted by the suitability, such that 
#' cells with lower suitabilities will further decrease the chance that the species
#' is detected when sampled.
#' @param error.probability \code{TRUE} or \code{FALSE}. Only useful if 
#' \code{type = "presence-absence"}. Probability to attribute an erroneous presence
#' in cells where the species is absent.
#' @param bias  \code{"no.bias"},  \code{"country"},  \code{"region"},  
#' \code{"extent"},  \code{"polygon"} or \code{"manual"}. The method used to 
#' generate a sampling bias: see details.
#' @param bias.strength a positive numeric value. The strength of the bias to be applied
#' in \code{area} (as a multiplier). Above 1, \code{area} will be oversampled. Below 1, \code{area}
#' will be undersampled.
#' @param bias.area \code{NULL}, a character string, a \code{polygon} or an \code{extent} object.
#' The area in which the sampling will be biased: see details. If \code{NULL}
#' and \code{bias = "extent"}, then you will be asked to draw an
#' extent on the map.
#' @param weights \code{NULL} or a raster layer. Only used if \code{bias = "manual"}.
#' The raster of bias weights to be applied to the sampling of occurrences.
#' Higher weights mean a higher probability of sampling.
#' @param sample.prevalence \code{NULL} or a numeric value between 0 and 1. Only useful if 
#' \code{type = "presence-absence"}. Defines the sample prevalence, i.e. the proportion of presences
#' sampled. Note that the probabilities of detection and error are applied AFTER this parameter,
#' so the final sample prevalence may not different if you apply probabilities of detection and/or error 
#' @param plot \code{TRUE} or \code{FALSE}. If \code{TRUE}, the sampled occurrence
#' points will be plotted.
#' @details
#' \bold{How the function works:}
#' 
#' The function randomly selects \code{n} cells in which samples occur. If a \code{bias}
#' is chosen, then the selection of these cells will be biased according to the type and
#' strength of bias chosen. If the sampling is of \code{type "presence only"}, then only
#' cells where the species is present will be chosen. If the sampling is of 
#' \code{type "presence-absence"}, then all non-NA cells can be chosen.
#' 
#' The function then samples the species inside the chosen cells. In cells 
#' where the species is present the species will always be sampled unless 
#' the parameter \code{detection.probability} is lower than 1. In that case the
#' species will be sampled with the associated probability of detection.
#' 
#' In cells where the species is absent (in case of a \code{"presence-absence"}
#' sampling), the function will always assign absence unless \code{error.probability} 
#' is greater than 1. In that case, the species can be found present with the 
#' associated probability of error. Note that this step happens AFTER the detection
#' step. Hence, in cells where the species is present but not detected, it can
#' still be sampled due to a sampling error.
#' 
#' \bold{How to restrict the sampling area:}
#' 
#' Use the argument \code{sampling.area}:
#' \itemize{
#' \item{Provide the name (s) (or a combination of names) of country(ies), region(s) or continent(s).
#' Examples:
#' \itemize{
#' \item{\code{sampling.area = "Africa"}}
#' \item{\code{sampling.area = c("Africa", "North America", "France")}}
#' }}
#' \item{Provide a polygon (\code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} 
#' of package \code{sp})}
#' \item{Provide an \code{extent} object}
#' }
#' 
#' \bold{How the sampling bias works:}
#' 
#' The argument \code{bias.strength} indicates the strength of the bias.
#' For example, a value of 50 will result in 50 times more samples within the
#'  \code{bias.area} than outside.
#' Conversely, a value of 0.5 will result in half less samples within the
#' \code{bias.area} than outside.
#' 
#' \bold{How to choose where the sampling is biased:}
#' 
#' You can choose to bias the sampling in:
#' \enumerate{
#' \item{a particular country, region or continent (assuming your raster has
#' the WGS84 projection): 
#' 
#' Set the argument
#' \code{bias} to \code{"country"}, \code{"region"} or
#' \code{"continent"}, and provide the name(s) of the associated countries,
#' regions or continents to \code{bias.area} (see examples). 
#' 
#' List of possible \code{bias.area} names:
#' \itemize{
#' \item{Countries: type \code{levels(getMap()@@data$SOVEREIGNT)} in the console}
#' \item{Regions: "Africa", "Antarctica", "Asia", "Australia", "Europe", 
#' "North America", "South America"}
#' \item{Continents: "Africa", "Antarctica", "Australia", "Eurasia", 
#' "North America", "South America"}}
#' }
#' \item{a polygon:
#' 
#' Set \code{bias} to \code{"polygon"}, and provide your
#' polygon to \code{area}.
#' }
#' \item{an extent object:
#' 
#' Set \code{bias} to \code{"extent"}, and either provide your
#' extent object to \code{bias.area}, or leave it \code{NULL} to draw an extent on
#' the map.}
#' } 
#' 
#' Otherwise you can define manually your sampling bias, \emph{e.g.} to create
#' biases along roads. In that case you have to provide to \code{weights} a raster layer in which
#' each cell contains the probability to be sampled.
#' @return a \code{list} with 3 (unbiased sampling) to 4 (biased sampling) elements:
#' \itemize{
#' \item{\code{sample.points}: the data.frame containing the coordinates of 
#' samples, the real presence-absences (or presence-only) and the sampled presence-
#' absences}
#' \item{\code{detection.probability}: the chosen probability of detection of
#' the virtual species}
#' \item{\code{error.probability}: the chosen probability to assign presence
#' in cells where the species is absent}
#' \item{\code{bias}: if a bias was chosen, then the type of bias and the
#' associated \code{area} will be included.}
#' }
#' @export
#' @import raster
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, C. Bellard & F. Courchamp
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
#' sp <- generateRandomSp(env, niche.breadth = "wide")
#' 
#' # Sampling of 25 presences
#' sampleOccurrences(sp, n = 25)
#' 
#' # Sampling of 30 presences and absebces
#' sampleOccurrences(sp, n = 30, type = "presence-absence")
#' 
#' # Reducing of the probability of detection
#' sampleOccurrences(sp, n = 30, type = "presence-absence", 
#'                   detection.probability = 0.5)
#'                   
#' # Further reducing in relation to environmental suitability
#' sampleOccurrences(sp, n = 30, type = "presence-absence", 
#'                   detection.probability = 0.5,
#'                   correct.by.suitability = TRUE)
#'                   
#' # Creating sampling errors (far too much)
#' sampleOccurrences(sp, n = 30, type = "presence-absence", 
#'                   error.probability = 0.5)
#'                   
#' # Introducing a sampling bias (oversampling)
#' biased.area <- extent(0.5, 0.7, 0.6, 0.8)
#' sampleOccurrences(sp, n = 50, type = "presence-absence", 
#'                   bias = "extent",
#'                   bias.area = biased.area)
#' # Showing the area in which the sampling is biased
#' plot(biased.area, add = TRUE)     
#' 
#' # Introducing a sampling bias (no sampling at all in the chosen area)
#' biased.area <- extent(0.5, 0.7, 0.6, 0.8)
#' sampleOccurrences(sp, n = 50, type = "presence-absence", 
#'                   bias = "extent",
#'                   bias.strength = 0,
#'                   bias.area = biased.area)
#' # Showing the area in which the sampling is biased
#' plot(biased.area, add = TRUE)             



sampleOccurrences <- function(x, n,
                              type = "presence only",
                              sampling.area = NULL,
                              detection.probability = 1,
                              correct.by.suitability = FALSE,
                              error.probability = 0,
                              bias = "no.bias", 
                              bias.strength = 50,
                              bias.area = NULL,
                              weights = NULL,
                              sample.prevalence = NULL,
                              plot = TRUE)
{
  results <- list()
  
  if("virtualspecies" %in% class(x))
  {
    if("RasterLayer" %in% class(x$occupied.area))
    {
      sp.raster <- x$occupied.area
    } else if("RasterLayer" %in% class(x$pa.raster))
    {
      sp.raster <- x$pa.raster
    } else stop("x must be:\n- a raster layer object\nor\n- the output list from functions 
                generateRandomSp(), convertToPA() or limitDistribution()")
  } else if ("RasterLayer" %in% class(x))
  {
    sp.raster <- x
  } else stop("x must be:\n- a raster layer object\nor\n- the output list from functions 
              generateRandomSp(), convertToPA() or limitDistribution()")
  
  if(sp.raster@data@max > 1 | sp.raster@data@min < 0)
  {
    stop("There are values above 1 or below 0 in your presence/absence raster. 
         Please make sure that the provided raster is a correct P/A raster and not a suitability raster.")
  }
  
  original.raster <- sp.raster
  
  if(!is.null(sample.prevalence))
  {
    if(sample.prevalence < 0 | sample.prevalence > 1)
    {
      stop("Sample prevalence must be a numeric between 0 and 1")
    }
  }
  
  if(!is.null(sampling.area))
  {
    if(is.character(sampling.area))
    {
      worldmap <- rworldmap::getMap()
      if (any(!(sampling.area %in% c(levels(worldmap@data$SOVEREIGNT),
                                     levels(worldmap@data$REGION),
                                     levels(worldmap@data$continent)))))
      {
        stop("The choosen sampling.area is incorrectly spelled.\n Type 'levels(getMap()@data$SOVEREIGNT)', 'levels(worldmap@data$REGION)' and levels(worldmap@data$continent) to obtain valid names.")
      }
      sampling.area <- worldmap[which(worldmap@data$SOVEREIGNT %in% sampling.area | 
                                               worldmap@data$REGION %in% sampling.area |
                                               worldmap@data$continent %in% sampling.area), ]
    } else if(!(class(sampling.area) %in% c("SpatialPolygons", "SpatialPolygonsDataFrame", "Extent")))
    {
      stop("Please provide to sampling.area either \n
           - the names of countries, region and/or continents in which to sample\n
           - a SpatialPolygons or SpatialPolygonsDataFrame\n
           - an extent\n
           in which the sampling will take place")
    }
    
    sample.area.raster1 <- rasterize(sampling.area,
                                            sp.raster, 
                                            field = 1,
                                            background = NA,
                                     silent = TRUE)
    sp.raster <- sp.raster * sample.area.raster1
  }
  
  
  if(correct.by.suitability)
  {
    if(!("virtualspecies" %in% class(x)) | !("suitab.raster" %in% names(x)))
    {
      stop("If you choose to weight the probability of detection by the suitability of the species (i.e., correct.by.suitability = TRUE),
           then you need to provide an appropriate virtual species containing a suitability raster to x.")
    }
  }
  
  if(!is.numeric(detection.probability) | detection.probability > 1 | detection.probability < 0)
  {
    stop("detection.probability must be a numeric value between 0 and 1")
  }
  
  if(!is.numeric(error.probability) | error.probability > 1 | error.probability < 0)
  {
    stop("error.probability must be a numeric value between 0 and 1")
  }
  
  if(length(bias) > 1)
  {
    stop('Only one bias can be applied at a time')
  }
  
  if (!(bias %in% c("no.bias", "country", "region", "continent", "extent", "polygon", "manual")))
  {
    stop('Argument bias must be one of : "no.bias", "country", "region", "continent", "extent", "polygon", "manual"')
  }
  
  if(!is.numeric(bias.strength) & bias != "no.bias")
  {
    stop("Please provide a numeric value for bias.strength")
  }
  
  if (bias %in% c("country", "region", "continent"))
  {
    if(!("rworldmap" %in% rownames(installed.packages())))
    {
      stop('You need to install the package "rworldmap" in order to use bias = "region" or bias = "country"')
    }
    worldmap <- rworldmap::getMap()
    
    if(bias == "country")
    {
      if (any(!(bias.area %in% levels(worldmap@data$SOVEREIGNT))))
      {
        stop("country name(s) must be correctly spelled. Type 'levels(getMap()@data$SOVEREIGNT)' to obtain valid names.")
      }    
      results$bias <- list(bias = bias,
                           bias.strength = bias.strength,
                           bias.area = bias.area)
    } else if (bias == "region")
    {
      if (any(!(bias.area %in% levels(worldmap@data$REGION))))
      {
        stop(paste("region name(s) must be correctly spelled, according to one of the following : ", 
                   paste(levels(worldmap@data$REGION), collapse = ", "), sep = "\n"))
      } 
      results$bias <- list(bias = bias,
                           bias.strength = bias.strength,
                           bias.area = bias.area)
    } else if (bias == "continent")
    {
      if (any(!(bias.area %in% levels(worldmap@data$continent))))
      {
        stop(paste("region name(s) must be correctly spelled, according to one of the following : ", 
                   paste(levels(worldmap@data$continent), collapse = ", "), sep = "\n"))
      } 
      results$bias <- list(bias = bias,
                           bias.strength = bias.strength,
                           bias.area = bias.area)
    }
  }
  if (bias == "polygon") # Projections are not checked here. Perhaps we should add projection check between raster & polygon in the future?
                                # This is especially important given that randomPoints weights samplings by the cell area (because cells closer to
                                # the equator are larger)
  {
    if(!(class(bias.area) %in% c("SpatialPolygons", "SpatialPolygonsDataFrame")))
    {
      stop("If you choose bias = 'polygon', please provide a polygon of class SpatialPolygons or SpatialPolygonsDataFrame to argument bias.area")
    }
    warning("Polygon projection is not checked. Please make sure you have the same projections between your polygon and your presence-absence raster")
    results$bias <- list(bias = bias,
                         bias.strength = bias.strength,
                         bias.area = bias.area)
  }

  if (bias == "extent")
  {
    results$bias <- list(bias = bias,
                         bias.strength = bias.strength,
                         bias.area = bias.area)
  }
  
  if(type == "presence-absence")
  {
    sample.raster <- sp.raster
    sample.raster[!is.na(sample.raster)] <- 1
  } else if (type == "presence only")
  {
    sample.raster <- sp.raster
  } else stop("type must either be 'presence only' or 'presence-absence'")
  
  if (bias == "manual")
  {
    if(!("RasterLayer" %in% class(weights)))
    {
      stop("You must provide a raster layer of weights (to argument weights) if you choose bias == 'manual'")
    }
    bias.raster <- weights * sample.raster
    results$bias <- list(bias = bias,
                         weights = weights)
  } else
  {
    bias.raster <- sample.raster    
  }
  
  if(bias == "country")
  {
    bias.raster1 <- rasterize(worldmap[which(worldmap@data$SOVEREIGNT %in% bias.area), ],
                              bias.raster, 
                              field = bias.strength,
                              background = 1,
                              silent = TRUE)
    bias.raster <- bias.raster * bias.raster1
  } else if(bias == "region")
  {
    bias.raster1 <- rasterize(worldmap[which(worldmap@data$REGION %in% bias.area), ],
                              bias.raster, 
                              field = bias.strength,
                              background = 1,
                              silent = TRUE)
    bias.raster <- bias.raster * bias.raster1
  } else if(bias == "continent")
  {
    bias.raster1 <- rasterize(worldmap[which(levels(worldmap@data$continent) %in% bias.area), ],
                              bias.raster, 
                              field = bias.strength,
                              background = 1,
                              silent = TRUE)
    bias.raster <- bias.raster * bias.raster1
  } else if(bias == "extent")
  {
    if(!(class(bias.area) == "Extent"))
    {
      message("No object of class extent provided: click twice on the map to draw the extent in which presence points will be sampled")
      plot(sp.raster)
      bias.area <- drawExtent(show = TRUE)      
    }
    bias.raster <- bias.raster * rasterize(bias.area, sp.raster, field = bias.strength, background = 1)
    results$bias <- list(bias = bias,
                         bias.area = bias.area)
  } else if(bias == "polygon")
  {
    bias.raster1 <- rasterize(bias.area,
                              bias.raster, 
                              field = bias.strength,
                              background = 1,
                              silent = TRUE)
    bias.raster <- bias.raster * bias.raster1
  }

  
  if(bias != "no.bias")
  {
    if(type == "presence only")
    {
      sample.points <- dismo::randomPoints(sample.raster * bias.raster, n = n, prob = TRUE, tryf = 2)
    } else
    {
      if(is.null(sample.prevalence))
      {
        sample.points <- dismo::randomPoints(sample.raster * bias.raster, n = n, prob = TRUE)
      } else
      {
        tmp1 <- sample.points
        tmp1[sp.raster != 1] <- NA
        sample.points <- dismo::randomPoints(tmp1 * bias.raster, n = sample.prevalence * n, prob = TRUE)
        tmp1 <- sample.raster
        tmp1[sp.raster != 0] <- NA
        sample.points <- rbind(sample.points,
                               dismo::randomPoints(tmp1 * bias.raster, n = (1 - sample.prevalence) * n, prob = TRUE))
        rm(tmp1)
      }
    }
  } else
  {
    if(type == "presence only")
    {
      sample.points <- dismo::randomPoints(sample.raster, n = n, prob = TRUE, tryf = 2)
    } else
    {
      if(is.null(sample.prevalence))
      {
        sample.points <- dismo::randomPoints(sample.raster, n = n, prob = TRUE)
      } else
      {
        tmp1 <- sample.raster
        tmp1[sp.raster != 1] <- NA
        sample.points <- dismo::randomPoints(tmp1, n = sample.prevalence * n, prob = TRUE)
        tmp1 <- sample.raster
        tmp1[sp.raster != 0] <- NA
        sample.points <- rbind(sample.points,
                               dismo::randomPoints(tmp1, n = (1 - sample.prevalence) * n, prob = TRUE))
        rm(tmp1)
      }
    }
  }
  
  if(type == "presence only")
  {
    if(error.probability != 0)
    {
      warning("The error probability has no impact when sampling 'presence only' 
               points, because these samplings occur only within the boundaries of the species range")
    }
    sample.points <- data.frame(sample.points,
                                Real = 1,
                                Observed = sample(c(NA, 1),
                                                  size = nrow(sample.points),
                                                  prob = c(1 - detection.probability,
                                                           detection.probability),
                                                  replace = TRUE))
  } else if(type == "presence-absence")
  {
    sample.points <- data.frame(sample.points,
                                Real = extract(sp.raster, sample.points))
    
    if(correct.by.suitability)
    {
      suitabs <- extract(x$suitab.raster, sample.points[, c("x", "y")])
    } else { suitabs <- rep(1, nrow(sample.points)) }
    
    sample.points$Observed <- NA
    if(correct.by.suitability)
    {
      sample.points$Observed[which(sample.points$Real == 1)] <-
        sapply(detection.probability * suitabs[which(sample.points$Real == 1)],
               function(y)
               {
                 sample(c(0, 1),
                        size = 1,
                        prob = c(1 - y, y))
               })
    } else
    {
      sample.points$Observed[which(sample.points$Real == 1)] <-
        sample(c(0, 1), size = length(which(sample.points$Real == 1)),
               prob = c(1 - detection.probability, detection.probability),
               replace = TRUE)
    }
    sample.points$Observed[which(sample.points$Real == 0 | sample.points$Observed == 0)] <-
      sample(c(0, 1), size = length(which(sample.points$Real == 0 | sample.points$Observed == 0)),
             prob = c(1 - error.probability, error.probability),
             replace = TRUE)
    
  }
  
  if(plot)
  {
    plot(original.raster)
    if(type == "presence only")
    {
      points(sample.points[, c("x", "y")], pch = 16, cex = .5)
    } else
    {
      points(sample.points[sample.points$Observed == 1, c("x", "y")], pch = 16, cex = .8)
      points(sample.points[sample.points$Observed == 0, c("x", "y")], pch = 1, cex = .8)
    }
  }


  results$sample.points <- sample.points
  results$detection.probability <- detection.probability
  results$error.probability <- error.probability
  if(type == "presence-absence")
  {

    true.prev <- length(sample.points$Real[which(sample.points$Real == 1)]) / nrow(sample.points)
    obs.prev <- length(sample.points$Real[which(sample.points$Observed == 1)]) / nrow(sample.points)
    
    results$sample.prevalence <- c(true.sample.prevalence = true.prev,
                                   observed.sample.prevalence = obs.prev)
  }

  
  class(results) <- append(class(results), "VSSampledPoints")
  return(results)
}