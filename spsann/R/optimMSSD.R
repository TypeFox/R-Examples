#' Optimization of sample configurations for spatial interpolation (I)
#'
#' Optimize a sample configuration for spatial interpolation. The criterion used is the mean squared shortest 
#' distance (\bold{MSSD}) between sample points and prediction points.
#'
#' @inheritParams spJitter
#' @template spSANN_doc
#' @template spJitter_doc
#' 
#' @details 
#' \subsection{Spatial coverage sampling}{
#' Spatial coverage sampling is based on the knowledge that the kriging variance depends upon the distance 
#' between sample points. As such, the better the spread of the sample points in the spatial domain, the 
#' smaller the kriging variance. This is similar to using a regular grid of sample points. However, 
#' a regular grid usually is suboptimal for irregularly shaped areas.
#' }
#'
#' @return
#' \code{optimMSSD} returns an object of class \code{OptimizedSampleConfiguration}: the optimized sample
#' configuration with details about the optimization.
#'
#' \code{objMSSD} returns a numeric value: the energy state of the sample configuration -- the objective
#' function value.
#' 
#' @note
#' This function was derived with modifications from the method known as \emph{spatial coverage sampling} 
#' originally proposed by Brus, de Gruijter and van Groenigen (2006), and implemented in the R-package 
#' \pkg{\link[spcosa]{spcosa}} by Dennis Walvoort, Dick Brus and Jaap de Gruijter.
#' 
#' @references
#' Brus, D. J.; de Gruijter, J. J.; van Groenigen, J.-W. Designing spatial coverage samples using the k-means
#' clustering algorithm. In: P. Lagacherie,A. M.; Voltz, M. (Eds.) \emph{Digital soil mapping -- an 
#' introductory perspective}. Elsevier, v. 31, p. 183-192, 2006.
#'
#' de Gruijter, J. J.; Brus, D.; Bierkens, M.; Knotters, M. \emph{Sampling for natural resource monitoring}.
#' Berlin: Springer, p. 332, 2006.
#'
#' Walvoort, D. J. J.; Brus, D. J.; de Gruijter, J. J. An R package for spatial coverage sampling and random 
#' sampling from compact geographical strata by k-means. \emph{Computers and Geosciences}. v. 36, p. 
#' 1261-1267, 2010.
#'
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso
#' \code{\link[raster]{distanceFromPoints}}, \code{\link[spcosa]{stratify}}.
#' @aliases optimMSSD objMSSD MSSD
#' @concept spatial interpolation
#' @export
#' @examples
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' schedule <- scheduleSPSANN(chains = 1, initial.temperature = 5000000,
#'                            x.max = 1540, y.max = 2060, x.min = 0, 
#'                            y.min = 0, cellsize = 40)
#' set.seed(2001)
#' res <- optimMSSD(points = 10, candi = candi, schedule = schedule)
#' objSPSANN(res) - objMSSD(candi = candi, points = res)
#' 
# FUNCTION - MAIN ##############################################################
optimMSSD <-
  function (points, candi,
            # SPSANN
            schedule = scheduleSPSANN(), plotit = FALSE, track = FALSE, 
            boundary, progress = "txt", verbose = FALSE) {
    
    # Objective function name
    objective <- "MSSD"
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Calculate the initial energy state. The distance matrix is calculated
    # using the SpatialTools::dist2(). The function .objMSSD() does the
    # squaring internaly.
    dm <- SpatialTools::dist2(candi[, 2:3], conf0[, 2:3])
    energy0 <- data.frame(obj = .objMSSD(x = dm))
    
    # Other settings for the simulated annealing algorithm
    old_dm <- dm
    best_dm <- dm
    old_energy <- energy0
    best_energy <- data.frame(obj = Inf)
    actual_temp <- schedule$initial.temperature
    k <- 0 # count the number of jitters
    
    # Set progress bar
    eval(.set_progress())
    
    # Initiate the annealing schedule
    for (i in 1:schedule$chains) {
      n_accept <- 0
      
      for (j in 1:schedule$chain.length) { # Initiate one chain
        
        for (wp in 1:n_pts) { # Initiate loop through points
          k <- k + 1
          
          # Plotting and jittering
          eval(.plot_and_jitter())
          
          # Update the matrix of distances and calculate the new energy state
          x2 <- matrix(new_conf[wp, 2:3], nrow = 1)
          new_dm <- .updateMSSDCpp(x1 = candi[, 2:3], x2 = x2, dm = old_dm, 
                                   idx = wp)
          new_energy <- data.frame(obj = .objMSSD(new_dm))
          
          # Update the matrix of distances and calculate the new energy state
          # x2 <- SpatialTools::dist2(coords = candi[, 2:3], coords2 = x2)
          # new_dm <- old_dm
          # new_dm[, wp] <- x2
          # new_energy <- mean(apply(new_dm, 1, min) ^ 2)
          
          # Evaluate the new system configuration
          accept <- .acceptSPSANN(old_energy[[1]], new_energy[[1]], actual_temp)
          if (accept) {
            old_conf <- new_conf
            old_energy <- new_energy
            old_dm <- new_dm
            n_accept <- n_accept + 1
          } else {
            new_energy <- old_energy
            new_conf <- old_conf
            # new_dm <- old_dm
          }
          if (track) energies[k, ] <- new_energy
          
          # Record best energy state
          if (new_energy[[1]] < best_energy[[1]] / 1.0000001) {
            best_k <- k
            best_conf <- new_conf
            best_energy <- new_energy
            best_old_energy <- old_energy
            old_conf <- old_conf
            best_dm <- new_dm
            best_old_dm <- old_dm
          }
          
          # Update progress bar
          eval(.update_progress())
          
        } # End loop through points
        
      } # End the chain
      
      # Check the proportion of accepted jitters in the first chain
      eval(.check_first_chain())
      
      # Count the number of chains without any change in the objective function.
      # Restart with the previously best configuration if it exists.
      if (n_accept == 0) {
        no_change <- no_change + 1
        if (no_change > schedule$stopping) {
          # if (new_energy[[1]] > best_energy[[1]] * 1.000001) {
            # old_conf <- old_conf
            # new_conf <- best_conf
            # old_energy <- best_old_energy
            # new_energy <- best_energy
            # new_dm <- best_dm
            # old_dm <- best_old_dm
            # no_change <- 0
            # new_sm <- best_sm
            # old_sm <- best_old_sm
            # cat("\nrestarting with previously best configuration\n")
          # } else { 
            break 
          # }
        }
        if (verbose) {
          cat("\n", no_change, "chain(s) with no improvement... stops at",
              schedule$stopping, "\n")
        }
      } else {
        no_change <-  0
      }
      
      # Update control parameters
      actual_temp <- actual_temp * schedule$temperature.decrease
      x.max <- x_max0 - (i / schedule$chains) * (x_max0 - x.min) + cellsize[1]
      y.max <- y_max0 - (i / schedule$chains) * (y_max0 - y.min) + cellsize[2]
      
    } # End the annealing schedule
    
    # Prepare output
    eval(.prepare_output())
    return (res)
  }
# FUNCTION - CALCULATE THE CRITERION VALUE #####################################
#' @rdname optimMSSD
#' @export
objMSSD <-
  function (points, candi) {
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Compute the matrix of distances
    dm <- SpatialTools::dist2(coords1 = candi[, 2:3], coords2 = points[, 2:3])
    
    # Calculate the energy state
    res <- .objMSSD(dm)
    
    # Output
    return (res)
  }
# INTERNAL FUNCTION - NEAREST POINT ############################################
# Function to identify the nearest point.
.nearestPoint <-
  function (dist.mat, which.pts) {
    if (dim(dist.mat)[1] != dim(dist.mat)[2] || !is.matrix(dist.mat)) {
      stop ("'dist.mat' should be a n x n matrix")
    }
    if (!is.vector(which.pts) || length(which.pts) >= dim(dist.mat)[2]) {
      stop (paste("'which.pts' should be a vector of length smaller than ",
                  dim(dist.mat)[2], sep = ""))
    }
    sub_dm <- dist.mat[, which.pts]
    min_dist_id <- apply(sub_dm, MARGIN = 1, FUN = which.min)
    return (min_dist_id)
  }
# 
.subsetCandi <-
  function (candi, cellsize) {
    x <- seq(min(candi[, 1]), max(candi[, 1]), by = cellsize * 2)
    y <- seq(min(candi[, 2]), max(candi[, 2]), by = cellsize * 2)
    xy <- expand.grid(x, y)
    colnames(xy) <- colnames(candi)
    dup <- rbind(xy, candi)
    dup <- which(duplicated(dup)) - nrow(xy)
    return (dup)
  }
# End!
