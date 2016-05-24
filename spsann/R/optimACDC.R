#' Optimization of sample configurations for spatial trend identification and estimation (III)
#'
#' Optimize a sample configuration for spatial trend identification and estimation. An utility function 
#' \emph{U} is defined so that the sample reproduces the bivariate association/correlation between the 
#' covariates, as well as their marginal distribution (\bold{ACDC}). The utility function is obtained 
#' aggregating two objective functions: \bold{CORR} and \bold{DIST}.
#'
#' @inheritParams spJitter
#' @template spSANN_doc
#' @template ACDC_doc
#' @template MOOP_doc
#' @template spJitter_doc
#' 
#' @details 
#' Visit the help pages of \code{\link[spsann]{optimCORR}} and \code{\link[spsann]{optimDIST}} to see the 
#' details of the objective functions that compose \bold{ACDC}.
#' 
#' @return
#' \code{optimACDC} returns an object of class \code{OptimizedSampleConfiguration}: the optimized sample
#' configuration with details about the optimization.
#' 
#' \code{objACDC} returns a numeric value: the energy state of the sample configuration -- the objective
#' function value.
#' 
#' @note
#' This function was derived with modifications from the method known as the \emph{conditioned Latin Hypercube 
#' sampling} originally proposed by Minasny and McBratney (2006), and implemented in the R-package 
#' \pkg{\link[clhs]{clhs}} by Pierre Roudier.
#' 
#' @references
#' Minasny, B.; McBratney, A. B. A conditioned Latin hypercube method for sampling in the presence of 
#' ancillary information. \emph{Computers & Geosciences}, v. 32, p. 1378-1388, 2006.
#'
#' Minasny, B.; McBratney, A. B. Conditioned Latin Hypercube Sampling for calibrating soil sensor data to soil
#' properties. Chapter 9. Viscarra Rossel, R. A.; McBratney, A. B.; Minasny, B. (Eds.) \emph{Proximal Soil
#' Sensing}. Amsterdam: Springer, p. 111-119, 2010.
#'
#' Roudier, P.; Beaudette, D.; Hewitt, A. A conditioned Latin hypercube sampling algorithm incorporating
#' operational constraints. \emph{5th Global Workshop on Digital Soil Mapping}. Sydney, p. 227-231, 2012.
#' 
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[clhs]{clhs}}, \code{\link[pedometrics]{cramer}}
#' @aliases optimACDC objACDC ACDC
#' @export
#' @examples
#' data(meuse.grid, package = "sp")
#' candi <- meuse.grid[1:1000, 1:2]
#' nadir <- list(sim = 10, seeds = 1:10)
#' utopia <- list(user = list(DIST = 0, CORR = 0))
#' covars <- meuse.grid[1:1000, 5]
#' schedule <- scheduleSPSANN(
#'   chains = 1, initial.temperature = 5, x.max = 1540, y.max = 2060, 
#'   x.min = 0, y.min = 0, cellsize = 40)
#' set.seed(2001)
#' res <- optimACDC(
#'   points = 10, candi = candi, covars = covars, nadir = nadir,
#'   use.coords = TRUE, utopia = utopia, schedule = schedule)
#' objSPSANN(res) - objACDC(
#'   points = res, candi = candi, covars = covars, 
#'   use.coords = TRUE, nadir = nadir, utopia = utopia)
# MAIN FUNCTION ################################################################
optimACDC <-
  function (points, candi, 
            # DIST and CORR
            covars, strata.type = "area", use.coords = FALSE, 
            # SPSANN
            schedule = scheduleSPSANN(), plotit = FALSE, track = FALSE,
            boundary, progress = "txt", verbose = FALSE,
            # MOOP
            weights = list(CORR = 0.5, DIST = 0.5),
            nadir = list(sim = NULL, seeds = NULL, user = NULL, abs = NULL),
            utopia = list(user = NULL, abs = NULL)) {
    
    # Objective function name
    objective <- "ACDC"
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Check other arguments
    check <- .optimACDCcheck(
      candi = candi, covars = covars, use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Prepare 'covars' and create the starting sample matrix 'sm'
    eval(.prepare_acdc_covars())
    
    # Base data and initial energy state
    pcm <- .corCORR(obj = covars, covars.type = covars.type)
    scm <- .corCORR(obj = sm, covars.type = covars.type)
    pop_prop <- .strataACDC(
      n.pts = n_pts, strata.type = strata.type, covars = covars, covars.type = covars.type)
    nadir <- .nadirACDC(
      n.pts = n_pts, n.cov = n_cov, n.candi = n_candi, pcm = pcm, nadir = nadir, covars.type = covars.type,
      covars = covars, pop.prop = pop_prop, candi = candi)
    utopia <- .utopiaACDC(utopia = utopia)
    energy0 <- .objACDC(
      sm = sm, n.cov = n_cov, pop.prop = pop_prop, pcm = pcm, scm = scm, nadir = nadir, weights = weights, 
      n.pts = n_pts, utopia = utopia, covars.type = covars.type)
    
    # Other settings for the simulated annealing algorithm
    old_sm <- sm
    new_sm <- sm
    best_sm <- sm
    old_scm <- scm
    new_scm <- scm
    best_scm <- scm
    old_energy <- energy0
    best_energy <- data.frame(obj = Inf, CORR = Inf, DIST = Inf)
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
          
          # Update sample and correlation matrices, and energy state
          new_sm[wp, ] <- covars[new_conf[wp, 1], ]
          new_scm <- .corCORR(obj = new_sm, covars.type = covars.type)
          new_energy <- .objACDC(
            sm = new_sm, pop.prop = pop_prop, scm = new_scm, nadir = nadir, weights = weights, pcm = pcm, 
            n.pts = n_pts, n.cov = n_cov, utopia = utopia, covars.type = covars.type)
          
          # Avoid the following error:
          # Error in if (new_energy[1] <= old_energy[1]) { : 
          #   missing value where TRUE/FALSE needed
          # Source: http://stackoverflow.com/a/7355280/3365410
          # ASR: The reason for the error is unknown to me.
          if (is.na(new_energy[1])) {
            new_energy <- old_energy
            new_conf <- old_conf
            new_sm <- old_sm
            new_scm <- old_scm
          }
          
          # Evaluate the new system configuration
          accept <- .acceptSPSANN(old_energy[[1]], new_energy[[1]], actual_temp)
          if (accept) {
            old_conf <- new_conf
            old_energy <- new_energy
            old_sm <- new_sm
            old_scm <- new_scm
            n_accept <- n_accept + 1
          } else {
            new_energy <- old_energy
            new_conf <- old_conf
            new_sm <- old_sm
            new_scm <- old_scm
          }
          if (track) energies[k, ] <- new_energy
          
          # Record best energy state
          if (new_energy[[1]] < best_energy[[1]] / 1.0000001) {
            best_k <- k
            best_conf <- new_conf
            best_energy <- new_energy
            best_old_energy <- old_energy
            old_conf <- old_conf
            best_sm <- new_sm
            best_old_sm <- old_sm
            best_scm <- new_scm
            best_old_scm <- old_scm
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
            # new_sm <- best_sm
            # new_scm <- best_scm
            # old_sm <- best_old_sm
            # old_scm <- best_old_scm
            # no_change <- 0
            # cat("\nrestarting with previously best configuration\n")
          # } else {
            break 
          # }
        }
        if (verbose) {
          cat("\n", no_change, "chain(s) with no improvement... stops at", schedule$stopping, "\n")
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
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
# candi: candidate locations
# covars: covariates
# use.coords: should the coordinates be used
# strata.type: type of stratification of numeric covariates
.optimACDCcheck <-
  function (candi, covars, use.coords, strata.type) {
    
    # covars
    if (is.vector(covars)) {
      if (use.coords == FALSE) {
        res <- "'covars' must have two or more columns"
        return (res)
      }
      if (nrow(candi) != length(covars)) {
        res <- "'candi' and 'covars' must have the same number of rows"
        return (res)
      }
    } else {
      if (nrow(candi) != nrow(covars)) {
        res <- "'candi' and 'covars' must have the same number of rows"
        return (res)
      }
    }
    
    # strata.type
    # aa <- match(strata.type, c("area", "range"))
    if (!strata.type %in% c("area", "range")) {
      res <- paste("'strata.type = ", strata.type, "' is not supported", sep = "")
      return (res)
    }
  }
# INTERNAL FUNCTION - BREAKS FOR NUMERIC COVARIATES ############################
# Now we define the breaks and the distribution, and return it as a list.
# Quantiles now honour the fact that the data are discontinuous.
# NOTE: there might be a problem when the number of unique values is small (3)
.strataACDC <-
  function (n.pts, covars, strata.type, covars.type) {
    
    n_cov <- ncol(covars)
    
    if (covars.type == "factor") {
      
      # Compute the proportion of population points per marginal factor level
      res <- lapply(covars, function(x) table(x) / nrow(covars))
      
    } else { # Numeric covariates
      
      # equal area strata
      if (strata.type == "area") {
        
        # Compute the break points (discrete sample quantiles)
        probs <- seq(0, 1, length.out = n.pts + 1)
        breaks <- lapply(covars, stats::quantile, probs, na.rm = TRUE, type = 3)
        
      } else { # equal range strata
        
        # Compute the break points
        breaks <- lapply(1:n_cov, function(i)
          seq(min(covars[, i]), max(covars[, i]), length.out = n.pts + 1))
        
        # Find and replace by the closest population value
        d <- lapply(1:n_cov, function(i)
          SpatialTools::dist2(matrix(breaks[[i]]), matrix(covars[, i])))
        d <- lapply(1:n_cov, function(i) apply(d[[i]], 1, which.min))
        breaks <- lapply(1:n_cov, function(i) breaks[[i]] <- covars[d[[i]], i])
      }
      
      # Keep only the unique break points
      breaks <- lapply(breaks, unique)
      
      # Compute the proportion of population points per marginal sampling strata
      count <- lapply(1:n_cov, function (i) 
        graphics::hist(covars[, i], breaks[[i]], plot = FALSE)$counts
      )
      prop <- lapply(1:n_cov, function (i) count[[i]] / sum(count[[i]]))
      
      # Output
      res <- list(breaks = breaks, prop = prop)  
    }
    
    return (res)
  }
# INTERNAL FUNCTION - COMPUTE THE NADIR VALUE ##################################
.nadirACDC <-
  function (n.pts, n.cov, n.candi, nadir, candi, covars, pcm, pop.prop, covars.type) {
    
    # Simulate the nadir point
    if (!is.null(nadir$sim) && !is.null(nadir$seeds)) { 
      m <- paste("simulating ", nadir$sim, " nadir values...", sep = "")
      message(m)
      
      # Set variables
      nadirDIST <- vector()
      nadirCORR <- vector()
      
      # Begin the simulation
      for (i in 1:nadir$sim) {
        set.seed(nadir$seeds[i])
        pts <- sample(1:n.candi, n.pts)
        sm <- covars[pts, ]
        scm <- .corCORR(obj = sm, covars.type = covars.type)
        nadirDIST[i] <- .objDIST(
          sm = sm, n.pts = n.pts, n.cov = n.cov, pop.prop = pop.prop, covars.type = covars.type)
        nadirCORR[i] <- .objCORR(scm = scm, pcm = pcm)
      }
      
      # Prepare output
      res <- list(DIST = mean(nadirDIST), CORR = mean(nadirCORR))
      
    } else {
      
      # User-defined nadir values
      if (!is.null(nadir$user)) { 
        res <- list(DIST = nadir$user$DIST, CORR = nadir$user$CORR)
        
      } else {
        if (!is.null(nadir$abs)) { 
          message("sorry but the nadir point cannot be calculated")
        }
      }
    }
    return (res)
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE ############################
# This function is used to calculate the criterion value of ACDC.
# It calculates, scales, weights, and aggregates the objective function values.
# Scaling is done using the upper-lower bound approach.
# Aggregation is done using the weighted sum method.
.objACDC <-
  function (sm, n.cov, nadir, weights, n.pts, utopia, pcm, scm, covars.type, pop.prop) {
    
    # DIST
    obj_dist <- .objDIST(
      sm = sm, n.pts = n.pts, n.cov = n.cov, pop.prop = pop.prop, covars.type = covars.type)
    obj_dist <- (obj_dist - utopia$DIST) / (nadir$DIST - utopia$DIST)
    obj_dist <- obj_dist * weights$DIST
    
    # CORR
    obj_cor <- .objCORR(scm = scm, pcm = pcm)
    obj_cor <- (obj_cor - utopia$CORR) / (nadir$CORR - utopia$CORR)
    obj_cor <- obj_cor * weights$CORR
    
    # Prepare output
    res <- data.frame(obj = obj_dist + obj_cor, CORR = obj_cor, DIST = obj_dist)
    return (res)
  }
# INTERNAL FUNCTION - PREPARE THE UTOPIA POINT #################################
.utopiaACDC <-
  function (utopia) {
    
    if (!is.null(unlist(utopia$user))) {
      list(CORR = utopia$user$CORR, DIST = utopia$user$DIST)
      
    } else {
      message("sorry but the utopia point cannot be calculated")
    }
  }
# CALCULATE OBJECTIVE FUNCTION VALUE ###########################################
#' @rdname optimACDC
#' @export
objACDC <-
  function (points, candi,
            # DIST and CORR
            covars, strata.type = "area", use.coords = FALSE,
            # MOOP
            weights = list(CORR = 0.5, DIST = 0.5),
            nadir = list(sim = NULL, seeds = NULL, user = NULL, abs = NULL),
            utopia = list(user = NULL, abs = NULL)) {
    
    # Check arguments
    check <- .optimACDCcheck(
      candi = candi, covars = covars, use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare 'covars' and create the starting sample matrix 'sm'
    eval(.prepare_acdc_covars())
    
    # Compute base data
    pcm <- .corCORR(obj = covars, covars.type = covars.type)
    scm <- .corCORR(obj = sm, covars.type = covars.type)
    pop_prop <- .strataACDC(
      n.pts = n_pts, strata.type = strata.type, covars = covars, covars.type = covars.type)
    nadir <- .nadirACDC(
      n.pts = n_pts, n.cov = n_cov, n.candi = n_candi, pcm = pcm, nadir = nadir, candi = candi, 
      covars = covars, pop.prop = pop_prop, covars.type = covars.type)
    utopia <- .utopiaACDC(utopia = utopia)
    
    # Compute the energy state
    energy <- .objACDC(
      sm = sm, pop.prop = pop_prop, covars.type = covars.type, weights = weights, pcm = pcm, scm = scm, 
      n.pts = n_pts, n.cov = n_cov, utopia = utopia, nadir = nadir)
    return (energy)
  }
