#' Optimization of sample configurations for variogram and spatial trend identification and estimation, and 
#' for spatial interpolation
#' 
#' Optimize a sample configuration for variogram and spatial trend identification and estimation, and for 
#' spatial interpolation. An utility function \emph{U} is defined so that the sample points cover, extend 
#' over, spread over, \bold{SPAN} the feature, variogram and geographic spaces. The utility function is 
#' obtained aggregating four objective functions: \bold{CORR}, \bold{DIST}, \bold{PPL}, and \bold{MSSD}.
#' 
#' @param x.max,x.min,y.max,y.min Numeric value defining the minimum and maximum quantity of random noise to 
#' be added to the projected x- and y-coordinates. The minimum quantity should be equal to, at least, the 
#' minimum distance between two neighbouring candidate locations. The units are the same as of the projected 
#' x- and y-coordinates. If missing, they are estimated from \code{candi}.
#' 
#' @inheritParams spJitter
#' @template spSANN_doc
#' @template ACDC_doc
#' @template MOOP_doc
#' @template PPL_doc
#' @template spJitter_doc
#' 
#' @details 
#' Visit the help pages of \code{\link[spsann]{optimCORR}}, \code{\link[spsann]{optimDIST}},
#' \code{\link[spsann]{optimPPL}}, and \code{\link[spsann]{optimMSSD}} to see the details of the objective
#' functions that compose \bold{SPAN}. 
#' 
#' @return
#' \code{optimSPAN} returns an object of class \code{OptimizedSampleConfiguration}: the optimized sample
#' configuration with details about the optimization.
#' 
#' \code{objSPAN} returns a numeric value: the energy state of the sample configuration -- the objective
#' function value.
#' 
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[spsann]{optimCORR}}, \code{\link[spsann]{optimDIST}}, \code{\link[spsann]{optimPPL}},
#' \code{\link[spsann]{optimMSSD}}
#' @aliases optimSPAN objSPAN SPAN
#' @export
#' @examples
#' \dontrun{
#' # This example takes more than 5 seconds to run!
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' nadir <- list(sim = 10, seeds = 1:10)
#' utopia <- list(user = list(DIST = 0, CORR = 0, PPL = 0, MSSD = 0))
#' covars <- meuse.grid[, 5]
#' schedule <- scheduleSPSANN(chains = 1, initial.temperature = 1,
#'                            x.max = 1540, y.max = 2060, x.min = 0, 
#'                            y.min = 0, cellsize = 40)
#' set.seed(2001)
#' res <- optimSPAN(points = 10, candi = candi, covars = covars, nadir = nadir,
#'                  use.coords = TRUE, utopia = utopia, schedule = schedule)
#' objSPSANN(res) -
#'   objSPAN(points = res, candi = candi, covars = covars, nadir = nadir,
#'           use.coords = TRUE, utopia = utopia)
#' }
# MAIN FUNCTION ################################################################
optimSPAN <-
  function(points, candi,
           # DIST and CORR
           covars, strata.type = "area", use.coords = FALSE,
           # PPL
           lags = 7, lags.type = "exponential", lags.base = 2, cutoff, 
           criterion = "distribution", distri, pairs = FALSE,
           # SPSANN
           schedule = scheduleSPSANN(), plotit = FALSE, track = FALSE,
           boundary, progress = "txt", verbose = FALSE,
           # MOOP
           weights = list(CORR = 1/6, DIST = 1/6, PPL = 1/3, MSSD = 1/3),
           nadir = list(sim = NULL, seeds = NULL, user = NULL, abs = NULL),
           utopia = list(user = NULL, abs = NULL)) {
    
    # Objective function name
    objective <- "SPAN"
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Check other arguments
    check <- .checkPPL(lags = lags, lags.type = lags.type, pairs = pairs,
                       lags.base = lags.base, cutoff = cutoff, 
                       criterion = criterion, distri = distri, fun = "optimPPL")
    if (!is.null(check)) stop(check, call. = FALSE)
    check <- .optimACDCcheck(candi = candi, covars = covars, 
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop(check, call. = FALSE)
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Prepare 'covars' and create the starting sample matrix 'sm'
    eval(.prepare_acdc_covars())
    
    # Base data
    # CORR
    pcm <- .corCORR(obj = covars, covars.type = covars.type)
    scm <- .corCORR(obj = sm, covars.type = covars.type)
    # DIST
    pop_prop <- .strataACDC(n.pts = n_pts, strata.type = strata.type, 
                            covars = covars, covars.type = covars.type)
    # PPL
    cutoff <- .cutoffPPL(cutoff = cutoff, x.max = x.max, y.max = y.max)
    lags <- .lagsPPL(lags = lags, lags.type = lags.type, cutoff = cutoff, 
                     lags.base = lags.base)
    n_lags <- length(lags) - 1
    dm_ppl <- SpatialTools::dist1(conf0[, 2:3])
    ppl <- .getPPL(lags = lags, n.lags = n_lags, dist.mat = dm_ppl, 
                   pairs = pairs)
    distri <- .distriPPL(n.lags = n_lags, n.pts = n_pts, criterion = criterion,
                         distri = distri, pairs = pairs)
    # MSSD
    dm_mssd <- SpatialTools::dist2(candi[, 2:3], conf0[, 2:3])
    
    # Nadir and utopia points
    nadir <- .nadirSPAN(n.pts = n_pts, n.cov = n_cov, n.candi = n_candi, 
                        nadir = nadir, candi = candi, covars = covars, 
                        pcm = pcm, pop.prop = pop_prop, lags = lags,
                        covars.type = covars.type, n.lags = n_lags, 
                        pairs = pairs, distri = distri, criterion = criterion)
    utopia <- .utopiaSPAN(utopia = utopia)
    
    # Energy state
    energy0 <- .objSPAN(sm = sm, n.cov = n_cov, nadir = nadir, utopia = utopia,
                        weights = weights, n.pts = n_pts, pcm = pcm, scm = scm,
                        covars.type = covars.type, pop.prop = pop_prop, 
                        ppl = ppl, n.lags = n_lags, criterion = criterion, 
                        distri = distri, pairs = pairs, dm.mssd = dm_mssd)
    
    # Other settings for the simulated annealing algorithm
    # DIST and CORR
    old_sm <- sm
    new_sm <- sm
    best_sm <- sm
    old_scm <- scm
    best_scm <- scm
    # PPL
    old_dm_ppl <- dm_ppl
    best_dm_ppl <- dm_ppl
    # MSSD
    old_dm_mssd <- dm_mssd
    best_dm_mssd <- dm_mssd
    # other
    old_energy <- energy0
    best_energy <- 
      data.frame(obj = Inf, CORR = Inf, DIST = Inf, PPL = Inf, MSSD = Inf)
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
          
          # Update base data and energy state
          # DIST and CORR
          new_sm[wp, ] <- covars[new_conf[wp, 1], ]
          new_scm <- .corCORR(obj = new_sm, covars.type = covars.type)
          # PPL
          new_dm_ppl <- 
            .updatePPLCpp(x = new_conf[, 2:3], dm = old_dm_ppl, idx = wp)
          ppl <- .getPPL(lags = lags, n.lags = n_lags, dist.mat = new_dm_ppl, 
                         pairs = pairs)
          # MSSD
          x2 <- matrix(new_conf[wp, 2:3], nrow = 1)
          new_dm_mssd <- .updateMSSDCpp(x1 = candi[, 2:3], x2 = x2, 
                                        dm = old_dm_mssd, idx = wp)
          # Energy state
          new_energy <- 
            .objSPAN(sm = new_sm, n.cov = n_cov, nadir = nadir, 
                     weights = weights, n.pts = n_pts, utopia = utopia, 
                     pcm = pcm, scm = new_scm, covars.type = covars.type, 
                     pop.prop = pop_prop, ppl = ppl, n.lags = n_lags,
                     criterion = criterion, distri = distri, 
                     pairs = pairs, dm.mssd = new_dm_mssd)
          
          # Avoid the following error:
          # Error in if (new_energy[1] <= old_energy[1]) { : 
          #   missing value where TRUE/FALSE needed
          # Source: http://stackoverflow.com/a/7355280/3365410
          # ASR: The reason for the error is unknown to me.
          if (is.na(new_energy[[1]])) {
            new_energy <- old_energy
            new_conf <- old_conf
            # DIST and CORR
            new_sm <- old_sm
            new_scm <- old_scm
            # PPL
            new_dm_ppl <- old_dm_ppl
            # MSSD
            new_dm_mssd <- old_dm_mssd
          }
          
          # Evaluate the new system configuration
          accept <- .acceptSPSANN(old_energy[[1]], new_energy[[1]], actual_temp)
          if (accept) {
            old_conf <- new_conf
            old_energy <- new_energy
            # DIST and CORR
            old_sm <- new_sm
            old_scm <- new_scm
            # PPL
            old_dm_ppl <- new_dm_ppl
            # MSSD
            old_dm_mssd <- new_dm_mssd
            n_accept <- n_accept + 1
          } else {
            new_energy <- old_energy
            new_conf <- old_conf
            # DIST and CORR
            new_sm <- old_sm
            new_scm <- old_scm
            # PPL
            new_dm_ppl <- old_dm_ppl
            # MSSD
            new_dm_mssd <- old_dm_mssd
          }
          if (track) energies[k, ] <- new_energy
          
          # Record best energy state
          if (new_energy[[1]] < best_energy[[1]] / 1.0000001) {
            best_k <- k
            best_conf <- new_conf
            best_energy <- new_energy
            best_old_energy <- old_energy
            old_conf <- old_conf
            # DIST and CORR
            best_sm <- new_sm
            best_old_sm <- old_sm
            best_scm <- new_scm
            best_old_scm <- old_scm
            # PPL
            best_dm_ppl <- new_dm_ppl
            best_old_dm_ppl <- old_dm_ppl
            # MSSD
            best_dm_mssd <- new_dm_mssd
            best_old_dm_mssd <- old_dm_mssd
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
            # DIST and CORR
            # new_sm <- best_sm
            # new_scm <- best_scm
            # old_sm <- best_old_sm
            # old_scm <- best_old_scm
            # PPL
            # new_dm_ppl <- best_dm_ppl
            # old_dm_ppl <- best_old_dm_ppl
            # MSSD
            # new_dm_mssd <- best_dm_mssd
            # old_dm_mssd <- best_old_dm_mssd
            # no_change <- 0
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
  }
# CALCULATE OBJECTIVE FUNCTION VALUE ###########################################
#' @rdname optimSPAN
#' @export
objSPAN <-
  function(points, candi, 
           # DIST and CORR
           covars, strata.type = "area", use.coords = FALSE,
           # PPL
           lags = 7, lags.type = "exponential", lags.base = 2, cutoff, 
           criterion = "distribution", distri, pairs = FALSE,
           # SPSANN
           x.max, x.min, y.max, y.min,
           # MOOP
           weights = list(CORR = 1/6, DIST = 1/6, PPL = 1/3, MSSD = 1/3),
           nadir = list(sim = NULL, seeds = NULL, user = NULL, abs = NULL),
           utopia = list(user = NULL, abs = NULL)) {
    
    # Check other arguments
    check <- .checkPPL(
      lags = lags, lags.type = lags.type, pairs = pairs, lags.base = lags.base, cutoff = cutoff, 
      criterion = criterion, distri = distri, fun = "optimPPL")
    if (!is.null(check)) stop(check, call. = FALSE)
    check <- .optimACDCcheck(
      candi = candi, covars = covars, use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop(check, call. = FALSE)
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    if (missing(cutoff)) {
      schedule <- scheduleSPSANN()
      eval(.prepare_jittering())
    }
    
    # Prepare 'covars' and create the starting sample matrix 'sm'
    eval(.prepare_acdc_covars())
    
    # Base data
    # CORR
    pcm <- .corCORR(obj = covars, covars.type = covars.type)
    scm <- .corCORR(obj = sm, covars.type = covars.type)
    # DIST
    pop_prop <- .strataACDC(
      n.pts = n_pts, strata.type = strata.type, covars = covars, covars.type = covars.type)
    # PPL
    cutoff <- .cutoffPPL(cutoff = cutoff, x.max = x.max, y.max = y.max)
    lags <- .lagsPPL(lags = lags, lags.type = lags.type, cutoff = cutoff, lags.base = lags.base)
    n_lags <- length(lags) - 1
    dm_ppl <- SpatialTools::dist1(conf0[, 2:3])
    ppl <- .getPPL(lags = lags, n.lags = n_lags, dist.mat = dm_ppl, pairs = pairs)
    distri <- .distriPPL(
      n.lags = n_lags, n.pts = n_pts, criterion = criterion, distri = distri, pairs = pairs)
    # MSSD
    dm_mssd <- SpatialTools::dist2(candi[, 2:3], conf0[, 2:3])
    
    # Nadir and utopia points
    nadir <- .nadirSPAN(
      n.pts = n_pts, n.cov = n_cov, n.candi = n_candi, nadir = nadir, candi = candi, covars = covars, 
      pcm = pcm, pop.prop = pop_prop, lags = lags, covars.type = covars.type, n.lags = n_lags, pairs = pairs,
      distri = distri, criterion = criterion)
    utopia <- .utopiaSPAN(utopia = utopia)
    
    # Energy state
    res <- .objSPAN(
      sm = sm, n.cov = n_cov, nadir = nadir, utopia = utopia, weights = weights, n.pts = n_pts, pcm = pcm, 
      scm = scm, covars.type = covars.type, pop.prop = pop_prop, ppl = ppl, n.lags = n_lags, 
      criterion = criterion, distri = distri, pairs = pairs, dm.mssd = dm_mssd)
    
    # Output
    return(res)
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE ############################
# This function is used to calculate the criterion value of SPAN.
# It calculates, scales, weights, and aggregates the objective function values.
# Scaling is done using the upper-lower bound approach.
# Aggregation is done using the weighted sum method.
.objSPAN <-
  function(sm, n.cov, nadir, weights, n.pts, utopia, pcm, scm, covars.type,
           pop.prop, ppl, n.lags, criterion, distri, pairs, dm.mssd) {
    
    # DIST
    obj_dist <- .objDIST(sm = sm, n.pts = n.pts, n.cov = n.cov, 
                         pop.prop = pop.prop, covars.type = covars.type)
    obj_dist <- (obj_dist - utopia$DIST) / (nadir$DIST - utopia$DIST)
    obj_dist <- obj_dist * weights$DIST
    
    # CORR
    obj_cor <- .objCORR(scm = scm, pcm = pcm)
    obj_cor <- (obj_cor - utopia$CORR) / (nadir$CORR - utopia$CORR)
    obj_cor <- obj_cor * weights$CORR
    
    # PPL
    obj_ppl <- .objPPL(ppl = ppl, n.lags = n.lags, n.pts = n.pts, 
                       criterion = criterion, distri = distri, pairs = pairs)
    obj_ppl <- (obj_ppl - utopia$PPL) / (nadir$PPL - utopia$PPL)
    obj_ppl <- obj_ppl * weights$PPL
    
    # MSSD
    obj_mssd <- .objMSSD(x = dm.mssd)
    obj_mssd <- (obj_mssd - utopia$MSSD) / (nadir$MSSD - utopia$MSSD)
    obj_mssd <- obj_mssd * weights$MSSD
    
    # Prepare output
    res <- data.frame(obj = obj_dist + obj_cor + obj_ppl + obj_mssd, 
                      CORR = obj_cor, DIST = obj_dist, PPL = obj_ppl,
                      MSSD = obj_mssd)
    return(res)
  }
# INTERNAL FUNCTION - COMPUTE THE NADIR VALUE ##################################
.nadirSPAN <-
  function(n.pts, n.cov, n.candi, nadir, candi, covars, pcm, pop.prop, 
           covars.type, lags, n.lags, pairs, distri, criterion) {
    
    # Simulate the nadir point
    if (!is.null(nadir$sim) && !is.null(nadir$seeds)) { 
      m <- paste("simulating ", nadir$sim, " nadir values...", sep = "")
      message(m)
      
      # Set variables
      nadirDIST <- vector()
      nadirCORR <- vector()
      nadirPPL <- vector()
      nadirMSSD <- vector()
      
      # Begin the simulation
      for (i in 1:nadir$sim) {
        set.seed(nadir$seeds[i])
        pts <- sample(1:n.candi, n.pts)
        sm <- covars[pts, ]
        # CORR
        scm <- .corCORR(obj = sm, covars.type = covars.type)
        nadirCORR[i] <- .objCORR(scm = scm, pcm = pcm)
        # DIST
        nadirDIST[i] <- .objDIST(sm = sm, n.pts = n.pts, n.cov = n.cov, 
                                 pop.prop = pop.prop, covars.type = covars.type)
        # PPL
        dm <- SpatialTools::dist1(candi[pts, 2:3])
        ppl <- .getPPL(lags = lags, n.lags = n.lags, dist.mat = dm, 
                       pairs = pairs)
        nadirPPL[i] <- .objPPL(ppl = ppl, n.lags = n.lags, n.pts = n.pts,
                               criterion = criterion, distri = distri, 
                               pairs = pairs)
        # MSSD
        dm <- SpatialTools::dist2(candi[, 2:3], candi[pts, 2:3])
        nadirMSSD[i] <- .objMSSD(x = dm)
      }
      
      # Prepare output
      res <- list(DIST = mean(nadirDIST), CORR = mean(nadirCORR), 
                  PPL = mean(nadirPPL), MSSD = mean(nadirMSSD))
      
    } else {
      
      # User-defined nadir values
      if (!is.null(nadir$user)) { 
        res <- list(DIST = nadir$user$DIST, CORR = nadir$user$CORR,
                    PPL = nadir$user$PPL, MSSD = nadir$user$MSSD)
        
      } else {
        if (!is.null(nadir$abs)) { 
          message("sorry but the nadir point cannot be calculated")
        }
      }
    }
    return(res)
  }
# INTERNAL FUNCTION - COMPUTE THE UTOPIA POINT #################################
.utopiaSPAN <-
  function(utopia) {
    
    if (!is.null(unlist(utopia$user))) {
      list(CORR = utopia$user$CORR, DIST = utopia$user$DIST,
           PPL = utopia$user$PPL, MSSD = utopia$user$MSSD)
      
    } else {
      message("sorry but the utopia point cannot be calculated")
    }
  }
