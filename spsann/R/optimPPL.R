#' Optimization of sample configurations for variogram identification and estimation
#'
#' Optimize a sample configuration for variogram identification and estimation. A criterion is defined so that
#' the optimized sample configuration has a given number of points or point-pairs contributing to each 
#' lag-distance class (\bold{PPL}).
#' 
#' @inheritParams spJitter
#' @template spSANN_doc
#' @template PPL_doc
#' @template spJitter_doc
#' 
#' @details 
#' \subsection{Lag-distance classes}{
#' Two types of lag-distance classes can be created by default. The first are evenly spaced lags 
#' (\code{lags.type = "equidistant"}). They are created by simply dividing the distance interval from 0.0001 
#' to \code{cutoff} by the required number of lags. The minimum value of 0.0001 guarantees that a point does 
#' not form a pair with itself. The second type of lags is defined by exponential spacings 
#' (\code{lags.type = "exponential"}). The spacings are defined by the base \eqn{b} of the exponential 
#' expression \eqn{b^n}, where \eqn{n} is the required number of lags. The base is defined using the argument 
#' \code{lags.base}. See \code{\link[pedometrics]{vgmLags}} for other details.
#'
#' Using the default uniform distribution means that the number of point-pairs per lag-distance class 
#' (\code{pairs = TRUE}) is equal to \eqn{n \times (n - 1) / (2 \times lag)}, where \eqn{n} is the total 
#' number of points and \eqn{lag} is the number of lags. If \code{pairs = FALSE}, then it means that the 
#' number of points per lag is equal to the total number of points. This is the same as expecting that each 
#' point contributes to every lag. Distributions other than the available options can be easily implemented 
#' changing the arguments \code{lags} and \code{distri}.
#' 
#' There are two optimizing criteria implemented. The first is called using \code{criterion = "distribution"}
#' and is used to minimize the sum of the absolute differences between a pre-specified distribution and the
#' observed distribution of points or point-pairs per lag-distance class. The second criterion is called using
#' \code{criterion = "minimum"}. It corresponds to maximizing the minimum number of points or point-pairs 
#' observed over all lag-distance classes.
#' }
#' 
#' @return
#' \code{optimPPL} returns an object of class \code{OptimizedSampleConfiguration}: the optimized sample
#' configuration with details about the optimization.
#'
#' \code{objPPL} returns a numeric value: the energy state of the sample configuration -- the objective 
#' function value.
#'
#' \code{countPPL} returns a data.frame with three columns: a) the lower and b) upper limits of each 
#' lag-distance class, and c) the number of points or point-pairs per lag-distance class.
#'
#' @references
#' Bresler, E.; Green, R. E. \emph{Soil parameters and sampling scheme for characterizing soil hydraulic
#' properties of a watershed}. Honolulu: University of Hawaii at Manoa, p. 42, 1982.
#'
#' Pettitt, A. N.; McBratney, A. B. Sampling designs for estimating spatial variance components. 
#' \emph{Applied Statistics}. v. 42, p. 185, 1993.
#'
#' Russo, D. Design of an optimal sampling network for estimating the variogram. \emph{Soil Science Society of
#' America Journal}. v. 48, p. 708-716, 1984.
#'
#' Truong, P. N.; Heuvelink, G. B. M.; Gosling, J. P. Web-based tool for expert elicitation of the variogram.
#' \emph{Computers and Geosciences}. v. 51, p.
#' 390-399, 2013.
#'
#' Warrick, A. W.; Myers, D. E. Optimization of sampling locations for variogram calculations. \emph{Water 
#' Resources Research}. v. 23, p. 496-500, 1987.
#'
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @aliases optimPPL countPPL objPPL PPL
#' @export
#' @examples
#' \dontrun{
#' # This example takes more than 5 seconds
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' schedule <- scheduleSPSANN(chains = 1, initial.temperature = 30,
#'                            x.max = 1540, y.max = 2060, x.min = 0, 
#'                            y.min = 0, cellsize = 40)
#' set.seed(2001)
#' res <- optimPPL(points = 10, candi = candi, schedule = schedule)
#' objSPSANN(res) - objPPL(points = res, candi = candi)
#' countPPL(points = res, candi = candi)
#' }
# FUNCTION - MAIN ##############################################################
optimPPL <-
  function (points, candi,
            # PPL
            lags = 7, lags.type = "exponential", lags.base = 2, cutoff, 
            criterion = "distribution", distri, pairs = FALSE,
            # SPSANN
            schedule = scheduleSPSANN(), plotit = FALSE, track = FALSE, 
            boundary, progress = "txt", verbose = FALSE) {
    
    # Objective function name
    objective <- "PPL"
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Check other arguments
    check <- .checkPPL(
      lags = lags, lags.type = lags.type, pairs = pairs, lags.base = lags.base, cutoff = cutoff, 
      criterion = criterion, distri = distri, fun = "optimPPL")
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Prepare cutoff and lags
    cutoff <- .cutoffPPL(cutoff = cutoff, x.max = x.max, y.max = y.max)
    lags <- .lagsPPL(lags = lags, lags.type = lags.type, cutoff = cutoff, lags.base = lags.base)
    n_lags <- length(lags) - 1
    
    # Initial energy state: points or point-pairs
    dm <- SpatialTools::dist1(conf0[, 2:3])
    ppl <- .getPPL(
      lags = lags, n.lags = n_lags, dist.mat = dm, pairs = pairs, n.pts = n_pts)
    distri <- .distriPPL(n.lags = n_lags, n.pts = n_pts, criterion = criterion,
                         distri = distri, pairs = pairs)
    energy0 <- data.frame(
      obj = .objPPL(ppl = ppl, n.lags = n_lags, n.pts = n_pts,
                    criterion = criterion, distri = distri, pairs = pairs))
    
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
          
          # Update the distance matrix using a Cpp function
          new_dm <- .updatePPLCpp(x = new_conf[, 2:3], dm = old_dm, idx = wp)
          
          # Recalculate the full distance matrix
          #new_dm <- SpatialTools::dist1(coords = new_conf[, 2:3])
          
          # Update the distance matrix in R
          #x2 <- matrix(new_conf[wp, 2:3], nrow = 1)
          #x2 <- SpatialTools::dist2(coords = new_conf[, 2:3], coords2 = x2)
          #new_dm <- old_dm
          #new_dm[wp, ] <- x2
          #new_dm[, wp] <- x2
          
          # Update the energy state: points or point-pairs?
          ppl <- .getPPL(
            lags = lags, n.lags = n_lags, dist.mat = new_dm, pairs = pairs, 
            n.pts = n_pts)
          new_energy <- data.frame(
            obj = .objPPL(n.lags = n_lags, n.pts = n_pts, pairs = pairs,
                          criterion = criterion, distri = distri, ppl = ppl))
          
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
# FUNCTION - CALCULATE THE CRITERION VALUE #####################################
#' @rdname optimPPL
#' @export
objPPL <-
  function (points, candi, 
            # PPL
            lags = 7, lags.type = "exponential", lags.base = 2, cutoff, distri,
            criterion = "distribution", pairs = FALSE,
            # SPSANN
            x.max, x.min, y.max, y.min) {
    
    # Check arguments
    check <- .checkPPL(lags = lags, lags.type = lags.type, pairs = pairs, 
                       lags.base = lags.base, cutoff = cutoff,
                       criterion = criterion, distri = distri, fun = "objPPL")
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare cutoff and lags
    if (missing(cutoff) || length(lags) == 1) {
      schedule <- scheduleSPSANN()
      eval(.prepare_jittering())
      cutoff <- .cutoffPPL(cutoff = cutoff, x.max = x.max, y.max = y.max)
      lags <- .lagsPPL(lags = lags, lags.type = lags.type, cutoff = cutoff, 
                       lags.base = lags.base)
    }
    n_lags <- length(lags) - 1
    
    # Initial energy state: points or point-pairs
    dm <- SpatialTools::dist1(points[, 2:3])
    ppl <- .getPPL(
      lags = lags, n.lags = n_lags, dist.mat = dm, pairs = pairs, n.pts = n_pts)
    distri <- .distriPPL(n.lags = n_lags, n.pts = n_pts, criterion = criterion,
                         distri = distri, pairs = pairs)
    res <- .objPPL(ppl = ppl, n.lags = n_lags, n.pts = n_pts, pairs = pairs,
                   criterion = criterion, distri = distri)
    
    # Output
    return (res)
  }
# FUNCTION - POINTS PER LAG-DISTANCE CLASS #####################################
#' @rdname optimPPL
#' @export
countPPL <-
  function (points, candi,
            # PPL
            lags = 7, lags.type = "exponential", lags.base = 2, cutoff, pairs = FALSE,
            # SPSANN
            x.max, x.min, y.max, y.min) {
    
    # Check arguments
    check <- .checkPPL(lags = lags, lags.type = lags.type, pairs = pairs,
                       lags.base = lags.base, cutoff = cutoff, fun = "countPPL")
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare cutoff and lags
    if (missing(cutoff)) {
      schedule <- scheduleSPSANN()
      eval(.prepare_jittering())
    }
    cutoff <- .cutoffPPL(cutoff = cutoff, x.max = x.max, y.max = y.max)
    lags <- .lagsPPL(lags = lags, lags.type = lags.type, cutoff = cutoff, 
                     lags.base = lags.base)
    n_lags <- length(lags) - 1
    
    # Distance matrix and counts
    dm <- SpatialTools::dist1(points[, 2:3])
    res <- .getPPL(
      lags = lags, n.lags = n_lags, dist.mat = dm, pairs = pairs, n.pts = n_pts)
    res <- data.frame(lag.lower = lags[-length(lags)], lag.upper = lags[-1],
                      ppl = res)
    return (res)
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
.checkPPL <-
  function (lags, lags.type, lags.base, cutoff, criterion, distri, pairs, fun) {
    
    # Argument 'pairs'
    if (!is.logical(pairs)) {
      res <- c("'pairs' must be a logical value")
      return (res)
    }
    
    # Arguments 'lags' and 'cutoff'
    if (length(lags) == 1) {
      #       if (fun != "optimPPL") {
      #         if (missing(cutoff)) {
      #           res <- c("'cutoff' is mandatory if the lag intervals are not set")
      #           return (res)
      #         }
      #       }
    } else {
      if (!missing(cutoff)) {
        res <- c("'cutoff' cannot be used when the lag intervals are set")
        return (res)
      }
      if (lags[1] == 0) {
        res <- c("lowest lag value must be larger than zero")
        return (res)
      }
    }
    
    # Argument 'lags.type'
    lt <- c("equidistant", "exponential")
    lt <- is.na(any(match(lt, lags.type)))
    if (lt) {
      res <- paste("'lags.type = ", lags.type, "' is not supported", sep = "")
      return (res)
    }
    
    # Argument 'lags.base'
    if (!is.numeric(lags.base) || length(lags.base) > 1) {
      res <- c("'lags.base' must be a numeric value")
      return (res)
    }
    
    if (fun != "countPPL") {
      
      # Argument 'criterion'
      cr <- is.na(any(match(criterion, c("distribution", "minimum"))))
      if (cr) {
        res <- paste("'criterion = ", criterion, "' is not supported", sep = "")
        return (res)
      }
      
      # Argument 'distri'
      if (!missing(distri)) {
        if (!is.numeric(distri)) {
          res <- c("'distri' must be a numeric vector")
          return (res)
        }
        if (length(lags) == 1) {
          if (length(distri) != lags) {
            res <- paste("'distri' must be of length ", lags, sep = "")
            return (res)
          }
        }
        if (length(lags) > 2) {
          nl <- length(lags) - 1
          if (length(distri) != nl) {
            res <- paste("'distri' must be of length ", nl, sep = "")
            return (res)
          }
        }
      }
    }
  }
# INTERNAL FUNCTION - COMPUTE THE DISTRIBUTION #################################
# If required and missing, the wanted distribution of poins (point-pairs) per
# lag distance class is computed internally.
.distriPPL <-
  function (n.lags, n.pts, criterion, distri, pairs) {
    
    if (criterion == "distribution" && missing(distri)) {
      
      if (pairs) { # Point-pairs per lag
        distri <- rep(n.pts * (n.pts - 1) / (2 * n.lags), n.lags)
        
      } else { # Point per lag
        distri <- rep(n.pts, n.lags)
      }
    }
    
    return (distri)
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE ############################
.objPPL <-
  function (ppl, n.lags, n.pts, criterion, distri, pairs) {
    
    if (pairs) { # Point-pairs per lag
      
      if (criterion == "distribution") {
        res <- sum(abs(distri - ppl))
        
      } else { # minimum
        a <- n.pts * (n.pts - 1) / (2 * n.lags)
        res <- a / (min(ppl) + 1)
      }
      
    } else { # Points per lag
      
      if (criterion == "distribution") {
        res <- sum(distri - ppl)
        
      } else { # minimum
        res <- n.pts / (min(ppl) + 1)
      }
    }
    
    # Output
    return (res)
  }
# INTERNAL FUNCTION - NUMBER OF POINTS (POINT-PAIRS) PER LAG-DISTANCE CLASS ####
# Note:
# It is 3 times faster to use the for loop with function 'which()' than using
# 'apply()' with functions 'table()' and 'cut()'.
# apply(X = dist.mat, 1, FUN = function (X) table(cut(X, breaks = lags)))
# apply(X = ppl, 1, FUN = function (X) sum(X != 0))
.getPPL <-
  function (lags, n.lags, dist.mat, pairs, n.pts) {
    
    # ppl <- vector()
    
    if (pairs) { # Point-pairs per lag
      # for (i in 1:n.lags) {
        # n <- which(dist.mat > lags[i] & dist.mat <= lags[i + 1])
        # ppl[i] <- length(n)
      # }
      ppl <- diff(sapply(1:length(lags), function (i) length(which(dist.mat <= lags[i]))) - n.pts) * 0.5
    } else { # Points per lag
      # for (i in 1:n.lags) {
        # n <- which(
          # dist.mat > lags[i] & dist.mat <= lags[i + 1], arr.ind = TRUE)
        # ppl[i] <- length(unique(c(n)))
      # }
      ppl <- sapply(1:n.lags, function (i)
        length(unique(c(which(dist.mat > lags[i] & dist.mat <= lags[i + 1], arr.ind = TRUE)))))
    }
    
    return (ppl)
  }
# INTERNAL FUNCTION - BREAKS OF THE LAG-DISTANCE CLASSES #######################
.lagsPPL <-
  function (lags, lags.type, cutoff, lags.base) {
    
    if (length(lags) == 1) {
      
      if (lags.type == "equidistant") {
        lags <- seq(0.0001, cutoff, length.out = lags + 1)
      }
      
      if (lags.type == "exponential") {
        idx <- lags.base ^ c(1:lags - 1)
        lags <- c(0.0001, rev(cutoff / idx))
      }
    }
    
    return (lags)
  }
# INTERNAL FUNCTION - DETERMINE THE CUTOFF #####################################
.cutoffPPL <- 
  function (cutoff, x.max, y.max) {
    
    if (missing(cutoff)) {
      cutoff <- sqrt(x.max ^ 2 + y.max ^ 2)
    }
    
    return (cutoff)
  }
