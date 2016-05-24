#' Optimization of sample configurations for spatial trend identification and estimation (I)
#'
#' Optimize a sample configuration for spatial trend identification and estimation. A criterion is defined so
#' that the sample reproduces the bivariate association/correlation between the covariates (\bold{CORR}).
#'
#' @inheritParams spJitter
#' @template spSANN_doc
#' @template ACDC_doc
#' @template spJitter_doc
#' 
#' @details 
#' \subsection{Association/Correlation between covariates}{
#' The \emph{correlation} between two numeric covariates is measured using the Pearson's \emph{r}, a 
#' descriptive statistic that ranges from \eqn{-1} to \eqn{+1}. This statistic is also known as the linear 
#' correlation coefficient.
#' 
#' When the set of covariates includes factor covariates, all numeric covariates are transformed into factor 
#' covariates. The factor levels are defined using the marginal sampling strata created from one of the two 
#' methods available (equal-area or equal-range strata).
#' 
#' The \emph{association} between two factor covariates is measured using the Cramér's \emph{V}, a descriptive 
#' statistic that ranges from \eqn{0} to \eqn{+1}. The closer to \eqn{+1} the Cramér's \emph{V} is, the 
#' stronger the association between two factor covariates.
#' 
#' The main weakness of using the Cramér's \emph{V} is that, while the Pearson's \emph{r} shows the degree 
#' and direction of the association between two covariates (negative or positive), the Cramér's \emph{V} only
#' measures the degree of association (weak or strong). The effect of replacing the Pearson's \emph{r} with 
#' the Cramér's \emph{V} on the spatial modelling outcome still is poorly understood.
#' }
#' 
#' @return
#' \code{optimCORR} returns an object of class \code{OptimizedSampleConfiguration}: the optimized sample
#' configuration with details about the optimization.
#' 
#' \code{objCORR} returns a numeric value: the energy state of the sample configuration -- the objective
#' function value.
#'
#' @references 
#' Cramér, H. \emph{Mathematical methods of statistics}. Princeton: Princeton University Press, p. 575, 1946.
#' 
#' Everitt, B. S. \emph{The Cambridge dictionary of statistics}. Cambridge: Cambridge University Press, p. 432,
#' 2006.
#'
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[clhs]{clhs}}, \code{\link[pedometrics]{cramer}}, \code{\link[spsann]{optimACDC}}
#' @aliases optimCORR objCORR CORR
#' @export
#' @examples
#' data(meuse.grid, package = "sp")
#' candi <- meuse.grid[1:1000, 1:2]
#' covars <- meuse.grid[1:1000, 5]
#' schedule <- scheduleSPSANN(
#'   initial.temperature = 5, chains = 1, x.max = 1540, y.max = 2060, 
#'   x.min = 0, y.min = 0, cellsize = 40)
#' set.seed(2001)
#' res <- optimCORR(
#'   points = 10, candi = candi, covars = covars, use.coords = TRUE, 
#'   schedule = schedule)
#' objSPSANN(res) - objCORR(
#'   points = res, candi = candi, covars = covars, use.coords = TRUE)
# MAIN FUNCTION ################################################################
optimCORR <-
  function (points, candi,
            # CORR
            covars, strata.type = "area", use.coords = FALSE,
            # SPSANN
            schedule = scheduleSPSANN(), plotit = FALSE, track = FALSE,
            boundary, progress = "txt", verbose = FALSE) {
    
    # Objective function name
    objective <- "CORR"
    
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
    energy0 <- data.frame(obj = .objCORR(scm = scm, pcm = pcm))
    
    # Other settings for the simulated annealing algorithm
    old_sm <- sm
    new_sm <- sm
    best_sm <- sm
    old_scm <- scm
    best_scm <- scm
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
          
          # Update sample and correlation matrices, and energy state
          new_sm[wp, ] <- covars[new_conf[wp, 1], ]
          new_scm <- .corCORR(obj = new_sm, covars.type = covars.type)
          new_energy <- data.frame(obj = .objCORR(scm = new_scm, pcm = pcm))
          
          # Avoid the following error:
          # Error in if (new_energy <= old_energy) { : 
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
# FUNCTION - CALCULATE ENERGY STATE ############################################
#' @rdname optimCORR
#' @export
objCORR <-
  function (points, candi,
            # CORR
            covars, strata.type = "area", use.coords = FALSE) {
    
    # Check other arguments
    check <- .optimACDCcheck(
      candi = candi, covars = covars, use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare 'covars' and create the starting sample matrix 'sm'
    eval(.prepare_acdc_covars())
    
    # Calculate the energy state
    pcm <- .corCORR(obj = covars, covars.type = covars.type)
    scm <- .corCORR(obj = sm, covars.type = covars.type)
    energy <- .objCORR(scm = scm, pcm = pcm)
    
    return (energy)
  }
# INTERNAL FUNCTION - CALCULATE ASSOCIATION/CORRELATION MATRIX ################
# obj: the matrix of covariates, for the population correlation matrix, and the
#      sample matrix, for the sample correlation matrix
# covars.type: the type of covariate
.corCORR <-
  function (obj, covars.type) {
    
    if (covars.type == "numeric") { 
      stats::cor(x = obj, use = "complete.obs")
      
    } else { # Factor covariates
      pedometrics::cramer(x = obj)
    }
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE ############################
# scm: sample correlation matrix
# pcm: population correlation matrix
.objCORR <-
  function (scm, pcm) {
    sum(abs(pcm - scm))
  }
