############################################################################################
# wpp.R
#
# The Witness Protection Program algorithm.
#
# Code by
#
#  - Ricardo Silva (ricardo@stats.ucl.ac.uk)
#  - Robin Evans (robin.evans@stats.ox.ac.uk)
#
# Current version: 28/04/2015
# First version: 31/03/2014

#' @title The Witness Protection Program for Causal Effect Estimation
#' 
#' @description
#' 
#' Perform a search for bounds on the average causal effect (ACE) of a given treatment variable
#' \eqn{X} on a given outcome \eqn{Y}. Bounds are based on finding conditional instrumental variables
#' using the faithfulness assumption relaxed to allow for a moderate degree of unfaithfulness.
#' Candidate models are generated from the method described in \code{\link{covsearch}}.
#' 
#' @references 
#' \url{http://papers.nips.cc/paper/5602-causal-inference-through-a-witness-protection-program}
#' 
#' @param problem a \code{\link{cfx}} problem instance for the ACE of a given treatment \eqn{X} on a given outcome \eqn{Y}.
#' @param epsilons an array of six positions corresponding to the relaxation parameters. In order: 
#'        (1) the maximum difference in the conditional probability of the outcome given everything else, as the witness changes levels;
#'        (2) the maximum difference in the conditional probability of the outcome given everything else, and the conditional distribution
#'            excluding latent variables for the witness set at 0;
#'        (3) the maximum difference in the conditional probability of the outcome given everything else, and the conditional distribution
#'            excluding latent variables for the witness set at 1;
#'        (4) the maximum difference in the conditional probability of the treatment given its causes, and the conditional distribution
#'            excluding latent variables
#'        (5) the maximum ratio between the conditional distribution of the latent variable given the witness and the 
#'            marginal distribution of the latent variable. This has to be greater than or equal to 1;
#'        (6) the minimum ratio between the conditional distribution of the latent variable given the witness and the 
#'            marginal distribution of the latent variable. This has to be in the interval (0, 1].
#' @param max_set maximum size of conditioning set. The cost of the procedure grows exponentially as a function of this, 
#'               so be careful when increasing the default value.
#' @param prior_ind prior probability of an independence.
#' @param prior_table effective sample size hyperparameter of a Dirichlet prior for testing independence with contingency tables.
#' @param cred_calc if \code{TRUE}, compute conditional credible intervals for the ACE of highest scoring model.
#' @param M if necessary to compute (conditional) credible intervals, use Monte Carlo with this number of samples.
#' @param analytical_bounds if \code{cred_calc} is \code{TRUE}, use the analytical method for computing bounds if this is also \code{TRUE}.
#' @param pop_solve if \code{TRUE}, assume we know the population graph in \code{problem} instead of data. Notice that data is
#'                 still used when computing posteriors over bounds.
#' @param verbose if \code{TRUE}, print out more detailed information while running the procedure.
#'
#' @return An object of class \code{wpp} containing the copies of the inputs \code{problem}, \code{epsilons},
#'         \code{prior_ind}, \code{prior_table}, \code{analytical_bounds}, plus the following fields:
#'   \item{\code{w_list}}{a list of arrays/lists, where each \code{w_list$witness[i]} is a witness, each
#'                        \code{w_list$Z[[i]]} is the corresponding admissible set, and each
#'                        \code{w_list$witness_score[i]} is the corresponding score for the witness/admissible set.}
#'   \item{\code{hw}}{witness corresponding to the highest scoring pair.}
#'   \item{\code{hZ}}{array containing admissible set corresponding to the highest scoring pair.}
#'   \item{\code{bounds}}{a two-column matrix where each row corresponds to a different witness/admissible set
#'                        combination, and the two columns correspond to an estimate of the lower bound and upper bound
#'                        as given by the posterior expected value given an inferred causal structure.}
#'   \item{\code{bounds_post}}{a two-column matrix, where rows correspond to different Monte carlo samples, and the two
#'                             columns correspond to lower and upper bounds on the ACE as implied by \code{epsilons} with
#'                             witness \code{hw} and admissible set \code{hZ}.}
#'
#' @details
#' Each pair of witness/admissible set found by \code{covsearch} will generate a corresponding lower bound and upper 
#' bound. The bounds reported in \code{bounds} are based on the posterior expected contingency table implied by 
#' \code{prior_table}, which uses a numerical method to optimize the bounds. Besides these point estimates, posterior
#' distributions on the lower and upper bound for the highest scoring witness/admissible set can also be computed if the
#' flag \code{cred_calc} is set to \code{TRUE}, and reported on \code{bounds_post}. If the option \code{analytical_bounds}
#' is set to \code{FALSE}, the posterior distribution calculation will use the numerical method. It provides tighter
#' bounds, but the computational cost is much higher. Please notice these posteriors are for the bounds conditional on the 
#' given choice of witness and admissible set: uncertainty on this choice is not taken into account.
#' 
#' A complete explanation of the method is given by Silva and Evans (2014, "Causal inference through a witness
#' protection program", \emph{Advances in Neural Information Processing Systems}, 27, 298--306).
#' 
#' Note: messages about numerical problems when calling the bound optimizer are not uncommon and are accounted for
#' within the procedure.
#' 
#' @examples
#' 
#' ## Generate a synthetic problem
#' problem <- simulateWitnessModel(p = 4, q = 4, par_max = 3, M = 200)
#' 
#' ## Calculate true effect for evaluation purposes
#' sol_pop <- covsearch(problem, pop_solve = TRUE)
#' effect_pop <- synthetizeCausalEffect(problem)
#' cat(sprintf("ACE (true) = %1.2f\n", effect_pop$effect_real))
#'
#' ## WPP search (with a small number of Monte Carlo samples)
#' epsilons <- c(0.2, 0.2, 0.2, 0.2, 0.95, 1.05)
#' sol_wpp <- wpp(problem, epsilons, M = 100) 
#' summary(sol_wpp)
#' 
#' @import rcdd
#' @export

wpp <- function(problem, epsilons, max_set = 12, prior_ind = 0.5, prior_table = 10, cred_calc = TRUE,
                M = 1000, analytical_bounds = TRUE, pop_solve = FALSE, verbose = FALSE) {  
  
  if (class(problem) != "cfx") {
    stop("a CausalFX object is necessary")
  }
  
  if (any(epsilons[1:5] <= 0) || any(epsilons[1:5] > 1) || epsilons[6] < 1) {
    stop("relaxation parameters are invalid")
  }
  
  w_list <- covsearch(problem, max_set = max_set, pop_solve = pop_solve, prior_ind = prior_ind, prior_table = prior_table,
                            cred_calc = cred_calc, verbose = verbose)
    
  if (length(w_list$witness) == 0) {
    cat("No solution found\n")
    return(list(w_list = w_list, bounds = list()))
  }
  
  N <- nrow(problem$data)
  bounds <- matrix(0, nrow = length(w_list$witness), ncol = 2)
  
  for (i in seq_along(w_list$witness)) {
    
    w <- w_list$witness[i]
    if (verbose) cat(i, ":: ################################# Bounding using witness", w, "\n")
    Z <- w_list$Z[[i]]
    num_states <- 2^length(Z)    
    
    sub_counts <- dtoc(problem$data[, c(w, problem$X_idx, problem$Y_idx, Z)])    
    param_mean <- wppParamPosteriorExpectation(sub_counts, 1, 2, 3, 3 + seq_along(Z), prior_table = prior_table)
    theta_mean <- param_mean$theta_mean
    P_Z_hat    <- param_mean$Z_mean
    
    for (j in seq_len(num_states)) {
      intervals <- wppIntervalGenerationNumerical(epsilons, theta_samples = theta_mean[[j]])
      bounds[i, ] <- bounds[i, ] + c(intervals[1], intervals[2]) * P_Z_hat[j]
    }
    
    if (verbose) cat("\n")
    
  }
  
  # Purge failed witnessess
  
  if (all(is.na(bounds[, 1] + bounds[, 2]))) {
    cat("No solution found\n")
    return(list(w_list = w_list, bounds = list()))
  }
  
  f_witness <- c()
  f_Z <- list()
  f_witness_score <- matrix(nrow = 0, ncol = 2)
  f_bounds <- matrix(nrow = 0, ncol = 2)
  for (i in 1:length(w_list$witness)) {
    if (!is.na(bounds[i, 1]) && !is.na(bounds[i, 2])) {
      f_witness <- c(f_witness, w_list$witness[i])
      f_Z[[length(f_witness)]] <- w_list$Z[[i]]
      f_witness_score <- rbind(f_witness_score, w_list$witness_score[i,])
      f_bounds <- rbind(f_bounds, bounds[i,])
    }
  }
  f_w_list <- list(witness = f_witness, Z = f_Z, witness_score = f_witness_score)

  idx_highest <- which.max(rowSums(f_witness_score))
  hw <- f_witness[idx_highest]
  hZ <- f_Z[[idx_highest]]
  
  # Credible interval calculation (slow)
  
  if (cred_calc) {
    
    bounds_post <- matrix(rep(0, 2 * M), ncol = 2)
    w <- hw
    Z <- hZ
    num_states <- 2^length(Z)        
    
    if (analytical_bounds) {
    
      sub_counts <- dtoc(problem$data[, c(w, problem$X_idx, problem$Y_idx, Z)])
      
      param_samples <- wppPosteriorSampling(sub_counts, 1, 2, 3, 3 + seq_along(Z), M = M, prior_table = prior_table)
      theta_samples <- param_samples$theta_samples
      P_Z_hat       <- param_samples$Z_samples
            
      for (j in seq_len(num_states)) {
        
        mX1W0 <- theta_samples[[j]]$sXW0; mX1W1 <- theta_samples[[j]]$sXW1
        mY1X0 <- theta_samples[[j]]$sYX0; mY1X1 <- theta_samples[[j]]$sYX1        
        P_YX.W0 <- cbind((1 - mY1X0) * (1 - mX1W0), (1 - mY1X1) * mX1W0, mY1X0 * (1 - mX1W0), mY1X1 * mX1W0)
        P_YX.W1 <- cbind((1 - mY1X0) * (1 - mX1W1), (1 - mY1X1) * mX1W1, mY1X0 * (1 - mX1W1), mY1X1 * mX1W1)
        
        intervals <- wppIntervalGenerationAnalytical(P_YX.W0, P_YX.W1, theta_samples[[j]]$sW, epsilons)
        bounds_post[, 1] <- bounds_post[, 1] + intervals[, 1] * P_Z_hat[, j]
        bounds_post[, 2] <- bounds_post[, 2] + intervals[, 2] * P_Z_hat[, j]
        
      }    
      
    } else {

      cat("Sampling posterior over highest scoring bound...\n")
      
      sub_counts <- dtoc(problem$data[, c(Z, w, problem$X_idx, problem$Y_idx)])
      dim(sub_counts) <- c(num_states, 2, 2, 2)
      if (length(Z) > 0) {  
        Z_totals <- marginTable(sub_counts, 1)
        P_Z_hat <- rdirichlet(M, Z_totals + prior_table / num_states)
      } else {
        P_Z_hat <- matrix(rep(1, M), ncol = 1)
      }
      
      for (j in seq_len(num_states)) {
      
        cat(sprintf("   Level %d out of %d\n", j, num_states))
        N_YX.W0 <- as.numeric(sub_counts[j, 1, , ])
        N_YX.W1 <- as.numeric(sub_counts[j, 2, , ])
        N_W     <- marginTable(sub_counts[j, , ,], 1)
        N <- list(YX.W0 = N_YX.W0, YX.W1 = N_YX.W1, W = N_W)
        
        intervals <- wppIntervalGenerationNumerical(epsilons = epsilons, N = N, prior_table = prior_table / num_states, M = M, verbose = verbose)
        bounds_post[, 1] <- bounds_post[, 1] + intervals[, 1] * P_Z_hat[, j]
        bounds_post[, 2] <- bounds_post[, 2] + intervals[, 2] * P_Z_hat[, j]
        
      }    
      
    }
    
  } else {
    
    bounds_post <- NULL
    
  }
  
  # Finalize
  
  if (length(f_witness) == 0) {
    cat("No solution found\n")
  }

  out <- list(w_list = f_w_list, hw = hw, hZ = hZ, bounds = f_bounds, bounds_post = bounds_post, 
              problem = problem, epsilons = epsilons, 
              prior_ind = prior_ind, prior_table = prior_table, 
              analytical_bounds = analytical_bounds)
  class(out) <- "wpp"
  return(out)
}

wppPosteriorSampling <- function(counts, w, x, y, Z, M, prior_table = 10) {
  # Generates a sample of the posterior distribution of the contingency table parameters of a model
  # for p(W, X, Y | Z).
  #
  # * Input:
  #
  # - counts: array of binary counts of the same dimensionality of (w, x, y, Z)
  # - w, x, y, Z: the indices of the corresponding variables (Z being a vector here)
  # - M: sample size for Monte Carlo simulations
  # - prior_table: Dirichlet hyperparameter for contingency tables
  #
  # * Output:
  #
  # - theta_samples: a list with fields, 
  #
  #        sW, posterior samples of P(W = 1), 
  #      sXW0, posterior samples of P(X = 1 | W = 0, Z = z));
  #      sXW1, analogously;
  #      sYX0, posterior samples of P(Y = 1 | X = 0, Z = z);
  #      sYX1, analogously
  #
  # - Z_samples: samples from the marginal of Z, assuming its dimensionality is small enough
  #              to that there is no memory overflow
  
  return(binParamPosteriorSampling(counts, w, x, y, Z, M, prior_table))
}

wppParamPosteriorExpectation <- function(counts, w, x, y, Z, prior_table = 10) {
  # Generates the posterior expected distribution of the contingency table parameters of a model
  # for p(W, X, Y | Z).
  #
  # * Input:
  #
  # - counts: array of binary counts of the same dimensionality of (w, x, y, Z)
  # - w, x, y, Z: the indices of the corresponding variables (Z being a vector here)
  # - prior_table: Dirichlet hyperparameter for contingency tables
  #
  # * Output:
  #
  # - theta_mean: a list with fields, 
  #
  #        sW, posterior expected value of P(W = 1), 
  #      sXW0, posterior expected value of P(X = 1 | W = 0, Z = z));
  #      sXW1, analogously;
  #      sYX0, posterior expected value of P(Y = 1 | X = 0, Z = z);
  #      sYX1, analogously
  #
  # - Z_mean: posterior expected values of the marginal of Z, 
  #
  # It is assumed the dimensionality of the problem is small so that there is no memory overflow.
  
  return(binParamPosteriorExpectation(counts, w, x, y, Z, prior_table))
}

wppIntervalGenerationNumerical <- function(epsilons, theta_samples = NULL, N = NULL, prior_table = 10, M = 1000, verbose = FALSE) {
  # Using a Monte Carlo procedure, we calculate posteriors over causal bounds.
  #
  # * Input
  #
  # - theta_samples:
  # - N: if provided instead of theta_samples, this is a list containing the fields
  #      YX.W0, X.W1: two matrices of four columns each, where column indices 1, 2, 3, 4 refer to
  #                   counts of (Y = 0, X = 0), (Y = 0, X = 1), (Y = 1, X = 0), (Y = 1, X = 1)
  #      W: array with counts for W = 0 and W = 1
  # - epsilons: vector containing the relaxations for the various constraints. In particular,
  #
  #     |eta_x0^star - eta_x1^star| <= epsilons[1]
  #     |eta_x0^star - P(Y = 1 | X = x, W = 0)| <= epsilons[2]
  #     |eta_x1^star - P(Y = 1 | X = x, W = 1)| <= epsilons[3]
  #     |delta_w0^star - P(X = 1 | W = 0, U)| <= epsilons[4]
  #     |delta_w1^star - P(X = 1 | W = 1, U)| <= epsilons[4]
  #     epsilons[5] * P(U) <= P(U | W = w) <= epsilons[6] * P(U)  
  #
  # - verbose: if TRUE, print messages
  
  # Unfortunately, I don't know how to get supress the "*Error: Possibly an LP cycling occurs.  Use the Criss-Cross method"
  # messages that appear when running rcdd. "sink", "capture.output", "invisible", nothing works.
  
  if (!is.null(theta_samples)) M <- length(theta_samples$sW)
  advance_m <- !is.null(theta_samples)
  
  ACE_lp <- matrix(NA, nrow = M, ncol = 2)
  m <- 1
  
  if (verbose) cat("Solving samples\n")
  attempts <- 1
  
  while ((m <= M && !advance_m) || (attempts <= M && advance_m)) {
    
    if (verbose) cat("Solving sample", m , "\n")
        
    if (advance_m) m <- attempts
    attempts <- attempts + 1
    
    if (!is.null(theta_samples)) {
      mX1W0 <- theta_samples$sXW0[m]; mX1W1 <- theta_samples$sXW1[m]
      mY1X0 <- theta_samples$sYX0[m]; mY1X1 <- theta_samples$sYX1[m]
      P_YX.W0 <- cbind((1 - mY1X0) * (1 - mX1W0), (1 - mY1X1) * mX1W0, mY1X0 * (1 - mX1W0), mY1X1 * mX1W0)
      P_YX.W1 <- cbind((1 - mY1X0) * (1 - mX1W1), (1 - mY1X1) * mX1W1, mY1X0 * (1 - mX1W1), mY1X1 * mX1W1)
      P_W     <- theta_samples$sW[m]
    } else if (!is.null(N)) {
      # Entries in YX correspond to (Y = 0, X = 0), (Y = 0, X = 1), (Y = 1, X = 0), (Y = 1, X = 1)
      # in that order    
      P_YX.W0 <- rdirichlet(1, N$YX.W0 + prior_table / 4)
      P_YX.W1 <- rdirichlet(1, N$YX.W1 + prior_table / 4)
      P_W     <- rdirichlet(1, N$W + prior_table / 2); P_W <- P_W[2]      
    } else {
      stop("Counts or probabilities should be provided in order to sample ACE bounds")
    }
    
    # First, find the extreme points in which eta_xw^star can vary.
    
    P_Y.00 <- P_YX.W0[3] / (P_YX.W0[1] + P_YX.W0[3]); if (P_YX.W0[1] + P_YX.W0[3] == 0) P_Y.00 <- 0
    P_Y.01 <- P_YX.W1[3] / (P_YX.W1[1] + P_YX.W1[3]); if (P_YX.W1[1] + P_YX.W1[3] == 0) P_Y.01 <- 0
    eta_space_0 <- getVEtaStar(epsilons, c(P_Y.00, P_Y.01))
    if (length(eta_space_0) == 0) { # Due to rounding
      next
    }
    P_Y.10 <- P_YX.W0[4] / (P_YX.W0[2] + P_YX.W0[4]); if (P_YX.W0[2] + P_YX.W0[4] == 0) P_Y.10 <- 0
    P_Y.11 <- P_YX.W1[4] / (P_YX.W1[2] + P_YX.W1[4]); if (P_YX.W1[2] + P_YX.W1[4] == 0) P_Y.11 <- 0
    eta_space_1 <- getVEtaStar(epsilons, c(P_Y.10, P_Y.11))
    if (length(eta_space_1) == 0) {
      next
    }
    
    # Obtain the extreme points of delta_w^star
    
    P_X1.W0 <- P_YX.W0[2] + P_YX.W0[4]
    P_X1.W1 <- P_YX.W1[2] + P_YX.W1[4]
    delta_space_0 <- c(max(0, P_X1.W0 - epsilons[4]), min(1, P_X1.W0 + epsilons[4]))
    delta_space_1 <- c(max(0, P_X1.W1 - epsilons[4]), min(1, P_X1.W1 + epsilons[4]))
    
    # Complete the operation to obtain extreme points of zeta
    
    all_extreme <- buildTableParameters(eta_space_0, eta_space_1, delta_space_0, delta_space_1, FALSE)
    
    # Run linear program
    
    P_ZETA <- c(P_YX.W0, P_YX.W1)
    P_matrix <- all_extreme$T2
    lpp_params <- buildLppScdd(P_matrix, c(1 - P_W, P_W), P_ZETA, epsilons[5], epsilons[6])
    if (length(lpp_params) == 0) {
      next
    }
    c_type <- rep(0, length(lpp_params$constr_type)); c_type[which(lpp_params$constr_type == "=")] <- 1
    H <- cbind(c_type, lpp_params$b, -lpp_params$A)
    result_max_relax <- NULL; result_min_relax <- NULL

    try(result_max_relax <- lpcdd(H, lpp_params$C, minimize = FALSE), silent = TRUE)
    if (length(result_max_relax) == 0) next    
    try(result_min_relax <- lpcdd(H, lpp_params$C, minimize = TRUE), silent = TRUE)
    if (length(result_min_relax) == 0) next        
    
    if (result_max_relax$solution.type != "Optimal" || result_min_relax$solution.type != "Optimal") next
    
    ACE_lp[m, 1] <- result_min_relax$optimal.value
    ACE_lp[m, 2] <- result_max_relax$optimal.value    
    m <- m + 1
    
  }
  
  return(ACE_lp)
}

getVEtaStar <- function(epsilon, Px) {
  # Get the feasible region for the local constraints on eta_xw^\star in terms of its extreme points.
  # Basically, the extreme points corresponding to
  #
  # |eta_x0^star - eta_x1^star| <= epsilon[1]
  # |eta_x0^star - P(Y = 1 | X = x, W = 0)| <= epsilon[2]
  # |eta_x1^star - P(Y = 1 | X = x, W = 1)| <= epsilon[3]
  # 0 <= eta_xw^star <= 1
  #
  # where Px[1] = P(Y = 1 | X = x, W = 0), Px[2] = P(Y = 1 | X = x, W = 1)
  #
  # We use rcdd here out of laziness, but it can be easily done analytically.
  
  A <- rbind(c(1, -1), c(-1, 1), c(1, 0), c(-1, 0), c(0, 1), c(0, -1), c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
  b <- c(epsilon[1], epsilon[1], 
         epsilon[2] + Px[1], epsilon[2] - Px[1],
         epsilon[3] + Px[2], epsilon[3] - Px[2],
         1, 0, 1, 0)
  H <- makeH(A, b, rep(0, ncol(A)), 0)
  H <- H[-1, ]
  V <- scdd(H)
  V <- V$output
  return(V[, 3:ncol(V), drop = FALSE])
}

buildTableParameters <- function(eta_space_0, eta_space_1, delta_space_0, delta_space_1, monotonic) {  
  # Build the mapping between the eta^star/delta^star parameters into the zeta^star/eta^star space.
  # For internal use of wppIntervalGeneration.
  #
  num_base_pair_eta_0 <- nrow(eta_space_0)
  num_base_pair_eta_1 <- nrow(eta_space_1)
  base_pair_delta <- rbind(c(delta_space_0[1], delta_space_1[1]), 
                           c(delta_space_0[1], delta_space_1[2]), 
                           c(delta_space_0[2], delta_space_1[1]), 
                           c(delta_space_0[2], delta_space_1[2]));
  num_base_pair_delta <- nrow(base_pair_delta)
  
  # Build T1: the extreme points in eta^star/delta^star space
  
  T1 <- matrix(rep(0, num_base_pair_eta_0 * num_base_pair_eta_1 * num_base_pair_delta * 6), ncol = 6)
  t1_row <- 1;
  for (i in 1:num_base_pair_eta_0)
    for (j in 1:num_base_pair_eta_1)
      for (k in 1:num_base_pair_delta) {
        T1[t1_row, ] <- c(eta_space_0[i, ], eta_space_1[j, ], base_pair_delta[k, ])
        t1_row <- t1_row + 1
      }
  if (monotonic) {
    remove_rows <- which(T1[, 6] < T1[, 5])
    if (length(remove_rows) > 0) T1 <- T1[-remove_rows, ]
  }
  
  # Now, map to T2: the extreme points in zeta^star/eta^star space
  
  T2 <- matrix(rep(0, nrow(T1) * 12), ncol = 12);
  
  T2[, 1] <- (1 - T1[, 1]) * (1 - T1[, 5])
  T2[, 2] <- (1 - T1[, 3]) * T1[, 5]
  T2[, 3] <- T1[, 1] * (1 - T1[, 5])
  T2[, 4] <- T1[, 3] * T1[, 5]
  T2[, 5] <- (1 - T1[, 2]) * (1 - T1[, 6])
  T2[, 6] <- (1 - T1[, 4]) * T1[, 6]
  T2[, 7] <- T1[, 2] * (1 - T1[, 6])
  T2[, 8] <- T1[, 4] * T1[, 6]
  T2[, 9:12] <- T1[, 1:4]
  
  return(list(T1 = T1, T2 = T2))
}

buildLppScdd <- function(M, P_W, P_ZETA, beta_lower, beta_upper) {
  # Build LPP solution using scdd as follows (example):
  #
  # > T <- buildTableParameters(epsilon, c(0, 1), c(0, 1))
  # > lpp_ace <- buildLppScdd(T$T2, LPP_SCALE, model$P_ZETA)
  #
  # 'scale' should be set in a way to make all coefficients reasonably representable with
  # not many digits, or representation can explode (or maybe I don't know how
  # to properly use scdd)
  # 
  # P_W and P_ZETA are vectors which respectively represent the  distributions P(W = w) and 
  # P(Y = y, X = x | W = w) as in
  #
  # P_W = c(P(W = 0), P(W = 1))
  # P_ZETA = c(P(Y = 0, X = 0 | W = 0), P(01 | 0), P(10 | 0), P(11 | 0), P(00 | 1), P(01 | 1) etc.)
  #
  # For internal use of wppIntervalGeneration.
  
  OMEGA_POS.00 <- 9
  
  ETA_POS.00 <- 13
  ETA_POS.01 <- 14
  ETA_POS.10 <- 15
  ETA_POS.11 <- 16
  
  NUM_VARS <- 16
  NUM_RELAX_CONSTRAINTS <- 40
  
  # Basic setup 
  
  V <- makeV(M)
  pre_H <- c()
  try(pre_H <- scdd(V), silent = TRUE)
  if (length(pre_H) == 0) return(list())  
  H_matrix <- pre_H$output
  ineq_idx <- which(H_matrix[,1] == 0)
  eq_idx <- which(H_matrix[,1] == 1)
  if (length(ineq_idx) == 1) {
    A_ineq0 <- matrix(c(-H_matrix[ineq_idx, 3:ncol(H_matrix)], rep(0, 4)), nrow = 1)
  } else {
    A_ineq0 <- cbind(-H_matrix[ineq_idx, 3:ncol(H_matrix)], matrix(rep(0, length(ineq_idx) * 4), ncol = 4))  
  }
  b_ineq0 <- H_matrix[ineq_idx, 2]
  
  A_eq <- matrix(0, ncol = ncol(H_matrix) + 2, nrow = 2)
  A_eq[1, 1:4] <- 1 # Sum k_{yx.w0} == 1
  A_eq[2, 5:8] <- 1 # Sum k_{yx.w1} == 1 
  b_eq <- c(1, 1)
  
  # Now add the constraints related to the randomization relaxation
  
  A <- matrix(rep(0, NUM_RELAX_CONSTRAINTS * NUM_VARS), ncol = NUM_VARS)
  b <- rep(0, NUM_RELAX_CONSTRAINTS)
  pc <- 1
  
  # This relates kappa to the observables
  
  for (w in 0:1)
    for (x in 0:1)
      for (y in 0:1) {
        var_pos <- 1 + w * 4  + y * 2 + x * 1;
        A[pc, var_pos] <- 1
        b[pc] <-  P_ZETA[var_pos] / beta_lower
        pc <- pc + 1
        A[pc, var_pos] <- -1
        b[pc] <- -P_ZETA[var_pos] / beta_upper
        pc <- pc + 1
      }
  
  # This relates eta_xw to omega_xw
  
  for (w in 0:1)
    for (x in 0:1) {
      var_pos_eta <- ETA_POS.00 + w * 2 + x * 1;
      var_pos_omega <- OMEGA_POS.00 + w * 2 + x * 1;
      A[pc, var_pos_eta] <- 1
      A[pc, var_pos_omega] <- -beta_upper 
      pc <- pc + 1
      A[pc, var_pos_eta] <- -1
      A[pc, var_pos_omega] <- beta_lower
      pc <- pc + 1
    }
  
  # 0 <= eta_xw <= 1
  
  for (i in 0:3) {
    A[pc, ETA_POS.00 + i] <-  1; b[pc] <- 1; pc <- pc + 1
    A[pc, ETA_POS.00 + i] <- -1;             pc <- pc + 1
  }
  
  # 0 <= omega_xw <= 1
  
  for (i in 0:3) {
    A[pc, OMEGA_POS.00 + i] <-  1; b[pc] <- 1; pc <- pc + 1
    A[pc, OMEGA_POS.00 + i] <- -1;             pc <- pc + 1
  }
  
  A_ineq <- rbind(A_ineq0, A)
  b_ineq <- c(b_ineq0, b)  
  
  # Build objective function
  
  C <- rep(0, NUM_VARS)
  C[ETA_POS.00] <- -P_W[1]
  C[ETA_POS.01] <- -P_W[2]
  C[ETA_POS.10] <-  P_W[1]
  C[ETA_POS.11] <-  P_W[2]
  
  # Constraint encoding
  
  constr_type <- c(rep("<=", nrow(A_ineq)), rep("=", nrow(A_eq)))
  
  # Done, return information
  
  return(list(A = rbind(A_ineq, A_eq), b = c(b_ineq, b_eq), C = C, constr_type = constr_type))
  
}

wppIntervalGenerationAnalytical <- function(P_YX.W0, P_YX.W1, P_W, epsilons) {
  # "Message passing" implementation of constraint validation. Constraints are not
  # as tight as in wppIntervalGeneration, but it doesn't require running a linear
  # programming solver.
  #
  # This function returns TRUE if parameters P_YX.W0 and P_YX.W1 satisfied the constraints
  # of the relaxed confounded model, FALSE otherwise.
  #
  # * Input
  #
  # - P_YX.W0, P_YX.W1: two matrices of four columns each, where column indices 1, 2, 3, 4 refer to
  #                     (Y = 0, X = 0), (Y = 0, X = 1), (Y = 1, X = 0), (Y = 1, X = 1)
  # - P_W: array with the entries for the marginal of W
  # - epsilons: vector containing the relaxations for the various constraints. In particular,
  #
  #     |eta_x0^star - eta_x1^star| <= epsilons[1]
  #     |eta_x0^star - P(Y = 1 | X = x, W = 0)| <= epsilons[2]
  #     |eta_x1^star - P(Y = 1 | X = x, W = 1)| <= epsilons[3]
  #     |delta_w0^star - P(X = 1 | W = 0, U)| <= epsilons[4]
  #     |delta_w1^star - P(X = 1 | W = 1, U)| <= epsilons[4]
  #     epsilons[5] * P(U) <= P(U | W = w) <= epsilons[6] * P(U)
  #
  # Notice that the first three constraints entail
  #
  # |P(Y = 1 | X = x, W = 0) - P(Y = 1 | X = x, W = 1)| <= epsilons[1] + epsilons[2] + epsilons[3]
  #
  # which has to be encoded manually.
  #
  # * Output:
  #
  # - intervals: if return_solutions is TRUE, this is passed back - a matrix where first column
  #              are lower bounds, and second column are upper bounds
  
  N <- length(P_W)
  P_YX.W0[which(P_YX.W0 == 0)] <- 1.e-15
  P_YX.W1[which(P_YX.W1 == 0)] <- 1.e-15
  
  ## indices
  
  i_x.w   <- 1:4
  i_xp.w  <- c(2, 1, 4, 3)
  i_x.wp  <- c(3, 4, 1, 2)
  i_xp.wp <- c(4, 3, 2, 1)
  
  p_xy.w <- matrix(c(P_YX.W0, P_YX.W1), nrow = N)
  p_y.xw <- conditionMatrix(p_xy.w, 2, c(1, 3))
  p_x.w  <- marginMatrix(p_xy.w, c(1, 3))
  
  eps_w      <- epsilons[1]
  eps_Y      <- epsilons[2]
  eps_X      <- epsilons[4]
  beta_lower <- epsilons[5]
  beta_upper <- epsilons[6]
  
  # Calculate auxiliary variables
  UK_XY.W <- pmin(p_xy.w / beta_lower, 1)
  LK_XY.W <- p_xy.w / beta_upper
  
  UChi <- beta_upper * p_x.w
  LChi <- beta_lower * p_x.w
  
  L_X <- pmax(p_x.w - eps_X, 0)
  U_X <- pmin(p_x.w + eps_X, 1)
  
  U_Y <- pmin(p_y.xw[, c(3, 4, 7, 8), drop = FALSE] + eps_Y, 1)
  L_Y <- pmax(p_y.xw[, c(3, 4, 7, 8), drop = FALSE] - eps_Y, 0)
  
  U_bar <- rowMaxs(U_Y)
  L_bar <- rowMins(L_Y)
  
  ## Derive the box constraints for omega_xw first
  ## Theorem 1 bounds
  upper <-                 UK_XY.W[,c(3,4,7,8), drop = FALSE] + U_Y * UChi[,i_xp.w, drop = FALSE]
  upper <- pmin(upper,     UK_XY.W[,c(3,4,7,8), drop = FALSE] / L_X)
  upper <- pmin(upper, 1 - LK_XY.W[,c(1,2,5,6), drop = FALSE] / U_X)
  
  lower <-                 LK_XY.W[,c(3,4,7,8), drop = FALSE] + L_Y * LChi[,i_xp.w, drop = FALSE]
  lower <- pmax(lower,     LK_XY.W[,c(3,4,7,8), drop = FALSE] / U_X)
  lower <- pmax(lower, 1 - UK_XY.W[,c(1,2,5,6), drop = FALSE] / L_X)
  
  ## Theorem 2 bounds
  upper <- pmin(upper,     (UK_XY.W[,c(7,8,3,4), drop = FALSE] + eps_w * UChi[, i_x.wp, drop = FALSE]) / L_X[, i_x.wp, drop = FALSE])
  upper <- pmin(upper, 1 - (LK_XY.W[,c(5,6,1,2), drop = FALSE] - eps_w * UChi[, i_x.wp, drop = FALSE]) / U_X[, i_x.wp, drop = FALSE])
  lower <- pmax(lower,     (LK_XY.W[,c(7,8,3,4), drop = FALSE] - eps_w * UChi[, i_x.wp, drop = FALSE]) / U_X[, i_x.wp, drop = FALSE])
  lower <- pmax(lower, 1 - (UK_XY.W[,c(5,6,1,2), drop = FALSE] + eps_w * UChi[, i_x.wp, drop = FALSE]) / L_X[, i_x.wp, drop = FALSE])
  
  ## Bounds from Theorem 3
  upper <- pmin(upper,
                (UK_XY.W[,c(8,7,4,3), drop = FALSE] + UK_XY.W[,c(7,8,3,4), drop = FALSE] + UK_XY.W[,c(3,4,7,8), drop = FALSE] +
                   - LK_XY.W[,c(4,3,8,7), drop = FALSE] + UChi[,i_xp.w, drop = FALSE]*(U_bar + L_bar + 2 * eps_w) - L_bar))
  upper <- pmin(upper,
                (UK_XY.W[,c(4,3,8,7), drop = FALSE] + UK_XY.W[, c(7,8,3,4), drop = FALSE] + UK_XY.W[, c(3,4,7,8), drop = FALSE] - LK_XY.W[, c(8,7,4,3), drop = FALSE] +
                   + 2 * UChi[,i_xp.w, drop = FALSE] * eps_w + UChi[,c(4,3,2,1), drop = FALSE] * (U_bar + L_bar) - L_bar))
  lower <- pmax(lower,
                (- UK_XY.W[,c(8,7,4,3), drop = FALSE] + LK_XY.W[,c(7,8,3,4), drop = FALSE] + LK_XY.W[,c(3,4,7,8), drop = FALSE] +
                   LK_XY.W[,c(4,3,8,7), drop = FALSE] + - 2 * UChi[,i_xp.w, drop = FALSE] * eps_w + LChi[,c(4,3,2,1), drop = FALSE] * (U_bar + L_bar) - U_bar))
  lower <-  pmax(lower,
                 (-UK_XY.W[,c(4,3,8,7), drop = FALSE] + LK_XY.W[,c(7,8,3,4), drop = FALSE] + LK_XY.W[,c(3,4,7,8), drop = FALSE] +
                    + LK_XY.W[,c(8,7,4,3), drop = FALSE] - UChi[,i_xp.w, drop = FALSE] * 2 * eps_w + LChi[,i_xp.w, drop = FALSE]*(U_bar + L_bar) - U_bar))
  
  upper[is.nan(upper)] <- 1
  upper <- pmin(upper, 1)
  lower[is.nan(lower)] <- 0
  lower <- pmax(lower, 0)
  
  ## bounds for omega_xw
  omega_upper <- upper
  omega_lower <- lower
  
  ## bounds for differences omega_xw - omega_xw'
  diff_upper <- matrix(eps_w, nrow = N, ncol = 4)
  diff_lower <- matrix(-eps_w, nrow = N, ncol = 4)
  
  ## Now, iterate over linear constraints
  for (iter in 1:4) {
    
    upper <- omega_upper
    lower <- omega_lower
    
    ## #############################################
    ## Iteration over the linear constraints of Theorem 2
    upper <- pmin(upper,
                 omega_upper[,i_x.wp, drop = FALSE]*U_X[,i_xp.w, drop = FALSE] +
                   UK_XY.W[,c(3,4,7,8), drop = FALSE] + eps_w*UChi[,i_xp.w, drop = FALSE])
    upper <- pmin(upper,
                 (omega_upper[,i_x.wp, drop = FALSE] - 1)*L_X[,i_xp.w, drop = FALSE] +
                   1 - LK_XY.W[,c(1,2,5,6), drop = FALSE]  + eps_w * UChi[,i_xp.w, drop = FALSE])
    upper <- pmin(upper, omega_upper[,i_x.wp, drop = FALSE] + eps_w)
    
    lower <- pmax(lower,
                 omega_lower[,i_x.wp, drop = FALSE]*L_X[,i_xp.w, drop = FALSE] +
                   LK_XY.W[,c(3,4,7,8), drop = FALSE] - eps_w * UChi[,i_xp.w, drop = FALSE])
    lower <- pmax(lower,
                 (omega_lower[,i_x.wp, drop = FALSE] - 1)*U_X[,i_xp.w, drop = FALSE] +
                   1 - UK_XY.W[,c(1,2,5,6), drop = FALSE] - eps_w * UChi[,i_xp.w, drop = FALSE])
    lower <- pmax(lower, omega_lower[,i_x.wp, drop = FALSE] - eps_w)
    
    omega_upper <- upper
    omega_lower <- lower
    
    ## #############################################
    ## Iteration over the linear constraints of Theorem 2 to bound omega_xw - omega_xw'
    ## equation (9)
    upper <- pmin(diff_upper,
                 omega_upper[,i_x.wp, drop = FALSE]*(U_X[,i_xp.w, drop = FALSE] - 1) +
                   UK_XY.W[,c(3,4,7,8), drop = FALSE] + eps_w * UChi[,i_xp.w, drop = FALSE])
    upper <- pmin(upper,
                 (omega_upper[,i_x.wp, drop = FALSE] - 1)*(L_X[,i_xp.w, drop = FALSE] - 1)  +
                   - LK_XY.W[,c(1,2,5,6), drop = FALSE]  + eps_w * UChi[,i_xp.w, drop = FALSE])
    upper <- pmin(upper, omega_upper - omega_lower[,i_x.wp, drop = FALSE])
    
    lower <- pmax(diff_lower,
                 omega_lower[,i_x.wp, drop = FALSE]*(L_X[,i_xp.w, drop = FALSE]-1) +
                   LK_XY.W[,c(3,4,7,8), drop = FALSE] - eps_w * UChi[,i_xp.w, drop = FALSE])
    lower <- pmax(lower,
                 (omega_lower[,i_x.wp, drop = FALSE] - 1)*(U_X[,i_xp.w, drop = FALSE]-1) +
                   - UK_XY.W[,c(1,2,5,6), drop = FALSE] - eps_w * UChi[,i_xp.w, drop = FALSE])
    lower <- pmax(lower, omega_lower - omega_upper[,i_x.wp, drop = FALSE])
    
    diff_upper <- upper
    diff_lower <- lower
    
    ## #############################################
    ## Iteration over the linear constraints of Theorem 3 to bound omega_xw
    upper <- pmin(omega_upper,
                 diff_upper[,i_xp.w, drop = FALSE] - LK_XY.W[,c(4,3,8,7), drop = FALSE] + UK_XY.W[,c(3,4,7,8), drop = FALSE] +
                   + UK_XY.W[,c(8,7,4,3), drop = FALSE] + UK_XY.W[,c(7,8,3,4), drop = FALSE] +
                   - LChi[,i_x.w, drop = FALSE]*(U_bar + L_bar) + 2*eps_w + UChi[,i_x.wp, drop = FALSE] + U_bar)
    upper <- pmin(upper,
                 diff_upper[,i_xp.wp, drop = FALSE] - LK_XY.W[,c(8,7,4,3), drop = FALSE] + UK_XY.W[,c(7,8,3,4), drop = FALSE] +
                   + UK_XY.W[,c(4,3,8,7), drop = FALSE] + UK_XY.W[,c(3,4,7,8), drop = FALSE] +
                   + 2*eps_w*UChi[,i_x.wp, drop = FALSE] - LChi[,i_x.wp, drop = FALSE]*(U_bar + L_bar) + U_bar)
    
    lower <- pmax(omega_lower,
                 - diff_lower[,i_xp.wp, drop = FALSE] - UK_XY.W[,c(8,7,4,3), drop = FALSE] + LK_XY.W[,c(3,4,7,8), drop = FALSE] +
                   + LK_XY.W[,c(4,3,8,7), drop = FALSE] + LK_XY.W[,c(7,8,3,4), drop = FALSE] +
                   - UChi[,i_x.wp, drop = FALSE] * (U_bar + L_bar + 2*eps_w) + L_bar)
    lower <- pmax(lower,
                 - diff_lower[,i_xp.w, drop = FALSE] - UK_XY.W[,c(4,3,8,7), drop = FALSE] + LK_XY.W[,c(7,8,3,4), drop = FALSE] +
                   + LK_XY.W[,c(8,7,4,3), drop = FALSE] + LK_XY.W[,c(3,4,7,8), drop = FALSE] +
                   - 2*eps_w*UChi[,i_x.wp, drop = FALSE] - UChi[,i_x.w, drop = FALSE]*(U_bar + L_bar) + L_bar)
    
    omega_upper <- upper
    omega_lower <- lower
    
  }
  
  intervals <- matrix(0, nrow = N, ncol = 2)
    
  alpha_upper <- beta_upper * pmin(omega_upper, 1)
  alpha_lower <- beta_lower * omega_lower
    
  intervals[, 2] <- (alpha_upper[, 4] - alpha_lower[, 3]) * P_W + (alpha_upper[, 2] - alpha_lower[, 1]) * (1 - P_W)
  intervals[, 1] <- (alpha_lower[, 4] - alpha_upper[, 3]) * P_W + (alpha_lower[, 2] - alpha_upper[, 1]) * (1 - P_W)
  
  intervals[which(intervals[, 1] < -1)] <- -1
  intervals[which(intervals[, 2] >  1)] <-  1
  
  return(intervals)
  
}

wppSummarizeBounds <- function(object, taboo_vars = c()) {
  # Gets the output of a wpp search and provides a handful of summaries.
  #
  # * Input:
  #
  # - object: the output of a run of wpp
  # - taboo_vars: an array of integers representing variable indices. 
  #               Don't include in the summary bounds that made use of a conditioning set
  #               that includes any of the variables in this list
  #
  # * Output:
  #
  # - min_l, max_u: minimum of all lower bounds, maximum of all upper bounds
  # - narrowest_bound: narrowest of all bounds
  # - hscoring_bound: bound associated with highest score
  # - chosen_w, chosen_Z: witness and background set of the hscoring_bound
 
  if (class(object) != "wpp") {
    stop("a wpp object is necessary")
  }
    
  if (length(object$w_list$witness) == 0) {
    return(list(min_l = -1, max_u = 1, narrowest_bound = c(-1, 1), highest_bound = c(-1, 1), chosen_w = 0, chosen_Z = c()))
  }
  
  b <- object$bounds
  b[taboo_vars, ] <- NA
  
  min_l <- min(b[, 1], na.rm = TRUE)
  max_u <- max(b[, 2], na.rm = TRUE)    
  narrowest_bound <- b[which.min(b[,2] - b[,1]), ]
  
  scores <- rowSums(object$w_list$witness_score)
  idx_highest <- which.max(scores)
  hscore_bound <- b[idx_highest, ]
  chosen_w <- object$w_list$witness[idx_highest]
  chosen_Z <- object$w_list$Z[[idx_highest]]
  
  return(list(min_l = min_l, max_u = max_u, narrowest_bound = narrowest_bound, 
              chosen_w = chosen_w, chosen_Z = chosen_Z, hscore_bound = hscore_bound))
}

#' @title Summarize Witness Protection Program Outputs
#' 
#' @description 
#' \code{summary} method for class "\code{wpp}".
#'
#' @param object an object of class "\code{wpp}", usually a result of a call to \code{\link{wpp}}.
#' @param ... other parameters, ignored.
#' @return Besides fields inherented from the \code{wpp} object, a list summary statistics is included:
#'  \item{\code{lq}}{an array of 5 entries corresponding to evenly space quantiles of the lower bound of the ACE of the highest scoring
#'                  witness/admissible set pair.}
#'  \item{\code{lci}}{95\% marginal posterior credible interval for the chosen ACE lower bound.}
#'  \item{\code{uq}}{an array of 5 entries corresponding to evenly space quantiles of the lower bound of the ACE of the highest scoring
#'                  witness/admissible set pair.}
#'  \item{\code{uci}}{95\% marginal posterior credible interval for the chosen ACE lower bound.}
#'  \item{\code{hscore_bound}}{score of the witness/admissible oair of highest score.}
#'  \item{\code{chosen_w}}{witness of the witness/admissible pair of highest score.}
#'  \item{\code{chosen_Z}}{admissible set of the witness/admissible set of highest score.}
#'  \item{\code{min_l}}{estimated minimum lower bound of all lower bounds found by the procedure.}
#'  \item{\code{max_u}}{estimated maximum upper bound of all upper bounds found by the procedure.}
#'  \item{\code{narrowest_bound}}{lower and upper bounds corresponding to the narrowest ACE interval found by the procedure.}
#'  
#' @seealso
#' The model fitting function \code{\link{wpp}}.
#' 
#' @export

summary.wpp <- function(object, ...) {
  
  if (class(object) != "wpp") {
    stop("a wpp object is necessary")
  }
  
  num_bounds <- nrow(object$bounds)
  
  if (num_bounds == 0 ) {
    sprintf("No witnesses found\n")
    return
  }

  s <- wppSummarizeBounds(object)
  if (is.null(object$bound_post)) {
    lq  <- quantile(object$bounds_post[, 1], probs = seq(0, 1, 0.25), na.rm = TRUE)
    lci <- quantile(object$bounds_post[, 1], probs = c(0.025, 0.975), na.rm = TRUE)
    uq  <- quantile(object$bounds_post[, 2], probs = seq(0, 1, 0.25), na.rm = TRUE)
    uci <- quantile(object$bounds_post[, 2], probs = c(0.025, 0.975), na.rm = TRUE)
  } else {
    lq  <- NULL
    lci <- NULL
    uq  <- NULL
    uci <- NULL
  }
  
  ans <- list(lq = lq, lci = lci, uq = uq, uci = uci, 
              hscore_bound = s$hscore_bound, chosen_w = s$chosen_w, chosen_Z = s$chosen_Z,
              min_l = s$min_l, max_u = s$max_u, narrowest_bound = s$narrowest_bound,
              X_idx = object$problem$X_idx, Y_idx = object$problem$Y_idx, 
              witness = object$w_list$witness, bounds = object$bounds,
              analytical_bounds = object$analytical_bounds, varnames = object$problem$varnames)
  class(ans) <- "summary.wpp"
  return(ans)
  
}

#' @title Print Summaries of Witness Protection Program Outputs
#' 
#' @description 
#' Print output of \code{summary} method for class "\code{wpp}".
#'
#' @param x an object of class "\code{summary.wpp}", usually a result of a call to \code{\link{summary.wpp}}.
#' @param ... other parameters, ignored.
#'  
#' @details
#' Variable names with more than 20 characters are truncated when printing.
#' 
#' @export

print.summary.wpp <- function(x, ...) {
  
  if (class(x) != "summary.wpp") {
    stop("a summary.wpp object is necessary")
  }
  
  num_bounds <- nrow(x$bounds)
  
  if (num_bounds == 0 ) {
    sprintf("No witnesses found\n")
    return
  }
  
  cat("\n")
  cat(sprintf("WITNESS PROTECTION PROGRAM (Treatment %s, Outcome %s)\n",x$varnames[[x$X_idx]], x$problem$varnames[[x$Y_idx]]))
  cat(sprintf("Found %d pairs of witness/admissible sets.\n\n", num_bounds))  
  
  cat("* All bounds:\n\n")
  
  max_name <- 20
  max_seen_name <- nchar("  Witness")
  for (i in seq_len(num_bounds)) {
    max_seen_name <- max(max_seen_name, nchar(x$varnames[[x$witness[i]]]))
  }
  
  cat(sprintf(paste("%", max_seen_name, "s: Interval\n", sep = ""), "  Witness"))
  for (i in 1:num_bounds) {
    nm <- x$varnames[[x$witness[i]]]
    nm <- substr(nm, 1, min(max_name, nchar(nm)))
    cat(sprintf(paste("%", max_seen_name, "s: [% 1.2f, % 1.2f]\n", sep = ""), nm, x$bounds[i, 1], x$bounds[i, 2]))
  }
  cat("\n")
  cat(sprintf("* Smallest lower bound: % 1.2f\n", x$min_l))
  cat(sprintf("   Largest upper bound: % 1.2f\n", x$max_u))
  cat(sprintf("    Narrowest interval: [% 1.2f, % 1.2f]\n", x$narrowest_bound[1], x$narrowest_bound[2]))
  
  cat("\n")
  cat("* Highest score model:\n\n")
  cat(sprintf("         Witness: %s\n", x$varnames[[x$chosen_w]]))
  if (length(x$chosen_Z) == 0) {
    cat(sprintf("  Admissible set: EMPTY\n"))
  } else {
    cat(sprintf("  Admissible set: %s\n", x$varnames[[x$chosen_Z[1]]]))
    for (i in (1 + seq_len(length(x$chosen_Z) - 1))) {
      cat(sprintf("                  %s\n", x$varnames[[x$chosen_Z[i]]]))
    }
  }    
  cat(sprintf("        Interval: [%1.2f, %1.2f]\n", x$hscore_bound[1], x$hscore_bound[2]))
  
  if (is.null(x$bound_post)) {
    
    lq  <- x$lq
    lci <- x$lci
    uq  <- x$uq
    uci <- x$uci
    
    cat("\n")
    
    if (x$analytical_bounds) {
      cat(sprintf("          Lower bound quantiles (relaxed): (% 1.2f, % 1.2f, % 1.2f, % 1.2f, % 1.2f)\n", lq[[1]], lq[[2]], lq[[3]], lq[[4]], lq[[5]]))
      cat(sprintf("  95 per cent credible interval (relaxed): [% 1.2f, % 1.2f]\n", lci[[1]], lci[[2]]))    
      cat(sprintf("          Upper bound quantiles (relaxed): (% 1.2f, % 1.2f, % 1.2f, % 1.2f, % 1.2f)\n", uq[[1]], uq[[2]], uq[[3]], uq[[4]], uq[[5]]))
      cat(sprintf("  95 per cent credible interval (relaxed): [% 1.2f, % 1.2f]\n", uci[[1]], uci[[2]]))    
    } else {
      cat(sprintf("          Lower bound quantiles: (% 1.2f, % 1.2f, % 1.2f, % 1.2f, % 1.2f)\n", lq[[1]], lq[[2]], lq[[3]], lq[[4]], lq[[5]]))
      cat(sprintf("  95 per cent credible interval: [% 1.2f, % 1.2f]\n", lci[[1]], lci[[2]]))    
      cat(sprintf("          Upper bound quantiles: (% 1.2f, % 1.2f, % 1.2f, % 1.2f, % 1.2f)\n", uq[[1]], uq[[2]], uq[[3]], uq[[4]], uq[[5]]))
      cat(sprintf("  95 per cent credible interval: [% 1.2f, % 1.2f]\n", uci[[1]], uci[[2]]))    
      
    }
    
  }
  
}
