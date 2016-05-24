############################################################################################
# search.R
#
# Search procedures for finding causal effects. In particular, the main search algorithm is
# the one by 
#
# * D. Entner, P.O. Hoyer and P. Spirtes (2013)
#   "Data-driven covariate selection for nonparametric estimation of causal effects"
#   AISTATS 2013
#
# Code by
#
#  - Ricardo Silva (ricardo@stats.ucl.ac.uk)
#  - Robin Evans (robin.evans@stats.ox.ac.uk)
#
# Current version: 28/04/2015
# First version: 31/03/2014

#' @title Search for Causal Effect Covariate Adjustment
#' 
#' @description 
#' Find the witnesses and adjustment sets (if any) for the average causal effect (ACE) between a given treatment
#' variable \eqn{X} on a given outcome \eqn{Y}. This is done by an exhaustive search on a (reduced) set of possible 
#' candidates. Currently, only binary data is supported.
#'
#' @references 
#' \url{http://jmlr.org/proceedings/papers/v31/entner13a.html}
#' 
#' \url{http://papers.nips.cc/paper/5602-causal-inference-through-a-witness-protection-program}
#'
#' @param problem a \code{\link{cfx}} problem instance for the ACE of a given treatment \eqn{X} on a given outcome \eqn{Y}.
#' @param max_set maximum size of conditioning set. The cost of the procedure grows exponentially as a function of this, 
#'                so be careful when increasing the default value.
#' @param min_only for each witness, once a set of a particular size is found, don't look for larger ones.
#' @param prior_ind prior probability of an independence.
#' @param prior_table effective sample size hyperparameter of a Dirichlet prior for testing independence with contingency tables.
#' @param cred_calc if \code{TRUE}, compute conditional credible intervals for the ACE of highest scoring model.
#' @param M if necessary to compute (conditional) credible intervals, use Monte Carlo with this number of samples.
#' @param stop_at_first if \code{TRUE}, stop as soon as some witness is found.
#' @param pop_solve if \code{TRUE}, assume we know the population graph in \code{problem} instead of data.
#' @param verbose if \code{TRUE}, print out more detailed information while running the procedure.
#'
#' @return A list containing \code{problem} plus the following items:
#'   \item{\code{witness}}{array containing the indices of the witness variables.}
#'   \item{\code{Z}}{a list, where \code{Z[[i]]} is the \eqn{i}-th array containing the indices of the variables in 
#'                   the admissible set corresponding to witness \code{witness[i]}.}
#'   \item{\code{witness_score}}{array containing the scores of each witness/admissible set pair.}
#'   \item{\code{hw}}{witness corresponding to the highest scoring pair.}
#'   \item{\code{hZ}}{array containing admissible set corresponding to the highest scoring pair.}
#'   \item{\code{ACEs}}{array of average causal effects corresponding to each witness/admissible pair.}
#'   \item{\code{ACEs_post}}{array of samples corresponding to the posterior distribution of the ACE associated implied by
#'                           \code{hW} and \code{hZ}.}
#'                          
#' @details
#' The method assumes that the variables given in \code{problem} (other than \code{problem$X_idx} and
#' \code{problem$Y_idx}) are covariates which causally precede treatment and outcome. It then applies the 
#' faithfulness condition of Spirtes, Glymour and Scheines (2000, \emph{Causation, Prediction and Search}, MIT Press)
#' to derive an \emph{admissible set}: a set of covariates which removes all confounding between treatment and outcome
#' when adjusted for.
#' The necessary and sufficient conditions for finding an admissible set using the faithfulness assumption were
#' discussed by Enter, Hoyer and Spirtes (2013, \emph{JMLR W&CP}, vol. 31, 256--264). In order for a set to be proved
#' an admissible set, some auxiliary variable in the covariate set is necessary - we call this variable a "witness."
#' See Entner et al. for details. It is possible that no witness exists, which in this case the function returns an
#' empty solution. Multiple witness/admissible sets might exist. The criterion for finding a witness/admissible set
#' pair requires the testing of conditional independence constraints. The test is done by performing Bayesian model selection
#' with a Dirichlet prior over the contingency table of the variables in \code{problem} using the effective sample size
#' hyperparameter \code{prior_table}, and a prior probability of the independence hypothesis using the hyperparameter
#' \code{prior_ind}. 
#' 
#' For each witness/admissible set that passes this criterion, the function reports the posterior expected value
#' of the implied ACE for each pair, by first plugging-in the posterior expected value of the contingency table as
#' an estimate of the joint distribution. For a particular pair of witness/admissible set, chosen according to the
#' best fit to the conditional independencies required by the criterion of Enter et al. (see also Silva and Evans,
#' 2014, NIPS 298-306), we calculate the posterior distribution of the ACE. This posterior does not take into account
#' the uncertainty on the choice of witness/admissible set, but instead is the conditional posterior given this choice.
#' 
#' The search for a witness/admissible set is by brute-force: for each witness, evaluate all subsets of the remaining
#' covariates as candidate admissible sets. If there are too many covariates (more than \code{max_set}), only a filtered set
#' of size \code{max_set} is considered for each witness. The set is chosen by first scoring each covariate by its empirical mutual 
#' information with the witness given \code{problem$X_idx} and picking the top \code{max_set} elements, to which a brute-force
#' search is then applied.
#' 
#' @examples
#' 
#' ## Generate a synthetic problem
#' problem <- simulateWitnessModel(p = 4, q = 4, par_max = 3, M = 1000)
#' 
#' ## Idealized case: suppose we know the true distribution, 
#' ## get "exact" ACE estimands for different adjustment sets
#' sol_pop <- covsearch(problem, pop_solve = TRUE)
#' effect_pop <- synthetizeCausalEffect(problem)
#' cat(sprintf(
#'   "ACE (true) = %1.2f\nACE (adjusting for all) = %1.2f\nACE (adjusting for nothing) = %1.2f\n", 
#'    effect_pop$effect_real, effect_pop$effect_naive, effect_pop$effect_naive2))
#'
#' ## Perform inference and report results
#' covariate_hat <- covsearch(problem, cred_calc = TRUE, M = 1000)
#' summary(covariate_hat)
#' 
#' @export                          

covsearch <- function(problem, max_set = 12, min_only = TRUE, prior_ind = 0.5, prior_table = 10, 
                            cred_calc = FALSE, M = 1000, stop_at_first = FALSE, pop_solve = FALSE, verbose = FALSE) {
  
  if (class(problem) != "cfx") {
    stop("a CausalFX object is necessary")
  }
  
  if (!validateData(problem$data)) {
    stop("Data must be binary, encoded numerically as 0s and 1s")
  }
    
  x <- problem$X_idx
  y <- problem$Y_idx
  latents <- problem$latent_idx
  if (pop_solve) {
    num_v <- ncol(problem$graph)
    dep_params <- list(pop_solve = pop_solve, graph = problem$graph, 
                       ancestrals = problem$ancestrals)
  } else {
    num_v <- ncol(problem$data)
    dep_params <- list(pop_solve = pop_solve, prior_ind = prior_ind, prior_table = prior_table)
  }
  
  witness_pool <- seq_len(num_v)[-c(x, y, latents)]
  witness_choice <- c()
  Z_choice <- list()
  witness_score <- matrix(nrow = 0, ncol = 2)
  
  for (w in witness_pool) {
    if (length(witness_choice) > 0 && stop_at_first) break
    
    if (verbose) cat("################################# Testing witness", w, "\n")
    Z_set <- filterAdjustmentSet(problem, max_set, w, setdiff(witness_pool, w))
    num_combo <- 2^length(Z_set)
    found_size <- Inf
    
    Zs <- combinations(rep(2, length(Z_set)))
    if (!pop_solve) {
      problem$sub_counts <- dtoc(problem$data[, c(problem$X_idx, problem$Y_idx)])
    }
    
    for (i in seq_len(num_combo)) {
      if (length(witness_choice) > 0 && stop_at_first) break
      Z <- Z_set[Zs[i,] > 0]
      if (length(Z) > found_size) next

      # Test w independent of y given Z
      if (pop_solve) {
        dep_params$x <- w; dep_params$y <- y; dep_params$Z <- Z
      } else {
        sub_table <- dtoc(problem$data[, c(w, y, Z, x)])
        dep_params$sub_counts <- marginTable(sub_table, 1:(2 + length(Z)))
        dep_params$x <- 1; dep_params$y <- 2; dep_params$Z <- 2 + seq_along(Z)
      }
      D1 <- inferDep(dep_params)
      # Test w independent of y given (Z, x)
      if (pop_solve) {
        dep_params$x <- w; dep_params$y <- y; dep_params$Z <- c(Z, x)
      } else {
        dep_params$sub_counts <- sub_table
        dep_params$Z <- 2 + seq_along(c(Z, x))
      }
      D2 <- inferDep(dep_params)
      
      # Decide on witness based on evidence above
      if (!D1$decision && D2$decision) {
        witness_choice <- c(witness_choice, w)
        Z_choice[[length(witness_choice)]] <- Z
        witness_score <- rbind(witness_score, c(D1$scores[2] - D1$scores[1], D2$scores[1] - D2$scores[2]))
        if (verbose) cat("Accepting witness", w, " set", Z, "\n")
        found_size <- length(Z)
      }
    }
  }
  
  if (length(witness_choice) == 0) {
    out <- list(witness = witness_choice, Z = Z_choice, problem = problem)
    class(out) <- "covsearch"
    return(out)
  }
  
  if (!pop_solve) {
    ACEs <- rep(0, length(witness_choice))
    for (i in seq_along(witness_choice)) {
      ACEs[i] <- bindagCausalEffectBackdoor(problem, prior_table, Z_choice[[i]])
    }
  } else {
    ACEs <- NULL
  }

  idx_highest <- which.max(rowSums(witness_score))
  chosen_w <- witness_choice[idx_highest]
  chosen_Z <- Z_choice[[idx_highest]]
  
  if (cred_calc) {
        
    sub_counts <- dtoc(problem$data[, c(chosen_w, problem$X_idx, problem$Y_idx, chosen_Z)])
    param_samples <- binParamPosteriorSampling(sub_counts, 1, 2, 3, 3 + seq_along(chosen_Z), M = M)
    theta_samples <- param_samples$theta_samples
    P_Z_hat <- param_samples$Z_samples
    
    ACE_post <- rep(0, M)
    
    num_states <- 2^length(chosen_Z)         
    for (j in seq_len(num_states)) {
      ACE_post <- ACE_post + (theta_samples[[j]]$sYX1 - theta_samples[[j]]$sYX0) * P_Z_hat[, j]
    }
    
  } else {
    
    ACE_post <- NULL
    
  }
  
  out <- list(witness = witness_choice, Z = Z_choice, witness_score = witness_score, 
              hw = chosen_w, hZ = chosen_Z, ACEs = ACEs, ACE_post = ACE_post, problem = problem)
  class(out) <- "covsearch"
  
  return(out)
}

filterAdjustmentSet <- function(problem, max_set, w, Z) {
  # In case the set of possible candidates for adjustment is too large, we pre-select those
  # variables such that the estimated association between Z[i] and w given treatment X_idx is in the
  # top max_set candidates.
  #
  # The association is given by minus the KL divergence E_Z[KL(P(X, W, | Z) || P(X | Z)P(Y | Z)] as
  # given by the empirical distribution.
  #
  # * Input:
  #
  # - problem: a WPP problem instance for the ACE of some X on some Y
  # - max_set: maximum size of conditioning set
  # - w: witness candidate
  # - Z: set of candidate adjustment variables
  #
  # * Output:
  #
  # - Z_top: subset of Z corresponding to the top scoring candidates
  
  if (length(Z) <= max_set) {
    return(Z)
  }
  
  N <- nrow(problem$data)
  p_x1 <- sum(problem$data[, problem$X_idx] == 0) / N; p_x <- c(p_x1, 1 - p_x1)  
  scores <- rep(0, length(Z))
  for (i in seq_len(length(Z))) {
    probs <- dtoc(problem$data[, c(problem$X_idx, w, Z[i])]) / N    
    for (x in 1:2) {
      for (w in 1:2) {
        for (z in 1:2) {
          if (probs[x, w, z] > 0) {
            p_w_x <- (probs[x, w, 1] + probs[x, w, 2]) / p_x[x]
            p_z_x <- (probs[x, 1, z] + probs[x, 2, z]) / p_x[x]
            scores[i] <- scores[i] - probs[x, w, z] * (log(p_w_x) + log(p_z_x))
          }
        }
      }
    }
  }
  
  s <- sort(scores, decreasing = TRUE, index.return = TRUE)
  Z_top <- sort(Z[s$ix[seq_len(max_set)]])
  return(Z_top)
}

inferDep <- function(dep_params) {
  # Chech whether x and y are independent given Z, as estimated from data. 
  #
  # * Input:
  #
  # - x, y: the indices of the two variables to test
  # - Z: the indices of the conditioning set
  # - problem: a problem instance containing the counts
  # - dep_params: a list with several fields
  #       pop_solve: if TRUE, uses the problem instance's graph instead of the counts. So this
  #                  becomes an oracle
  #       graph, ancestrals: if pop_solve is TRUE, this is the graph and the ancestral sets
  #                          of the true model
  #       x, y, Z: these are the indices of x, y and Z in graph/data
  #       sub_counts: if pop_solve if FALSE, this is the contigency table for 
  #                   {x, y, Z}, in that order
  #
  # * Output:
  #
  # - TRUE, if x is judged to be independent of y given Z, FALSE otherwise
  
  if (dep_params$pop_solve) {
    decision <- dsep(dep_params$x, dep_params$y, dep_params$Z, dep_params$graph, dep_params$ancestrals)
    if (decision) {
      scores <- log(c(1, 0))
    } else {
      scores <- log(c(0, 1))
    }
    return(list(decision = decision, scores = scores))
  }
  return(inferDepBayes(dep_params$x, dep_params$y, dep_params$Z, dep_params$sub_counts, dep_params$prior_ind, dep_params$prior_table))
}

inferDepBayes <- function(x, y, Z, counts, prior_ind = 0.5, prior_table = 10) {
  # Checks whether x and y are independent given Z, as estimated from binary data. This uses
  # the saturated model with Heckerman's BDe score with an uniform prior over the two
  # hypotheses.
  #
  # * Input:
  #
  # - x, y: the indices of the two variables to test
  # - Z: the indices of the conditioning set
  # - counts: array of counts
  # - prior_ind: prior probability of independence
  # - prior_table: hyperparameter for the contingency table prior (effective sample size)
  #
  # * Output:
  #
  # - decision: TRUE, if x is judged to be independent of y given Z, FALSE otherwise
  # - scores: scores[1] is the log marginal likelihood of the model where x and y are independent given Z,
  #           while scores[2] is the case where x and y are dependent
  
  scores <- c(0, 0)
  counts2 <- marginTable(counts, c(y, x, Z))
    
  ## get score with x as a parent
  
  alpha     <- prior_table / 2^(length(Z) + 2)
  alpha_yj0 <- alpha
  alpha_yj1 <- alpha
  alpha_yj  <- alpha_yj0 + alpha_yj1

  N_yj1 <- subtable(counts2, 1, 2)
  N_yj0 <- subtable(counts2, 1, 1)    
  N_yj  <- N_yj0 + N_yj1
  
  scores[2] <- log(1 - prior_ind) + 
               sum(lgamma(alpha_yj) - lgamma(alpha_yj + N_yj) +
                   lgamma(alpha_yj0 + N_yj0) - lgamma(alpha_yj0) +
                   lgamma(alpha_yj1 + N_yj1) - lgamma(alpha_yj1))
  
  ## now without x as a parent
  
  alpha     <- prior_table / 2^(length(Z) + 1)
  alpha_yj0 <- alpha
  alpha_yj1 <- alpha
  alpha_yj  <- alpha_yj0 + alpha_yj1  
  
  N_yj1 <- marginTable(N_yj1, c(seq_along(Z) + 1))
  N_yj0 <- marginTable(N_yj0, c(seq_along(Z) + 1))
  N_yj  <- N_yj0 + N_yj1
  
  scores[1] <- log(prior_ind) +
               sum(lgamma(alpha_yj) - lgamma(alpha_yj + N_yj) +
                   lgamma(alpha_yj0 + N_yj0) - lgamma(alpha_yj0) +
                   lgamma(alpha_yj1 + N_yj1) - lgamma(alpha_yj1))
  
  return(list(decision = scores[1] > scores[2], scores = scores))
}

covariateSearchSummarizeACEs <- function(object, taboo_vars = c()) {
  # Gets the output of a wpp search and provides a handful of summaries.
  #
  # * Input:
  #
  # - object: the output of a run of covsearch
  # - taboo_vars: an array of integers representing variable indices. 
  #               Don't include in the summary ACEs that made use of a conditioning set
  #               that includes any of the variables in this list
  #
  # * Output:
  #
  # - chosen_w, chosen_Z: witness and background set of the highest scoring pair
  
  if (class(object) != "covsearch") {
    stop("a covsearch object is necessary")
  }
  
  if (length(object$witness) == 0) {
    return(list(chosen_w = 0, chosen_Z = c()))
  }
  
  scores <- rowSums(object$witness_score)
  scores[taboo_vars] <- -Inf
  idx_highest <- which.max(scores)
  chosen_w <- object$witness[idx_highest]
  chosen_Z <- object$Z[[idx_highest]]
  chosen_ACE <- object$ACEs[[idx_highest]]
  
  return(list(chosen_w = chosen_w, chosen_Z = chosen_Z, chosen_ACE = chosen_ACE))
}

#' @title Summarize Covariate Search Outputs
#' 
#' @description 
#' \code{summary} method for class "\code{covsearch}".
#'
#' @param object an object of class "\code{covsearch}", usually a result of a call to \code{\link{covsearch}}.
#' @param ... other parameters, ignored.
#' @return Besides fields inherented from the \code{covsearch} object, a list summary statistics is included:
#'  \item{\code{q}}{an array of 5 entries corresponding to evenly space quantiles of the ACE of the highest scoring
#'                  witness/admissible set pair.}
#'  \item{\code{ci}}{95\% marginal posterior credible interval for the chosen ACE.}
#'         
#' @seealso
#' The model fitting function \code{\link{covsearch}}.
#' 
#' @export

summary.covsearch <- function(object, ...) {

  if (class(object) != "covsearch") {
    stop("a covsearch object is necessary")
  }
  
  num_w <- length(object$witness)
  
  if (num_w == 0 ) {
    sprintf("No witnesses found\n")
    return
  }
  
  s <- covariateSearchSummarizeACEs(object)
  if (!is.null(object$ACE_post)) {    
    q <- quantile(object$ACE_post, probs = seq(0, 1, 0.25))
    ci <- quantile(object$ACE_post, probs = c(0.025, 0.975))
  } else {
    q <- NULL
    ci <- NULL
  }
  
  ans <- list(q = q, ci = ci, chosen_w = s$chosen_w, chosen_Z = s$chosen_Z, chosen_ACE = s$chosen_ACE,
              ACEs = object$ACEs, 
              varnames = object$problem$varnames, witness = object$witness,
              treatment = object$problem$varnames[[object$problem$X_idx]], outcome = object$problem$varnames[[object$problem$Y_idx]])
  class(ans) <- "summary.covsearch"
  return(ans)
  
}

#' @title Print Summaries of Covariate Search Outputs
#' 
#' @description 
#' Print output of \code{summary} method for class "\code{covsearch}".
#'
#' @param x an object of class "\code{summary.covsearch}", usually a result of a call to \code{\link{summary.covsearch}}.
#' @param ... other parameters, ignored.
#' 
#' @details
#' Variable names with more than 20 characters are truncated when printing.
#' 
#' @export

print.summary.covsearch <- function(x, ...) {
  
  if (class(x) != "summary.covsearch") {
    stop("a summary.covsearch object is necessary")
  }
  
  num_w <- length(x$witness)
  
  if (num_w == 0 ) {
    sprintf("No witnesses found\n")
    return
  }
  
  cat("\n")
  cat(sprintf("COVARIATE SEARCH BY INDEPENDENCE CONSTRAINTS (Treatment %s, Outcome %s)\n",x$treatment, x$outcome))
  cat(sprintf("Found %d pairs of witness/admissible sets.\n\n", num_w))  
  
  cat("* All ACEs:\n\n")
  
  max_name <- 20
  max_seen_name <- nchar("  Witness")
  for (i in seq_len(num_w)) {
    max_seen_name <- max(max_seen_name, nchar(x$varnames[x$witness[i]]))
  }
  
  cat(sprintf(paste("  %", max_seen_name, "s: ACE\n", sep = ""), "  Witness"))
  for (i in 1:num_w) {
    nm <- x$varnames[[x$witness[i]]]
    nm <- substr(nm, 1, min(max_name, nchar(nm)))
    cat(sprintf(paste("  %", max_seen_name, "s: % 1.2f\n", sep = ""), nm, x$ACEs[i]))
  }
  
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
    cat(sprintf("\n          Estimated ACE: % 1.2f\n", x$chosen_ACE))
  }    
  if (!is.null(x$q)) {    
    q <- x$q
    ci <- x$ci
    cat(sprintf("\n  Implied ACE quantiles: (% 1.2f, % 1.2f, % 1.2f, % 1.2f, % 1.2f)\n", q[[1]], q[[2]], q[[3]], q[[4]], q[[5]]))
    cat(sprintf("\n  95 per cent credible ACE interval: [% 1.2f, % 1.2f]\n", ci[[1]], ci[[2]]))
  }
  
}
