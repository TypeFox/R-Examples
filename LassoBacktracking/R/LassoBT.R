#' Fit linear models with interactions using the Lasso.
#' 
#' Computes a number of Lasso solution paths with increasing numbers of
#' interactions present in the design matrices corresponding to each path.
#' Previous paths are used to speed up computation of subsequent paths so the
#' process is very fast.
#' 
#' @param x input matrix of dimension nobs by nvars; each row is an observation
#'   vector.
#' @param y response variable; shoud be a numeric vector.
#' @param nlambda the number of lambda values. Must be at least 3.
#' @param iter_max the number of iterations of the Backtracking algorithm to
#'   run. \code{iter_max=1} corresponds to a single lasso or elasticnet fit.
#'   Values greater than 1 will fit interactions.
#' @param lambda.min.ratio smallest value in \code{lambda} as a fraction of the
#'   largest value at which all main effects coefficients are 0.
#' @param lambda user supplied \code{lambda} sequence of decreasing penalty
#'   parameters. Typical usage is to allow the function to compute its own
#'   \code{lambda} sequence. Inappropraite sequences may cause convergence
#'   problems.
#' @param thresh convergence threshold for coordinate descent. Each inner
#'   coordinate descent loop continues until either the maximum change in the
#'   objective after any coefficient update is less than \code{thresh} or
#'   \code{1E5} iterations have been performed.
#' @param verbose if \code{TRUE} will print iteration numbers.
#' @param inter_orig an optional 2-row matrix with each column giving
#'   interactions that are to be added to the design matrix before the algorithm
#'   begins.
#' @details
#'   The Lasso optimisations are performed using coordinate descent similarly to
#'   the \pkg{glmnet} package. An intercept term is always included. Variables
#'   are centred and scaled to have equal empirical variance. Interactions are
#'   constructed from these centred and scaled variables, and the interactions
#'   themselves are also centred and scaled. Note the coefficients are returned
#'   on the original scale of the variables. Coefficients returned for
#'   interactions are for simple pointwise products of the original variables
#'   with no scaling.
#' @return An object with S3 class "\code{BT}".
#'  \describe{ 
#'   \item{\code{call}}{the call that produced the object}
#'   \item{\code{a0}}{list of intercept vectors}
#'   \item{\code{beta}}{list of matrices of coefficients
#'   stored in sparse column format (\code{CsparseMatrix})} 
#'   \item{\code{fitted}}{list of fitted values}
#'   \item{\code{lambda}}{the sequence of \code{lambda} values used}
#'   \item{\code{nobs}}{the number of observations}
#'   \item{\code{nvars}}{the number of variables} 
#'   \item{\code{var_indices}}{the indices of the non-constant columns of the
#'   design matrix}
#'   \item{\code{interactions}}{a 2-row matrix with columns
#'   giving the interactions that were added to the design matrix} 
#'   \item{\code{path_lookup}}{a matrix with columns corresponding to iterations
#'   and rows to lambda values. Entry \eqn{ij} gives the component of the
#'   \code{a0} and \code{beta} lists that gives the coefficients for the
#'   \eqn{i}th \code{lambda} value and \eqn{j}th iteration} 
#'   \item{\code{l_start}}{a vector with component entries giving the minimimum
#'   \code{lambda} index in the corresponding copmonents of \code{beta} and
#'   \code{a0}}
#'  }
#' @references
#'   Shah, R. D. (2016) \emph{Modelling interactions in high-dimensional data with
#'   Backtracking. JMLR, to appear.}
#'   \url{http://www.statslab.cam.ac.uk/~rds37/papers/shah16.pdf}
#' @seealso
#'   \code{\link{predict.BT}}, \code{\link{coef.BT}} methods and the \code{\link{cvLassoBT}}
#'   function.
#' @examples
#' x <- matrix(rnorm(100*250), 100, 250)
#' y <- x[, 1] + x[, 2] - x[, 1]*x[, 2] + x[, 3] + rnorm(100)
#' out <- LassoBT(x, y, iter_max=10)
#' @export
#' @import Matrix
#' @useDynLib LassoBacktracking
#' @importFrom Rcpp sourceCpp
LassoBT <- function(x, y, nlambda=100L, iter_max=1L, lambda.min.ratio = ifelse(nobs < nvars, 0.01, 0.0001),
                    lambda=NULL, thresh=1e-07, verbose=FALSE, inter_orig) {
  # Arguments to add later
  standardize=TRUE; intercept=TRUE; alpha=1;

  maxit <- 1e5
  if (!is.matrix(x)) stop("x should be a matrix with two or more columns")
  np <- dim(x)
  if (is.null(np) | (np[2] < 1L)) 
    stop("x should be a matrix with at least one non-constant column")
  nobs <- n <- as.integer(np[1])
  p <- as.integer(np[2])
  
  inter_space <- as.integer(iter_max * (iter_max-1) / 2) 
  # Check inter_orig
  if (missing(inter_orig)) {
    inter_orig <- matrix(0L, nrow=2, ncol=0)
    n_inter_orig <- 0L
    interactions <- matrix(0L, nrow=2, ncol=inter_space)
  } else {
    if ((!is.matrix(inter_orig)) || (nrow(inter_orig)!=2) || any(is.na(inter_orig))) {
      stop("inter_orig should be an integer matrix with two rows and all entries must be between
         1 and ncol(x)")
    }
    n_inter_orig <- ncol(inter_orig)
    inter_orig <- as.integer(inter_orig)
    dim(inter_orig) <- c(2L, n_inter_orig)
    if ((max(inter_orig) > p) || (min(inter_orig) < 1L) )
      stop("inter_orig should be an integer matrix with two rows and all entries must be between
       1 and ncol(x)")
    duplicated_cols <- duplicated(t(inter_orig))
    inter_orig <- inter_orig[, !duplicated_cols, drop=FALSE]
    rm(duplicated_cols)
    # sort inter_orig and remove pure quadratic effects
    inter_orig <- apply(inter_orig, 2, sort)
    no_pure_quad_effs <- (inter_orig[2, ] - inter_orig[1, ]) > 0
    inter_orig <- inter_orig[, no_pure_quad_effs, drop=FALSE]
    inter_orig <- inter_orig[, order(inter_orig[1, ], inter_orig[2, ]), drop=FALSE]
    interactions <- cbind(inter_orig, matrix(0L, nrow=2, ncol=inter_space))
    inter_orig <- inter_orig - 1L # so it is zero-indexed
  }
    
  inter_max <- inter_space + n_inter_orig
  p_max <- p + inter_max

  this.call <- match.call()
  
  # Check y ##
  y <- as.numeric(y)
  if (length(y) != n) stop("y does not have the same number of observations as x")
  if (any(is.na(y)) || any(is.na(x))) stop("No missing values in x or y allowed") 
  
  if (alpha > 1) {
    warning("alpha > 1; set to 1")
    alpha <- 1
  }
  if (alpha < 0) {
    warning("alpha < 0; set to 0")
    alpha <- 0
  }
  alpha <- as.double(alpha)
  
  # Intercept
  if (intercept == TRUE) {
    mean_y <- mean(y)
    y_c <- y - mean_y
    null_dev <- 2*sum(y_c^2) # should change to mean later
    resid_cur <- matrix(y_c, ncol=nlambda, nrow=length(y_c)) # deep copy
  }
  
  # Standardize x
  centres <- rep(0.0, p_max)
  scales <- rep(0.0, p_max)
  var_names_main <- integer(p)
  # enlarged x matrix to which we will add interactions
  X <- scale_cen(x, scales, centres, p, p_max, var_names_main) # scales, centres and p can change
  if (p < 1L)
    stop("x should be a matrix with at least one non-constant column")
  length(var_names_main) <- nvars <- p_eff_old <- p_eff <- p
  p_max <- p + inter_max
  length(centres) <- length(scales) <- p_max
  length(X) <- n*p_max
  dim(X) <- c(n, p_max)
  
  # Potentially add interactions to x
  p_eff_old <- p_eff <- add_inter_orig(X, scales, centres, p_eff, inter_orig)
  
  # active_sets
  strong_act_comp_set <- c(rep(TRUE, p_eff), rep(FALSE, p_max - p_eff))
  # complement of the union of strong and active sets except we set the interaction variables to FALSE
  strong_setdiff <- integer(0) # strong set minus active set
  active_set <- integer(0) # will always contain all non-zeroes
  violations <- rep(FALSE, p_max) # initially no violations
  active_used <- vector(mode="list", length=nlambda) # list of the acive sets used in the previous iteration
  for (j in seq_len(nlambda)) active_used[[j]] <- integer(0)
  rm(j)
    
  # objects for determining what to watch
  interacting_mains <- integer(0)
  l_watch <- 0L # the point where active_mains had elements outside interacting mains, and
  # so a new watch set was put in place
  terminate <- FALSE

  # Create lambda sequence
  inner_prod_abs <- rep(0, p_max)
  if (is.null(lambda)) {
    nlambda <- as.integer(nlambda)
    # need to check nlambda ##
    if (nlambda < 3L) {
      warning("Setting nlambda to 3")
      nlambda <- 3L
    }
    inner_prod_abs[1:p_eff] <- colSums(X[, 1:p_eff] * y_c)
    # inner_prod_abs <- inner_prod_abs_comp2(X, strong_act_comp_set, resid_cur, p, 0L)
    active_set <- which.max(inner_prod_abs[1:p_eff])
    # active_set contains the non-zeroes but may contain extra vars
    strong_act_comp_set[active_set] <- FALSE
    lambda_max <- inner_prod_abs[active_set]
    lambda_min <- lambda_max * lambda.min.ratio
    lambda <- exp(seq(from=log(lambda_max), to=log(lambda_min), length.out=nlambda))
  } else {
    # check lambda
    nlambda <- length(lambda)
    if (!all(diff(lambda) < 0)) stop("lambda must be a decreasing sequence")
    if (nlambda < 3L) stop("lambda must have at least 3 components")
    # inner_prod_abs <- rep(0, p_max)
  }
  alpha_lam_vec <- alpha * lambda
  alpha_lam_div_vec <- 1 + lambda*(1-alpha)
  
  
  # pat_lookup, iterations etc.
  path_lookup <- matrix(NA, nrow=nlambda, ncol=iter_max)
  path_lookup[, 1L] <- 1L
  # path_lookup[i, j] gives the the component of beta_all that will contain the coef for l=j, iter=i
  # first row will be 1's
  
  l <- 2L
  l0 <- 1L
  l_start <- 1L
  l_start_vec <- integer(iter_max)
  n_interactions <- integer(iter_max)
  n_interactions[1L] <- p_eff # we subtract p from this later
  l_start_vec[1L] <- 1L
  iter <- 1L
  cur_path <- 1L
  
  # coefficients
  beta_all <- vector(iter_max, mode="list")
  beta_cur <- matrix(0, nrow=p_max, ncol=nlambda)
  beta_all[[1]] <- matrix(0, nrow=p_max, ncol=nlambda)
  a0 <- vector(iter_max, mode="list")
  a0[[1]] <- numeric(nlambda)
  
  repeat {
    repeat {
      lambda_cur <- lambda[l]
      repeat {
        beta_active(X, beta_cur, resid_cur, active_set, l0, thresh, maxit, lambda[l], alpha_lam_vec[l],
                    alpha, alpha_lam_div_vec[l], null_dev) # will modify resid_cur and beta
        inner_prod_abs[strong_setdiff] <- inner_prod_abs_comp(X, strong_setdiff, resid_cur, l0)
        violations[strong_setdiff] <- inner_prod_abs[strong_setdiff] > lambda_cur
        if (any(violations[strong_setdiff])) {
          # strong predictors violated
          # update active_set and strong_setdiff
          strong_set_violators <- which(violations[strong_setdiff])
          # strong_set_violators indexes strong_setdiff
          active_set <- c(active_set, strong_setdiff[strong_set_violators])
          violations[strong_setdiff[strong_set_violators]] <- FALSE
          
          strong_setdiff <- strong_setdiff[-strong_set_violators]
        } else {
          # no strong predictors violating
          # check other predictors
          inner_prod_abs[strong_act_comp_set] <- inner_prod_abs_comp2(X, strong_act_comp_set, resid_cur, p_eff, l0)
          violations[strong_act_comp_set] <- inner_prod_abs[strong_act_comp_set] > lambda_cur
          if (any_indmax(violations, p_eff)) {
            # Weak predictors violated
            w_violations <- which_indmax(violations, p_eff)
            active_set <- c(active_set, w_violations)
            ##
            strong_act_comp_set[w_violations] <- FALSE
            violations[w_violations] <- FALSE # note equivalent of violations <- rep(FALSE, p)
          } else {
            # beta[, l] is correct. Now update strong set
            if (l < nlambda) {
              strong_setdiff_add1 <- in_log(active_used[[l+1L]], strong_act_comp_set)
              # used active set that is in strong_act_comp_set
              strong_act_comp_set[strong_setdiff_add1] <- FALSE
              
              strong_setdiff_add2 <- which(strong_act_comp_set & (inner_prod_abs > 2*lambda[l+1L] - lambda[l]))
              strong_act_comp_set[strong_setdiff_add2] <- FALSE
              # note 1 and 2 are disjoint
              
              strong_setdiff <- c(strong_setdiff, strong_setdiff_add1, strong_setdiff_add2)
            }
            break
          }
        }
      }

      
      if (l >= nlambda) break
      l <- l + 1L
      l0 <- l0 + 1L
      beta_cur[active_used[[l]], l] <- 0 # set beta_cur[, l] to zero
      if (l >= 3L) {
        # l >= 3 unless one of l_start = 1
        beta_cur[active_set, l] <- beta_cur[active_set, l-1L] + ((lambda[l-1L]-lambda[l])/(lambda[l-2L]-lambda[l-1L]))*
          (beta_cur[active_set, l-1L] - beta_cur[active_set, l-2L])
      } else if (l == 2L) {
        # l >= 2
        beta_cur[active_set, l] <- beta_cur[active_set, l-1L]
      }
      resid_cur[, l] <- as.numeric(y_c - X[, active_set, drop=FALSE] %*% beta_cur[active_set, l, drop=FALSE])
      # Make use of previous active sets
      active_used[[l]] <- active_set
    }
    
    # copy beta_cur to beta_all
    beta_all[[cur_path]] <- Matrix(beta_cur[, l_start:nlambda])
    
    if (iter >= iter_max) {
      terminate <- TRUE
      break
    }

    # Backtracking part: find next l_start
    repeat {
      if (l_watch >= nlambda) {
        terminate <- TRUE
        break
      }
      l_watch <- l_watch + 1L
      
      # True non-zeroes
      #active_true <- active_used[[l_watch]][which(beta_all[[path_lookup[l_watch, iter]]][active_used[[l_watch]], l_watch]!= 0)]
      active_true <- active_used[[l_watch]][which(beta_cur[active_used[[l_watch]], l_watch] != 0)]
      active_mains_ind <- active_true <= p
      active_mains <- active_true[active_mains_ind]
      
      if (length(active_mains) > length(interacting_mains)) {
        active_main_new <- active_mains %w/o% interacting_mains
        # check size
        extra_inter <- as.integer(length(active_main_new)*(length(interacting_mains)+
                                                             (length(active_main_new)-1)/2))
        if ((inter_space >= extra_inter) && (extra_inter > 0L)) {
          iter <- iter+1L
          if (verbose) cat(paste("Iteration", iter, "\n"))
          
          p_eff_old <- p_eff
          # add interactions to X (since we know there is space)
          p_eff <- add_inter(X, interacting_mains, active_main_new, scales, centres, p_eff, interactions, p, inter_orig)
          
          inter_space <- p_max-p_eff
          interacting_mains <- c(interacting_mains, active_main_new)
          # find smallest l where KKT conditions are violated or return l+1 if no violations
          l0_start <- find_l0(X, p_eff_old, p_eff, resid_cur, nlambda-1L, lambda)
          l_start <- l0_start+1L # this already increments l
          if (l0_start > 0L) {
            path_lookup[1L:l0_start, iter] <- path_lookup[1L:l0_start, iter-1L]
          } else {
            path_lookup[, iter] <- cur_path + 1L
          }
          if (l0_start < nlambda) {
            l_watch <- min(l0_start, l_watch) # i.e. l_watch <- l-1 after l has been updated
            cur_path <- cur_path + 1L
            
            # reset active sets etc.
            strong_act_comp_set[(p_eff_old+1L):p_eff] <- TRUE
            strong_act_comp_set[active_set] <- TRUE
            active_set <- active_used[[l_start]]
            strong_act_comp_set[active_set] <- FALSE
            strong_setdiff <- integer(0)
            
            n_interactions[cur_path] <- p_eff
            l_start_vec[cur_path] <- l_start
            path_lookup[l_start:nlambda, iter] <- cur_path # changed l0_start to l_start
            break
          }
        } else if (inter_space < extra_inter) {
          terminate <- TRUE
          break
        }
      }
    }
    if (terminate) break
    l <- l_start
    l0 <- l0_start
  }
    
  if (verbose) cat("Finished all iterations\nNow rescaling coefficients and computing fitted values...\n")
  # shorten objects
  length(beta_all) <- cur_path
  length(a0) <- cur_path
  length(l_start_vec) <- cur_path
  length(n_interactions) <- cur_path
  n_interactions <- n_interactions-p
  length(path_lookup) <- nlambda*iter
  dim(path_lookup) <- c(nlambda, iter)
  length(scales) <- p_eff
  length(centres) <- p_eff
  length(interactions) <- 2L*(p_eff-p)
  length(X) <- n*p_eff
  dim(X) <- c(n, p_eff)

  fit_all <- vector(mode="list", length=cur_path)
  beta_add <- matrix(0.0, nrow=p, ncol=nlambda)
  
  # rename interactions
  interactions <- var_names_main[interactions]
  dim(interactions) <- c(2L, p_eff-p)
  rownames(interactions) <- c("var1", "var2")
  
  var_names <- c(as.character(var_names_main), apply(interactions, 2, function(x) paste(x, collapse=":")))
  
  for (k in seq_along(beta_all)) {
    #beta_all[[k]] <- Matrix(beta_all[[k]][, l_start_vec[k]:nlambda])
    colnames(beta_all[[k]]) <- l_start_vec[k]:nlambda
    change_dim(beta_all[[k]]@Dim, p_eff)
    rownames(beta_all[[k]]) <- var_names
    fit_all[[k]] <- as.matrix(X %*% beta_all[[k]]) + mean_y
    # rescale beta
    beta_all[[k]] <- beta_all[[k]] / scales
    a0[[k]] <- mean_y - colSums(beta_all[[k]]*centres)
    # This loop can be a little slow
    # Can make this part faster if necessary
    for (inter in seq_len(n_interactions[k])) {
      cur_inter <- as.numeric(beta_all[[k]][inter+p, ])
      beta_add[interactions[1L, inter], l_start_vec[k]:nlambda] <-
        beta_add[interactions[1L, inter], l_start_vec[k]:nlambda] -
        centres[interactions[2L, inter]]*cur_inter
      beta_add[interactions[2L, inter], l_start_vec[k]:nlambda] <-
        beta_add[interactions[2L, inter], l_start_vec[k]:nlambda] -
        centres[interactions[1L, inter]]*cur_inter
    }
    beta_all[[k]][1L:p, ] <- beta_all[[k]][1L:p, ] + beta_add[, l_start_vec[k]:nlambda]
    # Reset relevant part of beta_add to be zero
    if (k < length(l_start_vec)) zero(beta_add, l_start_vec[k+1L])
  }
  out <- list("call"=this.call,
              "a0"=a0,
              "beta"=beta_all,
              "fitted"=fit_all,
              "lambda"=lambda,
              "nulldev"=null_dev,
              "nobs"=n,
              "nvars"=p,
              "var_indices"=var_names_main,
              "interactions"=interactions,
              "path_lookup"=path_lookup,
              "l_start"=l_start_vec,
              "n_interactions"=n_interactions)
  class(out) <- "BT"
  return(out)
}