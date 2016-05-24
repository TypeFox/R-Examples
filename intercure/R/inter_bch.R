# inmost function from ICE package
inmost2 <- function(data, eps = 1e-04){
  lr <- data.frame(data)
  names(lr) <- c("l", "r")
  lr$l <- lr$l * (1 + eps)
  lr.stk <- utils::stack(lr, select = c("l", "r"))
  lr.order <- lr.stk[order(lr.stk[, 1]), ]
  n <- length(lr.order[, 1])
  lr1 <- data.frame(lr.order[-n, 1],
                    lr.order[-1, 1],
                    paste(lr.order[-n, 2],
                          lr.order[-1, 2])
                    )
  names(lr1) <- c("q", "p", "label")
  innermost <- lr1[lr1$label == "l r", ]
  innermost[, 1] <- innermost[, 1] / (1 + eps)
  innermost
}

# Create global variables s, r, m and sm1 for main function
create_sr <- function(L,R){
  eq_int <- inmost2(data.frame(L,R), 0)
  s <- eq_int$q
  r <- eq_int$p
  if(!(max(L[R == Inf]) > max(r[r != Inf]))) {
    stop("There is no event-free case with L_i > r_m")
  }
  m <- length(r[r != Inf])
  if(as.logical(sum(!s %in% L) | sum(!r %in% R)))
    cat("Warning: a member of s or r is not in L or R sets (check dataframe)")
  s <- s
  r <- r
  m <- m
  s0 <- (-Inf)
  sm1 <- c(s0,s)
  sm1 <- sm1[-length(sm1)]
  inmost_list <- list(s = s, r = r, m = m, sm1 = sm1)
  inmost_list
}

# Computes the X matrix for the k-th iteraction, proposed by Shen and Hao
compute_Xk <- function(theta, vector_prob,
                       dados, L, R, covariates,
                       F_hat, inmost_list){
  s <- inmost_list$s
  r <- inmost_list$r
  m <- inmost_list$m

  # Matrix with Intercepts
  Z <- (data.frame(1,dados[,covariates]))
  colnames(Z) <- c("intercept",covariates)

  # Evaluating gamma defined by Shen
  Xk <- exp(-exp(theta %*% t(Z)) %x% sapply(s[1:m], F_hat))
  Xk <- Xk * (1 - exp(-exp(theta %*% t(Z)) %x% vector_prob))
  Xk <- rbind(Xk, exp(-exp(theta %*% t(Z))))

  # Called X_k on Shen (a_ij on Turnbull's article)
  for(i in 1:nrow(dados)){
    Xk[,i] <- as.numeric(L[i] <= s & r <= R[i]) * Xk[,i]
  }
  # Defining  Xk
  Xk <- Matrix::Matrix((1 / colSums(Xk)) * t(Xk))
  return(Xk)
}

# Computes auxiliar matrix C given Xk
compute_C <- function(Xk, inmost_list){
  m <- inmost_list$m
  sm1 <- inmost_list$sm1
  s <- inmost_list$s
  aux_C <- diag(m) * NA
  for (j in 1:m){
    aux_C[,j] <- (s[j] <= sm1[-length(sm1)])
  }
  C <- Matrix::Matrix( (Xk[,1:m]) %*% aux_C)
  return(C)
}

# Expected log-likelihood
log_lik <- function(theta = (1:(1 + length(covariates))) * 0,
                  vector_prob,
                  dados, Xk, C,
                  covariates, F_hat, m){
  Z <- data.frame(1, dados[,covariates])
  colnames(Z) <- c("intercept", covariates)
  l <- as.numeric(exp(theta %*% t(Z)))
  soma_1 <- as.numeric(-t(l) %*% (C[][,1:m] %*% vector_prob[1:m]))
  soma_2 <- sum((Matrix::diag(log(1 - exp(-vector_prob %*% t(l)))
                      %*% (Xk[,1:m]))))
  soma_3 <- as.numeric(-Xk[, (m + 1)] %*% l)
  soma_lvero <- as.numeric(soma_1 + soma_2 + soma_3)
  return(soma_lvero)
}


#####################################################################################################################
#  WARNING: IF THETA BIG ENOUGH -> NANs ARE OBTAINED WITH THE R PRECISION IN Xk[,1:m]/(exp(t(b_kp1%x%p))-1))        #
#####################################################################################################################

# First order derivative for p
d_fp <- function(p, theta, Xk, C, Z, m){
  b_kp1 <- exp(theta %*% t(Z))
  dfp <- b_kp1 %*%
    ((Xk[][,1:m] / (exp(t(b_kp1 %x% p)) - 1))
     - C[][,1:m])
  return(dfp)
}

# Second order derivative for p
dd_fp <- function(p, theta, Xk, C, Z, m){
  b_kp1 <- exp(theta %*% t(Z))
  ddfp <- p * 0
  for(i in 1:nrow(Xk)){
    ddfp <- ddfp - (Xk[][i,1:m] * exp(2 * theta %*% t(Z[i,])) *
                   exp(b_kp1[i] * p)) / ((exp(b_kp1[i] * p) - 1) *
                                           (exp(b_kp1[i] * p) - 1))
  }
  return(ddfp)
}

# Updates v
v_plus <- function(p, tau = (1 / p), dfp, ddfp, sigma){
  return(sum( (1 / (sigma / mean(tau * p)) + p * dfp) /
               (tau - p * ddfp)) /
           sum( (p / (tau - p * ddfp))))
}

# Delta p
delta_p <- function(p, tau = (1 / p), v_plus, dfp, ddfp, sigma){
  return((p * (dfp - v_plus) +
            (1 / (sigma / mean(tau * p)))) /
           (tau - p * ddfp))
}

# Euclidean Distance
norm_vec <- function(x) sqrt(sum(x * x))

# Proposed function to evaluate p's convergence
Fn <- function(p, tau, v, theta, Xk, C, eta, m, Z){
  dfp <- d_fp(p,theta,Xk[],C[],Z, m)
  v1 <- dfp + tau - v
  v2 <- (diag(tau) %*% p) - (1 / eta)
  v3 <- 1 - sum(p)
  return(c(v1, v2, v3))
}

# Backtracking for p maximizing
backtracking <- function(p, delta_p, tau,
                         delta_tau, v, delta_v,
                         theta, Xk, C, sigma, m, Z){
  # Defining eta
  eta <- sigma / mean(tau * p)
  pi <- 0.01
  rho <- 0.5
  psi_max <- min(1, min(- (tau / delta_tau)[delta_tau < 0]))
  psi <- psi_max * 0.99
  while(min(p + psi * delta_p) <= 0) psi <- psi * rho
  while((norm_vec(Fn(as.vector(p + psi * delta_p),
                     as.vector(tau + psi * delta_tau),
                     as.vector(v + psi * delta_v),
                     theta, as.matrix(Xk[]),
                     as.matrix(C[]), eta, m, Z))) >
        ( (1 - pi * psi) * norm_vec(Fn(p,
                                       tau,
                                       v,
                                       theta,
                                       as.matrix(Xk[]),
                                       as.matrix(C[]),
                                       eta, m, Z)))) {
    psi <- psi * rho
  }
  return(psi)
}

# Estimates p maximizing the expected log-likelihood
max_p <- function(p, theta, Xk, eps2=0.001, MAXITER=500, sigma, inmost_list, Z){
  # Stop criteria
  m <- inmost_list$m
  C <- compute_C(Xk, inmost_list)
  CRIT1 <- FALSE
  CRIT2 <- FALSE
  it <- 1

  # First and Second order derivatives
  dfp <- d_fp(p, theta, Xk, C, Z, m)
  ddfp <- dd_fp(p, theta, Xk, C, Z, m)

  # Defining initial tau
  ini_tau <- 1 / p
  tau <- ini_tau

  # Initial v
  v <- 0.1

  while( ( !CRIT1 | !CRIT2 ) & it <= MAXITER) {
    # New v (used on delta_v)
    new_v <- v_plus(p, tau, dfp, ddfp, sigma)

    # Delta_p
    delp <- delta_p(p, tau, new_v, dfp, ddfp, sigma)

    # tau_p
    tau_plus <- new_v - ddfp * delp - dfp

    # Delta_tau, Delta_v
    deltau <- tau_plus - tau
    delv <- new_v - v

    # Backtracking
    psi <- backtracking(p, delp, tau,
                        deltau, v, delv, theta,
                        Xk, C, sigma, m, Z)
    p <- as.vector(p + psi * delp)
    tau <- as.vector(tau + psi * deltau)
    v <- as.vector(v + psi * delv)
    dfp <- d_fp(p, theta, Xk, C, Z, m)
    ddfp <- dd_fp(p, theta, Xk, C, Z, m)

    # Updating convergence parameters
    CRIT1 <- (sum(p * tau) < eps2)
    CRIT2 <- (sqrt(sum( (dfp + tau - v) *
                       (dfp + tau - v))) < eps2)
    it <- it + 1
  }
  return(p)
}

#' Fits promotion time cure rate model for interval censored data
#'
#' \code{inter_bch} returns a list with the estimated parameters \code{par} and
#' their asymptotic covariance matrix \code{mcov}. The list also contains a
#' dummy variable \code{stop_c} assuming 0 if algorithm converged and 1 if a
#' stop criteria ended the process.
#'
#' @param dataset Dataset used to fit the model.
#' @param left Vector containing the last check times before event.
#' @param right Vector containing the first check times after event.
#' @param cov String vector containing the column names to be used on the
#'   cure rate predictor.
#' @param sigma Parameter for the primal-dual interior-point algorithm used on
#'   the maximization process. Default value set to 10.
#' @param crit_theta The effects minimum error for convergence purposes.
#' @param crit_p Minimum error of the non-parametric cumulative distribution function.
#' @param max_n Maximum number of iterations of the ECM algorithm.
#' @param output_files Boolean indicating if text outputs for the estimates and
#' variances should be generated.
#' @return The \code{inter_bch} function returns an list containing the
#'   following outputs:
#'   \item{\code{par}}{estimates of theta parameters.}
#'   \item{\code{mcov}}{estimates for the asymptotic covariance
#'   matrix of theta parameters.}
#'   \item{\code{stop_c}}{stop criteria
#'   indicator assuming 1 when process is stopped for a non-convergence
#'   criteria. Assumes 0 when convergence is reached.}
#' @examples
#' set.seed(3)
#' sample_set <- sim_bch(80)
#'
#' ## few iterations just to check how to use the function
#'
#' inter_bch(sample_set, sample_set$L,
#' sample_set$R, c("xi1","xi2"), max_n = 5)
#'
#' ## precise estimate (computationally intensive)
#' \dontrun{
#'
#' inter_bch(sample_set, sample_set$L, sample_set$R, c("xi1","xi2"))
#' }
#' @export
inter_bch <- function(dataset, left, right,
                      cov, sigma = 10,
                      crit_theta = 0.001, crit_p=0.005,
                      max_n = 100, output_files = FALSE){
  # Initializing output files
  if (output_files) {
    est_file_name <- paste("BCH_Estimates.txt", sep="")
    var_file_name <- paste("BCH_Variances.txt", sep="")
    fileconn <- file(est_file_name, "w")
    write(
      paste("PAR:\t",paste0(c("Intercept",cov), collapse="\t")),
      file=fileconn, append=T, sep="")
  }

  # Check if right >= left
  if(any(right < left)) stop("inserted left > right as input")

  # Check that there are no negative times
  if(any(left<0 | right<0)) stop("there are negative times as inputs on left or right vector")

  # Check if cov is in dataset
  if( any( !(cov %in% names(dataset) )) ) stop("specified covariate name not found on the dataset")

  # Specifying the covariates
  dataset <- as.data.frame(dataset)
  Z <- data.frame(1, dataset[,cov])
  colnames(Z) <- c("intercept",cov)

  # Defining variables s, r, m and sm1
  inmost_list <- create_sr(left,right)
  m <- inmost_list$m
  r <- inmost_list$r

  # Initial vector of probabilities (Turnbull "jumps")
  prob <- 1 / m #Uniform on first iteraction
  p <- rep.int(prob, m)

  # F initial estimator
  F_hat_aux <- data.frame(r[1:m], cumsum(p))
  colnames(F_hat_aux) <- c("time", "cum")
  F_hat <- stats::stepfun(F_hat_aux$time, c(0, F_hat_aux$cum))

  # Initial theta
  theta_k <- c(1:ncol(Z)) * 0

  # Convergence flags
  CONV <- FALSE
  CONV2 <- FALSE

  # Iteration counter
  it <- 1

  # ECM Loop
  while(!CONV | !CONV2){
    # Computing matrix X and C
    if ( (it) %% 10 == 0 ) cat("Iteration:", (it))
    Xk <- compute_Xk(theta_k, p, dataset, left, right, cov, F_hat, inmost_list)
    C <- compute_C(Xk, inmost_list)

    # Obtaining theta and it's variance with MLE
    llk <- function(theta) log_lik(theta,p,dataset,Xk,C,cov,F_hat, m)
    fit_theta <- stats::optim(theta_k, llk, method = "BFGS",
                       control = list(fnscale = -1), hessian = T)
    theta_knew <- fit_theta$par
    theta_var <- solve(-fit_theta$hessian)
    if(output_files) {
      utils::write.table(theta_var, file=var_file_name,
                  row.names=FALSE, col.names=FALSE)
    }

    # Checking theta convergence
    CONV <- (max(abs(theta_knew - theta_k)) < crit_theta)
    if ( (it) %% 10 == 0 )
      cat("\nTheta Max Difference:",max(abs(theta_knew - theta_k)))

    # Updating Xk for new theta
    Xk <- compute_Xk(theta_knew, p, dataset,
                     left, right, cov, F_hat, inmost_list)

    # Computing p vector for new Xk and theta
    novo_p <- max_p(p, theta_knew, Xk,
                    sigma = sigma, inmost_list = inmost_list, Z = Z)

    # Checking p convergence
    CONV2 <- (max(abs(novo_p - p)) < crit_p)
    if ( (it) %% 10 == 0 )
      cat("\np Max Difference:", max(abs(novo_p - p)), "\n\n")

    # Updating parameters
    p <- novo_p
    F_hat_aux <- data.frame(r[1:m],cumsum(p))
    colnames(F_hat_aux) <- c("time","cum")
    F_hat <- stats::stepfun(F_hat_aux$time,c(0,F_hat_aux$cum))
    theta_k <- theta_knew

    # Writing new estimates on file
    if (output_files) {
      write(paste("IT",it,":\t",
                  paste0(theta_k, collapse="\t")),
            file=fileconn, append=T, sep="")
    }




    it <- it + 1

    # Check if reached iteration limit defined by user
    if (it >= max_n) {
      if (output_files) {
        write("WARNING: CONVERGENCE NOT REACHED",
              file = fileconn, append=T, sep=" ")
      }
      cat("\n Convergence criteria not met.
          Estimates given for max_n=", max_n)
      break
    }
  }

  if (output_files) close(fileconn)

  # Check if reached the it limit
  crit_stop <- as.numeric(it >= max_n)

  # Outputs an list with useful metrics
  return(list("par" = theta_k, "p" = p, "mcov"=theta_var, "stop_c" = crit_stop))
}
