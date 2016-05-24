# Fast search for interval containing x
findInterval2 <- function(x,v) {
  n <- length(v)
  if (x < v[1])
    return (0)
  if (x >= v[n])
    return (n)
  i <- 1
  k <- n
  while ({
    j <- (k - i) %/% 2 + i; !(v[j] <= x && x < v[j + 1])
    }
    ) {
    if (x < v[j])
      k <- j
    else
      i <- j + 1
  }
  return (j)
}


#Creates a table with time and respective Nelson-Aalen estimates for Cum. Hazard Function
#including latent variables u_i
nelson_aalen_table <- function(data_set,
                               event_times,
                               delta,
                               beta,
                               covariates,
                               u) {
  factors <- u * exp(beta %*% t(data_set[,covariates]))
  event_times_relev <- unique(sort(event_times[delta == 1]))
  parc_f <- function(t)
    (sum(event_times[delta == 1] == t) / sum(factors[event_times >= t]))
  parcels <- sapply(event_times_relev  , parc_f )
  cum_haz <- cumsum(parcels)
  c_haz_f <- data.frame(event_times_relev, cum_haz)
  colnames(c_haz_f) <- c("time", "hazard")
  return(c_haz_f)
}

#Auxiliar function to Nelson-Aalen estimates table
aux_naalen <- function(tempos, naalen_f, par_cl = NULL) {
  z <- unique(sort(tempos))
  if (!is.null(par_cl)) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("The parallel package is needed for parallel version
           of inter_frailty or inter_frailty_cl. Please install it.",
           call. = FALSE)
    }
    vector_naalen <- parallel::parSapply(par_cl, z, naalen_f)
  } else {
    vector_naalen <- sapply(z, naalen_f)
  }
  piece_aalen <- data.frame(z, vector_naalen)
  colnames(piece_aalen) <- c("time", "hazard")
  return(piece_aalen)
}

# Survival Function given x(0) , x(1), Theta, Beta
# Note: Precision problems may occur on the linear predictors
surv_lam <- function(t, xi_0, xi_1, theta, beta, nelson_aalen_function){
  exp(-(exp(theta %*% (xi_0)) / 2) *
        ((1 - (1 / (1 + 2 * nelson_aalen_function(t) * exp(beta %*% (xi_1)))))))
}

# Inverse F function conditioned on l_i and r_i (Method 2 of Lam et al 2007)
inverse_lam_f <- function(w, l_i, r_i, xi_0, xi_1,
                          theta, beta, naalen_original) {
  surv_left <- surv_lam(l_i, xi_0, xi_1, theta, beta, naalen_original)
  surv_right <- surv_lam(r_i, xi_0, xi_1, theta, beta, naalen_original)
  if (surv_left == surv_right) {
    cat(" Warning! Precision problem on l_inear predictor.
        (surv_left = surv_right)")
    stop
  }
  k <- (1 - w) * surv_left + w * surv_right
  k_l_inha <- -2 * log(k) / ((exp(theta %*% (xi_0)) +
                               2 * log(k)) * 2 * exp(beta %*% (xi_1)))
  if (k_l_inha == 0) {
    cat(" Warning. Division by 0 due to precision problems for w= ", w)
    stop
  }
  naalen_l_i <- naalen_original(l_i)
  naalen_r_i <- naalen_original(r_i)
  a <- (naalen_r_i * l_i - naalen_l_i * r_i) / (naalen_r_i - naalen_l_i)
  b <- (r_i - l_i) / (naalen_r_i - naalen_l_i)
  a + k_l_inha * b
}

# Generates n observations y using the inverse transformation
gera_yh <- function(data_set, l_vector, r_vector, delta,
                    cov_theta, cov_beta, theta, beta, naalen_original) {
  tam <- length(l_vector)
  new_y <- rep(NA, tam)
  intercept <- 1
  x_theta <- data.frame(data_set[,cov_theta])
  x_beta <- data.frame(data_set[,cov_beta])
  u_h <- stats::runif(tam)
  for (i in 1:tam) {
    if ( (delta[i] == 0) | (l_vector[i] == r_vector[i]) ){
      new_y[i] <- l_vector[i]
      }
    else {
      xi_0 <- as.numeric(cbind(intercept, x_theta[i,]))
      xi_1 <- as.numeric(x_beta[i,])
      new_y[i] <- inverse_lam_f(u_h[i], l_vector[i], r_vector[i],
                                xi_0, xi_1, theta, beta, naalen_original)
    }
  }
  return(new_y)
}

# Generates latent variables k for n observations
gera_kh <- function(y_h, data_set, delta,
                    cov_theta, cov_beta, theta, beta, nelson_aalen_f) {
  k_h <- y_h * NA
  intercept <- 1
  xi_0 <- cbind(intercept, data_set[,cov_theta])
  xi_1 <- data_set[,cov_beta]
  num <- exp(as.vector(theta %*% t(xi_0)))
  den <- 2 + 4 * sapply(y_h, nelson_aalen_f) * exp(as.vector(beta %*% t(xi_1)))
  k_h <- stats::rpois(length(num), num / den) + delta
  return(k_h)
}

# Generates latent variables u for n observations
gera_uh <- function(y_h, k_h, data_set, r_vector, delta,
                    cov_beta, beta, nelson_aalen_f) {
  u_h <- y_h * NA
  r_estrela <- max(r_vector[delta == 1])
  xi_1 <- data_set[,cov_beta]
  alpha_gamma <- k_h + delta
  beta_gamma <- 1 / (0.5 + sapply(y_h, nelson_aalen_f) *
                       exp(as.vector(beta %*% t(xi_1))))
  cond_u <- (k_h == 0 | y_h > r_estrela)
  u_h <- ifelse(cond_u, 0, stats::rgamma(length(alpha_gamma),
                                  alpha_gamma, scale=beta_gamma))
  return(u_h)
}


# Returns covariance matrix from (5) of (Lam et al 2007)
var_matrix <- function(sum_var, alpha_matrix) {
  M <- nrow(alpha_matrix)
  alpha_j <- Matrix::colMeans(alpha_matrix)
  matriz_soma_2 <- diag(ncol(alpha_matrix)) * 0
  for (h in 1:M) {
    matriz_soma_2 <- matriz_soma_2 +
      ( (alpha_matrix[h,] - alpha_j) %*%
         t(alpha_matrix[h,] - alpha_j) ) / (M - 1)
  }
  sigma_alpha <- (sum_var / M) + (1 + 1 / M) * matriz_soma_2
  return(sigma_alpha)
}

# Convergence criteria returning TRUE of FALSE
convergence_lam <- function(alpha_new, alpha_old, tol = 0.001) {
  conv <- FALSE
  new_val <- c(alpha_new)
  old_val <- c(alpha_old)
  max_error <- max(abs(new_val - old_val))
  if (max_error < tol) conv <- TRUE
  return(conv)
}


#' Fits cure rate frailty model for interval censored data
#'
#' \code{inter_frailty} returns a list with the estimated parameters \code{par}
#' and their asymptotic covariance matrix \code{mcov}. The list also contains a
#' dummy variable \code{stop_c} assuming 0 if algorithm converged and 1 if a
#' stop criteria ended the process.
#'
#' @param dataset Dataset used to fit the model.
#' @param left Vector containing the last check times before event.
#' @param right Vector containing the first check times after event.
#' @param delta Flag vector indicating failure inside interval.
#' @param cov_theta String vector containing the column names to be used on the
#'   cure rate predictor.
#' @param cov_beta String vector containing the column names to be used on the
#'   predictor associated with the hazard function.
#' @param M Number of replicates generated by each iteration on the ANDA
#'   (Asymptotic Normal Data Augmentation) algorithm.
#' @param b Parameter for initial theta and beta variances.
#' @param tol Numeric for tolerance of convergence.
#' @param max_n Maximum number of algorithm's iterations without the burn
#'   in.
#' @param par_cl Registered SOCK cluster for parallel process. If NULL (default)
#'   the program loops are executed sequentially.
#' @param burn_in Number of burn in iterations.
#' @param output_files Boolean indicating if text outputs for the estimates and
#' variances should be generated.
#' @return The \code{inter_frailty} function returns an list containing the
#'   following outputs:
#'   \item{\code{par}}{estimates of theta and beta parameters.}
#'   \item{\code{mcov}}{estimates for the covariance matrix of theta and beta
#'   parameters.}
#'   \item{\code{stop_c}}{stop criteria indicator assuming 1 when process is
#'   stopped for a non-convergence criteria. Assumes 0 when convergence is
#'   reached.}
#' @examples
#' ## few iterations just to check how to use the function
#' set.seed(3)
#' sample_set <- sim_frailty(80)
#'
#' inter_frailty(sample_set, sample_set$L, sample_set$R, sample_set$delta,
#' c("xi1","xi2"), c("xi1","xi2"), M = 10, max_n = 3, burn_in = 0)
#'
#' ## precise estimate (computationally intensive)
#' \dontrun{
#'
#' inter_frailty(sample_set, sample_set$L, sample_set$R, sample_set$delta,
#' c("xi1"), c("xi2"), M = 50, max_n = 50, burn_in = 10)
#' }
#' @export
inter_frailty <- function(dataset, left, right, delta,
                          cov_theta, cov_beta,
                          M, b=0.001, tol = 0.001, max_n=100,
                          par_cl = NULL,
                          burn_in=30, output_files = FALSE) {
# Defining output variables
  if (output_files) {
    est_file_name <- paste("LAM_Estimates.txt", sep="")
    var_file_name <- paste("LAM_Variances.txt", sep="")
    fileconn <- file(est_file_name, "w")
    write(
      paste("PAR:\t",paste0(c("Intercept",cov_theta,cov_beta), collapse="\t")),
      file=fileconn, append=T, sep="")
  }

  # Checks delta
  if(any(!(delta %in% c(0,1)))) stop("delta vector should contain 0 or 1 only")

  # Check if right >= left
  if(any(right < left)) stop("inserted left > right as input")

  # Check that there are no negative times
  if(any(left<0 | right<0)) stop("there are negative times as inputs on left or right vector")

  # Check if cov is in dataset
  if( any( !(cov_theta %in% names(dataset) )) ) stop("specified incidence covariate name not found on the dataset")
  if( any( !(cov_beta %in% names(dataset) )) ) stop("specified latency covariate name not found on the dataset")

  # Initial values for y
  y_nxm <- u_nxm <- matrix(NA, nrow=M, ncol = nrow(dataset))
  for (i in 1:M) y_nxm[i,] <- ifelse (delta == 1,
                                    (left + right) / 2,
                                    left)

  # Initial values for the parameters
  compr_theta <- 1 + length(cov_theta); compr_beta <- length(cov_beta)
  compr_alpha <- compr_theta + compr_beta
  a_M <- matrix(NA, nrow=M, ncol=compr_alpha)
  lbls <- c("intercept", cov_theta, cov_beta)
  colnames(a_M) <- lbls
  alpha <- c(1:compr_alpha) * 0
  sigma_alpha <- b * diag(compr_alpha)
  beta <- alpha[(compr_theta + 1):compr_alpha]

  # Initial values for latent vector u
  u <- delta

  # Initial Nelson-Aalen estimator
  Vetores_NAalen <- nelson_aalen_table(dataset, y_nxm[1,],
                                       delta, beta, cov_beta, u)
  mean_naalen <- stats::stepfun(Vetores_NAalen$time,
                                c(0,Vetores_NAalen$hazard))

  # Initializing convergence criteria and parameters
  conv <- FALSE; n <- 0
  a_M_NEW <- a_M

  # Iterative process (with parallel computing)
  while(!conv | n <= burn_in) {
    if ( (n + 1) %% 10 == 0 ) cat("Iteration:", (n + 1),"\n")
    #iter_time <- system.time({
    if(!is.null(par_cl)){
      list_reg <- foreach::foreach(iterators::icount(M),
                          .packages=c("MASS","Matrix","survival"),
                          .export=c("surv_lam","inverse_lam_f", "gera_yh",
                                    "gera_kh", "gera_uh"),
                          .inorder=F) %dopar% {
        a_M <- MASS::mvrnorm(n=1, alpha, sigma_alpha)
        theta_M <- a_M[1:compr_theta]
        beta_M <- a_M[(compr_theta + 1):compr_alpha]
        y <- gera_yh(dataset, left, right, delta,
                     cov_theta, cov_beta, as.numeric(theta_M),
                     as.numeric(beta_M), mean_naalen)
        k <- gera_kh(y, dataset, delta, cov_theta, cov_beta,
                     theta_M, beta_M, mean_naalen)
        u <- gera_uh(y, k, dataset, right, delta, cov_beta,
                     beta_M, mean_naalen)

        # Poisson Regression for Theta
        o_set <- k * 0 - log(2)
        expression_theta <- paste("dataset$", cov_theta[1:length(cov_theta)],
                                  sep = "", collapse="+")
        formula_theta <- stats::formula(paste0("k~",expression_theta,"+offset(o_set)"))
        fit_theta <- stats::glm(formula_theta, family = stats::poisson)

        # Cox Regression for Beta
        expression_beta <- paste("dataset$", cov_beta[1:length(cov_beta)] ,
                                 sep = "", collapse="+")
        formula_beta <- stats::formula(paste0("Surv(y,delta)~",expression_beta,"+ offset(ifelse(log(u)==-Inf, -200,log(u)))"))
        fit_beta <- survival::coxph(formula_beta, method = 'breslow')

        # Outputs of Parallel Computing
        out <- list(fit_theta$coef, stats::vcov(fit_theta),
                    fit_beta$coef, stats::vcov(fit_beta), y, u)
        out
      }
    } else {
      list_reg <- foreach::foreach(iterators::icount(M),
                          .packages=c("MASS","Matrix","survival"),
                          .export=c("surv_lam","inverse_lam_f", "gera_yh",
                                    "gera_kh", "gera_uh"),
                          .inorder=F) %do% {
        a_M <- MASS::mvrnorm(n=1, alpha, sigma_alpha)
        theta_M <- a_M[1:compr_theta]
        beta_M <- a_M[(compr_theta + 1):compr_alpha]
        y <- gera_yh(dataset, left, right, delta,
                     cov_theta, cov_beta, as.numeric(theta_M),
                     as.numeric(beta_M), mean_naalen)
        k <- gera_kh(y, dataset, delta,
                     cov_theta, cov_beta,
                     theta_M, beta_M, mean_naalen)
        u <- gera_uh(y, k, dataset, right, delta,
                     cov_beta, beta_M, mean_naalen)

        # Poisson Regression for Theta
        o_set <- k * 0 - log(2)
        expression_theta <- paste("dataset$", cov_theta[1:length(cov_theta)],
                                  sep = "", collapse="+")
        formula_theta <- stats::formula(paste0("k~",expression_theta,"+offset(o_set)"))
        fit_theta <- stats::glm(formula_theta, family = stats::poisson)

        # Cox Regression for Beta
        expression_beta <- paste("dataset$", cov_beta[1:length(cov_beta)] ,
                                 sep = "", collapse="+")
        formula_beta <- stats::formula(paste0("Surv(y,delta)~",expression_beta,"+offset(ifelse(log(u)==-Inf, -200,log(u)))"))
        fit_beta <- survival::coxph(formula_beta, method = 'breslow')


        # Outputs of Parallel Computing
        out <- list(fit_theta$coef, stats::vcov(fit_theta),
                    fit_beta$coef, stats::vcov(fit_beta), y, u)
        out
      }
    }



    # Allocating M new and auxiliary parameter vectors
    sum_var_theta <- sum_var_beta <- 0
    for (h in 1:M) {
      sum_var_theta <- sum_var_theta + list_reg[[h]][[2]]
      sum_var_beta <- sum_var_beta + list_reg[[h]][[4]]
      a_M_NEW[h,1:compr_theta] <- list_reg[[h]][[1]]
      a_M_NEW[h, (compr_theta + 1):compr_alpha] <- list_reg[[h]][[3]]
      y_nxm[h,] <- list_reg[[h]][[5]]
      u_nxm[h,] <- list_reg[[h]][[6]]
    }

    # Matrix of the M beta vectors
    beta_M <- matrix(a_M_NEW[, (compr_theta + 1):compr_alpha], nrow=M)

    # Obtaining new Nelson-Aalen estimator for Cum. Hazard function
    if (!is.null(par_cl)) {
      step_list <- foreach::foreach(h=1:M, .export="nelson_aalen_table",
                           .inorder=F) %dopar% {
        V_NAalen <- nelson_aalen_table(dataset, y_nxm[h,],
                                       delta, beta_M[h,], cov_beta,
                                       u_nxm[h,])
        step_list <- stats::stepfun(V_NAalen$time,
                                    c(0, V_NAalen$hazard))
        step_list
      }
    } else {
      step_list <- foreach::foreach(h=1:M, .export="nelson_aalen_table",
                           .inorder=F) %do% {
        V_NAalen <- nelson_aalen_table(dataset, y_nxm[h,],
                                       delta, beta_M[h,], cov_beta,
                                       u_nxm[h,])
        step_list <- stats::stepfun(V_NAalen$time,
                                    c(0, V_NAalen$hazard))
        step_list
      }
    }

    expression <- paste("step_list[[", 1:M,"]](x)", sep = "", collapse = "+")
    eq4 <- paste("mean_naalen <- function(x) (",expression,")/", M)
    eval(parse(text = eq4))

    # Creating new times/survival table and a more efficient estimator
    V_NAalen <- aux_naalen(sort(y_nxm[, delta == 1]), mean_naalen, par_cl)
    new_mean_naalen <- stats::stepfun(V_NAalen$time,
                                      c(0, V_NAalen$hazard))

    #Calculating the new covariance matrix
    SUM_VAR <- as.matrix(Matrix::bdiag(list(sum_var_theta, sum_var_beta)))
    cov_matrix <- var_matrix(SUM_VAR, a_M_NEW)
    sigma_alpha <- cov_matrix

    #New vector of estimates
    alpha_new <- Matrix::colMeans(a_M_NEW)
    #})
    #print(iter_time)


    #Checking convergence
    conv <- convergence_lam(alpha_new ,alpha, tol)

    #Setting new alpha as old one for iteractive process
    alpha <- alpha_new

    #Writing alpha values
    if (output_files) {
      write(paste("IT",n + 1,":\t",
                  paste0(alpha, collapse="\t")),
            file=fileconn, append=T, sep="")
      utils::write.table(cov_matrix, file=var_file_name,
                  row.names=FALSE, col.names=FALSE)
    }

    # Setting new baseline cum. hazard estimator as old one for iteractive
    # process
    mean_naalen <- new_mean_naalen

    #Updating the iteration counter
    n <- n + 1

    #Checking if iteration counter reached max_n
    if (n == (max_n + burn_in)) {
      if (output_files) {
        write("Warning: Iteration Number achieved but convergence criteria
              not met.", file=fileconn, append=T, sep=" ")
        close(fileconn)
      }
      cat("Convergence criteria not met.
          Estimates given for max_n=", max_n)
      cat("\n")
      break
    }
  }
  #Kills the parallel proccess
  crit_stop <- as.numeric(n == (max_n + burn_in))
  alpha_list <- list(par = alpha, mcov = cov_matrix,
                    stop_c = crit_stop)
  if (output_files) close(fileconn)
  return(alpha_list)
}
