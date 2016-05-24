# Generates a vector of n latent variables ksi
gera_ksih_effect <- function(y_h, data_set, delta, cov_theta, cov_beta,
                             grp, theta, beta, gamma, nelson_aalen_function){
  m <- length(unique(grp))
  a_h <- b_h <- rep(NA, m)
  for (i in 1:m) {
    a_h[i] <- exp(gamma) + sum(delta[(grp == unique(grp)[i])])
    b_h[i] <- exp(gamma) +
      sum((exp(theta %*%
                 t(cbind(1,data_set[grp == unique(grp)[i],
                                    cov_theta]))) / 2) *
            (1 - 1 / (2 * exp(beta %*%
                                t(data_set[grp == unique(grp)[i],
                                           cov_beta])) *
                        (sapply(y_h[grp == unique(grp)[i]],
                                nelson_aalen_function)) + 1)))
  }
  stats::rgamma(m, a_h, b_h)
}

# Generates a vector of n latent variables k
gera_kh_effect <- function(y_h, ksi_h, data_set, delta,
                           cov_theta, cov_beta, grp, theta, beta,
                           nelson_aalen_function) {
  corresp <- data.frame(clu = unique(grp), effect = ksi_h)
  ksi_ij <- corresp$effect[base::match(grp, corresp$clu)]
  intercept <- 1
  xi_0 <- cbind(intercept, data_set[,cov_theta])
  xi_1 <- data_set[,cov_beta]
  num <- exp(as.vector(theta %*% t (xi_0))) * ksi_ij
  den <- 2 + 4 * sapply(y_h, nelson_aalen_function) *
    exp(as.vector(beta %*% t(xi_1)))
  k_h <- stats::rpois(length(num), num / den) + delta
  return(list(k = k_h, ksi = ksi_ij))
}

# Log-likelihood of gamma
log_vero_gamma <- function(gamma_par, ksi_h){
  logvero <- sum((stats::dgamma(ksi_h,
                         shape=exp(gamma_par),
                         rate=exp(gamma_par),
                         log=TRUE)))
  return(-logvero)
}


# Log-likelihood of theta given y_h, k_h, covariates and theta
log_vero_theta <- function(theta, cov_theta, y_h, k_h, data_set){
  total <- 0
  intercept <- 1
  cov <- cbind(intercept, data_set[,cov_theta])
  total <- sum((-exp(theta %*% t(cov)) / 2) + k_h * (theta %*% t(cov)))
  return(total)
}

# Log-likelihood of beta given y_h, u_h, covariates and beta
log_vero_beta <- function(beta, cov_beta, y_h, u_h, data_set){
  data_set_completa <- data_set
  data_set_completa$u <- u_h
  delta <- data_set_completa$delta
  factors <- data_set_completa$u * exp(beta %*% t(data_set_completa[,cov_beta]))
  total <- 0
  for(i in 1:length(y_h)) {
    xi_1 <- data_set_completa[i,cov_beta]
    soma_risco_i <- sum(factors[y_h >= y_h[i]])
    if (delta[i] == 1) parcel <- delta[i] *
      (sum(beta * xi_1) - log(soma_risco_i))
    else parcel <- 0
    total <- total + parcel
  }
  return(total)
}


# w may need more precision (rmpfr)

surv_cl <- function(tempo, i, j,
                    theta, beta, gamma,
                    cov_theta, cov_beta, grp,
                    nelson_aalen_function, dataset,
                    partial_times,
                    partial_l,
                    partial_r,
                    naalen_mod) {

  w <- exp(gamma)
  base_cl <- dataset[grp == i,]
  etas_i <- as.numeric(exp(theta %*% t(cbind(1, base_cl[,cov_theta]))))
  mus_i <- as.numeric(exp(beta %*% t(base_cl[,cov_beta])))

  if (j == 1) {
    surv_1_cl_i <- (1 + (etas_i[1] / (2 * w)) *
                      (1 - 1 / (2 * mus_i[1] *
                                  naalen_mod(tempo, partial_l[j], partial_r[j]) + 1))) ^ (-w)
    return(surv_1_cl_i)
  }

  num <- (etas_i[j] / 2) *
    (1 - (1 / (2 * mus_i[j] * naalen_mod(tempo, partial_l[j], partial_r[j]) + 1)))
  den <- w + sum(base_cl$delta[1:(j - 1)] * (etas_i[1:(j - 1)] / 2) *
                   (1 - (1 / (2 * mus_i[1:(j - 1)] *
                                (mapply(naalen_mod, partial_times[1:(j - 1)],
                                        partial_l[1:(j - 1)], partial_r[1:(j - 1)])) +
                                1))), na.rm = TRUE)
  surv_j_cl_i <- (1 + num / den) ^ (-w - sum(base_cl$delta[1:(j - 1)]))
  if (is.nan(surv_j_cl_i) | surv_j_cl_i == -Inf) {
    stop("Estimated survival is NaN or -Inf")
  }
  return(surv_j_cl_i)
}

gera_y_cl <- function(u_gen, l, r,
                          i, j, theta, beta, gamma,
                          cov_theta, cov_beta, grp,
                          nelson_aalen_function, dataset,
                          partial_times,
                          partial_l, partial_r){

  # consider how to use this to generate y (each lambda envolved?)
  naalen_mod <- function(t_gen, l_i, r_i) {
    res <- (nelson_aalen_function(l_i) * (r_i - t_gen) +
       nelson_aalen_function(r_i) * (t_gen - l_i)) / (r_i - l_i)
#     if(res < 0 ) { # | is.nan(res)
#       stop ("naalen_mod < 0") # or NaN
#     }
    res
  }

  s_cl_mod <- function(t_gen) surv_cl(t_gen, i, j,
                                      theta, beta, gamma,
                                      cov_theta, cov_beta, grp,
                                      nelson_aalen_function, dataset,
                                      partial_times,
                                      partial_l,
                                      partial_r,
                                      naalen_mod)
  w <- exp(gamma)
  base_cl <- dataset[grp == i,]
  etas_i <- as.numeric(exp(theta %*% t(cbind(1, base_cl[,cov_theta]))))
  mus_i <- as.numeric(exp(beta %*% t(base_cl[,cov_beta])))
  n_l <- nelson_aalen_function(l)
  n_r <- nelson_aalen_function(r)


  s_l <- s_cl_mod(l)

  s_r <- s_cl_mod(r)

  k_1 <- (s_l - s_r)*(1 - u_gen) + s_r

  if (j ==1) {
    k_2 <- (etas_i[1] / ((etas_i[1] - (k_1^(-1/w) - 1) * 2 * w) * 2 * mus_i[1])) - 1/(2 * mus_i[1])
  } else {
    k_1_5 <- (k_1 ^ ( 1 / (-w - sum(base_cl$delta[1:(j - 1)]) ) ) - 1)*(w + sum(base_cl$delta[1:(j - 1)] * (etas_i[1:(j - 1)] / 2) *
                                                                                   (1 - (1 / (2 * mus_i[1:(j - 1)] *
                                                                                                (mapply(naalen_mod, partial_times[1:(j - 1)],
                                                                                                        partial_l[1:(j - 1)],
                                                                                                        partial_r[1:(j - 1)])) +
                                                                                                1))), na.rm = TRUE))
    k_2 <- (etas_i[j] / ((etas_i[j] - 2*k_1_5)*2*mus_i[j])) - 1/(2*mus_i[j])
  }

  gen_time <- (k_2 * (r - l) - n_l*r + n_r*l) / (n_r - n_l)
  return(gen_time)
}

#Generates a vector of n observations using the previous function
gera_yh_effect2 <- function(left, right, delta, cov_theta, cov_beta,
                           grp, theta, beta, gamma,
                           nelson_aalen_function, dataset){
  new_yh <- NA
  for(i in unique(grp)) {
    times_cl_i <- rep(NA, nrow(dataset[grp == i,]))
    delta_cl_i <- delta[grp == i]
    l_cl_i <- left[grp == i]
    r_cl_i <- right[grp == i]
    for(j in c(1:length(times_cl_i))) {
      if(delta_cl_i[j] == 0) {
        times_cl_i[j] <- l_cl_i[j]
      } else {
        if (l_cl_i[j] == r_cl_i[j]) {
          times_cl_i[j] <- l_cl_i[j]
        } else {
          u_var <- stats::runif(1)

          times_cl_i[j] <- gera_y_cl(u_var, l_cl_i[j], r_cl_i[j],
                                     i, j,
                                     theta, beta, gamma,
                                     cov_theta,
                                     cov_beta,
                                     grp,
                                     nelson_aalen_function,
                                     dataset, times_cl_i,
                                     l_cl_i,
                                     r_cl_i)

          if(times_cl_i[j] > r_cl_i[j] | times_cl_i[j] < l_cl_i[j]) {
            cat("\n l_cl_i", l_cl_i[j])
            cat("\n r_cl_i", r_cl_i[j])
            cat("\n gen", times_cl_i[j])
            stop("generated time out of interval")
          }
        }
      }
    }
    new_yh <- c(new_yh, times_cl_i)
  }
  new_yh <- new_yh[-1]
  return(new_yh)
}






#' Cure rate frailty model for interval censored clustered data
#'
#' \code{inter_frailty_cl} returns a list with the estimated parameters
#' \code{par} and their covariance matrix \code{mcov}. The list also contains a
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
#' @param grp Vector containing cluster identifier (numeric or string).
#' @param M Number of replicates generated by each iteration on the ANDA
#'   (Asymptotic Normal Data Augmentation) algorithm.
#' @param b Parameter for initial theta and beta variances.
#' @param tol Numeric for tolerance of convergence.
#' @param max_n Maximum number of algorithm's iterations without the burn in.
#' @param par_cl Registered SOCK cluster for parallel process. If NULL (default)
#'   the program loops are executed sequentially.
#' @param burn_in Number of burn in iterations.
#' @param output_files Boolean indicating if text outputs for the estimates and
#'   variances should be generated.
#' @return The \code{inter_frailty_cl} function returns an list containing the
#'   following outputs: \item{\code{par}}{estimates of theta and beta
#'   parameters.} \item{\code{mcov}}{estimates for the covariance matrix of
#'   theta and beta parameters.} \item{\code{stop_c}}{stop criteria indicator
#'   assuming 1 when process is stopped for a non-convergence criteria. Assumes
#'   0 when convergence is reached.}
#' @examples
#' set.seed(3)
#' sample_set <- sim_frailty_cl(80, nclus = 3)
#'
#' ## few iterations just to check how to use the function
#'
#' inter_frailty_cl(sample_set, sample_set$L, sample_set$R, sample_set$delta,
#' c("xi1"), c("xi1"), grp = sample_set$clus, M = 20,
#' max_n = 2, burn_in = 0)
#'
#' ## precise estimate (computationally intensive)
#' \dontrun{
#'
#' inter_frailty_cl(sample_set, sample_set$L, sample_set$R, sample_set$delta,
#' c("xi1", "xi2"), c("xi1", "xi2"), grp = sample_set$clus, M = 20,
#' max_n = 30, burn_in = 10)
#' }
#' @export
inter_frailty_cl <- function(dataset, left, right, delta, cov_theta, cov_beta,
                             grp, M, b = 0.001, tol = 0.001, max_n=100,
                             par_cl = NULL,
                             burn_in=50, output_files = FALSE) {
  # arranges dataset by clusters
  order_aux <- order(grp)
  dataset <- as.data.frame(dataset[order_aux,])
  left <- left[order_aux]
  right <- right[order_aux]
  delta <- delta[order_aux]
  grp <- grp[order_aux]

  if (output_files) {
    est_file_name <- paste("LAM_Estimates.txt", sep="")
    var_file_name <- paste("LAM_Variances.txt", sep="")
    fileconn <- file(est_file_name, "w")
    write(paste("PAR:\t",paste0(c("Intercept", cov_theta, cov_beta,
                                  "log(gamma)"), collapse="\t")),
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

  #Initial values for y
  y_n_M <- k_n_M <- u_n_M <- ksi_n_M <- matrix(NA, nrow = M,
                                               ncol = nrow(dataset))
  for (i in 1:M) y_n_M[i,] <- ifelse (delta == 1, (left + right) / 2, left)

  #Auxiliar matrix for different clusters
  m <- length(unique(grp))
  z <- matrix(0, ncol = m, nrow = nrow(dataset))
  for (j in 1:m) {
    z[,j] <- ifelse (grp == unique(grp)[j], 1, 0)
  }

  #Initial values for the parameters
  compr_theta <- 1 + length(cov_theta)
  compr_beta <- length(cov_beta)
  compr_alpha <- compr_theta + compr_beta + 1
  a_M <- matrix(NA, nrow = M, ncol = compr_alpha)
  lbls <- c("intercept", cov_theta, cov_beta, "log_w")
  colnames(a_M) <- lbls
  alpha <- c(1:compr_alpha) * 0
  sigma_alpha <- b * diag(compr_alpha)
  beta <- alpha[(compr_theta + 1):(compr_alpha - 1)]

  #Initial values for latent vector u
  u <- delta

  #Initial Nelson-Aalen estimator
  Vetores_NAalen <- nelson_aalen_table(dataset, y_n_M[1,],
                                       delta, beta, cov_beta, u)
  naalen_avg <- stats::stepfun(Vetores_NAalen$time,
                               c(0, Vetores_NAalen$hazard))

  #Initializing convergence criteria and parameters
  conv <- FALSE
  n <- 0
  a_M_NEW <- a_M

  #Iterative process (with parallel computing)
  while(!conv | n <= burn_in) {
    if ( (n + 1) %% 10 == 0 ) cat("Iteration:", (n + 1),"\n")
    #iter_time <- system.time({
    if(!is.null(par_cl)){
      list_reg <- foreach::foreach(iterators::icount(M),
                          .packages=c("MASS","Matrix",
                                      "survival","stats4"),
                          .export=c("surv_cl",
                                    "gera_yh_effect2",
                                    "gera_ksih_effect", "gera_kh_effect",
                                    "log_vero_gamma", "gera_uh"),
                          .inorder=F) %dopar% {
        a_M <- MASS::mvrnorm(n=1, alpha, sigma_alpha)
        theta_M <- a_M[1:compr_theta]
        beta_M <- a_M[(compr_theta + 1):(compr_alpha - 1)]
        gamma_M <- a_M[compr_alpha]
        y <- as.numeric(gera_yh_effect2(left, right, delta, cov_theta,
                                          cov_beta, grp, theta_M,
                                          beta_M, gamma_M,
                                          naalen_avg, dataset))
        ksi <- gera_ksih_effect(y, dataset, delta,
                                cov_theta, cov_beta, grp,
                                theta_M, beta_M, gamma_M,
                                naalen_avg)
        k_aux <- gera_kh_effect(y, ksi, dataset, delta,
                                cov_theta, cov_beta, grp,
                                theta_M, beta_M,
                                naalen_avg)
        k <- k_aux$k # vetorzao de k para algum h
        ksi_geral <- k_aux$ksi # vetor de efeitos de grupo
        u <- gera_uh(y, k, dataset, right, delta, cov_beta, beta_M, naalen_avg) # checked!

        #Gamma Regression for w
        fit_gamma <- stats4::mle(log_vero_gamma,
                                 start = list(gamma_par = as.numeric(gamma_M)),
                                 fixed = list(ksi_h = ksi))

        #Poisson Regression for Theta
        o_set <- k * 0 - log(2) + log(ksi_geral)
        expression_theta <- paste("dataset$",
                                  cov_theta[1:length(cov_theta)] ,
                                  sep = "", collapse="+")
        formula_theta <- stats::formula(paste0("k~",expression_theta,"+offset(o_set)"))
        fit_theta <- stats::glm(formula_theta, family = stats::poisson)

        #Cox Regression for Beta
        expression_beta <- paste("dataset$",
                                 cov_beta[1:length(cov_beta)] ,
                                 sep = "", collapse="+")
        formula_beta <- stats::formula(paste0("Surv(y,delta)~",expression_beta,"+offset(ifelse(log(u)==-Inf, -200,log(u)))"))
        fit_beta <- survival::coxph(formula_beta, method = 'breslow')

        #Outputs of Parallel Computing
        out <- list(fit_theta$coef, stats::vcov(fit_theta),
                    fit_beta$coef, stats::vcov(fit_beta),
                    stats4::coef(fit_gamma)[1], stats4::vcov(fit_gamma),
                    y, u, k, ksi_geral)
        out
      }
    } else {
      list_reg <- foreach::foreach(iterators::icount(M),
                          .packages=c("MASS","Matrix","survival",
                                      "stats4"),
                          .export=c("surv_cl",
                                    "gera_yh_effect2",
                                    "gera_ksih_effect",
                                    "gera_kh_effect","log_vero_gamma",
                                    "gera_uh"),
                          .inorder=F) %do% {
        a_M <- MASS::mvrnorm(n=1, alpha, sigma_alpha)
        theta_M <- a_M[1:compr_theta]
        beta_M <- a_M[(compr_theta + 1):(compr_alpha - 1)]
        gamma_M <- a_M[compr_alpha]
        y <- as.numeric(gera_yh_effect2(left, right, delta, cov_theta,
                                          cov_beta, grp, theta_M,
                                          beta_M, gamma_M,
                                          naalen_avg, dataset))
        ksi <- gera_ksih_effect(y, dataset, delta,
                                cov_theta, cov_beta, grp,
                                theta_M, beta_M, gamma_M,
                                naalen_avg)
        k_aux <- gera_kh_effect(y, ksi, dataset, delta,
                                cov_theta, cov_beta, grp,
                                theta_M, beta_M,
                                naalen_avg)
        k <- k_aux$k
        ksi_geral <- k_aux$ksi
        u <- gera_uh(y, k, dataset, right, delta, cov_beta, beta_M, naalen_avg)

        #Gamma Regression for w
        fit_gamma <- stats4::mle(log_vero_gamma,
                                 start = list(gamma_par = as.numeric(gamma_M)),
                                 fixed = list(ksi_h = ksi))




        #Poisson Regression for Theta
        o_set <- k * 0 - log(2) + log(ksi_geral)
        expression_theta <- paste("dataset$",
                                  cov_theta[1:length(cov_theta)],
                                  sep = "",
                                  collapse="+")
        formula_theta <- stats::formula(paste0("k~",expression_theta,"+offset(o_set)"))
        fit_theta <- stats::glm(formula_theta, family = stats::poisson)



        #Cox Regression for Beta
        expression_beta <- paste("dataset$",
                                 cov_beta[1:length(cov_beta)],
                                 sep = "",
                                 collapse="+")
        formula_beta <- stats::formula(paste0("Surv(y,delta)~",expression_beta,"+offset(ifelse(log(u)==-Inf, -200,log(u)))"))
        fit_beta <- survival::coxph(formula_beta, method = 'breslow')

        #Outputs of Parallel Computing
        out <- list(fit_theta$coef, stats::vcov(fit_theta), fit_beta$coef,
                    stats::vcov(fit_beta), stats4::coef(fit_gamma)[1], stats4::vcov(fit_gamma),
                    y, u, k, ksi_geral)
        out
      }
    }

    #Allocating M new and auxiliary parameter vectors
    sum_var_gamma <- sum_var_theta <- sum_var_beta <- 0
    for (h in 1:M) {
      sum_var_theta <- sum_var_theta + list_reg[[h]][[2]]
      sum_var_beta <- sum_var_beta + list_reg[[h]][[4]]
      a_M_NEW[h,1:compr_theta] <- list_reg[[h]][[1]]
      a_M_NEW[h, (compr_theta + 1):(compr_alpha - 1)] <- list_reg[[h]][[3]]
      a_M_NEW[h,compr_alpha] <- list_reg[[h]][[5]]
      sum_var_gamma <- sum_var_gamma + list_reg[[h]][[6]]
      y_n_M[h,] <- list_reg[[h]][[7]]
      u_n_M[h,] <- list_reg[[h]][[8]]
      k_n_M[h,] <- list_reg[[h]][[9]]
      ksi_n_M[h,] <- list_reg[[h]][[10]]
    }

    #Matrix of the M beta vectors
    beta_M <- matrix(a_M_NEW[, (compr_theta + 1):(compr_alpha - 1)], nrow=M)

    #Obtaining new Nelson-Aalen estimator for Cum. Hazard function
    if(!is.null(par_cl)){
      step_list <- foreach::foreach(h=1:M,
                                    .export="nelson_aalen_table",
                                    .inorder=F) %dopar% {
        V_NAalen <- nelson_aalen_table(dataset, y_n_M[h,],
                                       delta, beta_M[h,],
                                       cov_beta, u_n_M[h,])
        step_list <- stats::stepfun(V_NAalen$time,
                                    c(0, V_NAalen$hazard))
        step_list
      }
    } else {
      step_list <- foreach::foreach(h=1:M,
                                    .export="nelson_aalen_table",
                                    .inorder=F) %do% {
        V_NAalen <- nelson_aalen_table(dataset, y_n_M[h,],
                                       delta, beta_M[h,],
                                       cov_beta, u_n_M[h,])
        step_list <- stats::stepfun(V_NAalen$time,
                                    c(0, V_NAalen$hazard))
        step_list
      }
    }

    expression <- paste("step_list[[", 1:M,"]](x)", sep = "", collapse="+")
    eq4 <- paste("naalen_avg<-function(x) (",expression,") / ", M)
    eval(parse(text = eq4))

    # Creating new times/survival table and a more efficient estimator
    V_NAalen <- aux_naalen(sort(y_n_M[, delta == 1]), naalen_avg, par_cl)
    naalen_avg_new <- stats::stepfun(V_NAalen$time,
                                     c(0, V_NAalen$hazard))

    #Calculating the new covariance matrix
    SUM_VAR <- as.matrix(Matrix::bdiag(list(sum_var_theta,
                                    sum_var_beta,
                                    sum_var_gamma)))
    cov_matrix <- var_matrix(SUM_VAR, a_M_NEW)
    sigma_alpha <- cov_matrix

    #New vector of estimates
    alpha_new <- Matrix::colMeans(a_M_NEW)


    #Checking convergence
    conv <- convergence_lam(alpha_new, alpha, tol)

    #Setting new alpha as old one for iteractive process
    alpha <- alpha_new
    #cat("\n alpha:", alpha)
    #})
    #print(iter_time)

    #Writing alpha values
    if (output_files) {
      write(paste("IT",n + 1,":\t",
                  paste0(alpha, collapse="\t")),
            file=fileconn, append=T, sep="")
      utils::write.table(cov_matrix, file=var_file_name,
                  row.names=FALSE, col.names=FALSE)
    }

    # Setting new baseline cum. hazard estimator as
    # old one for iteractive process
    naalen_avg <- naalen_avg_new

    #Updating the iteration counter
    n <- n + 1

    #Checking if iteration counter reached max_n
    if(n == (max_n + burn_in)) {
      if (output_files) {
        write("\nWarning: Iteration Number achieved but
              convergence criteria not met.",file=fileconn,append=T, sep=" ")
      }
      cat("\n Convergence criteria not met. Estimates given
          for max_n=", max_n)
      cat("\n")
      break
    }
  }
  crit_stop <- as.numeric(n == (max_n + burn_in))
  alpha_list <- list(par = alpha, mcov = cov_matrix, stop_c = crit_stop)
  if (output_files) close(fileconn)
  return(alpha_list)
}
