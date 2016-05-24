#'Testing pathway risk association
#'
#'This functions computes p-values frm score tests of genetic pathway risk association in 5 different models
#'
#'@param data a \code{data.frame} of \code{N} rows and set up as the output from \code{\link{sim_SCR_data}} with columns:\itemize{
#'\item{\code{XR}:}{ time to recurrence / death / censoring}
#'\item{\code{XD}:}{ time to death / censoring}
#'\item{\code{DeltaR}:}{ Indicator of censoring (0), recurrence (1), or death (2) for this earliest time XR}
#'\item{\code{DeltaD}:}{ Indicator of censoring (0) or death (1)}
#'\item{\code{XPFS}:}{ time to recurrence / death / censoring (=XR)}
#'\item{\code{DeltaPFS}:}{ Indicator of censoring (0) or recurrence or death, whichever came first (1)}
#'\item{\code{Z_1,...,Z_P}:}{ genomic variables}
#'}
#'@param ind_gene columns indices of genes in the pathway of interest. Default is \code{7:ncol(data)}).
#'@param num_perts number of perturbations used. Default is \code{1000}.
#'@param Ws optional inputed perturbations, should be a vector of length \code{N x num_perts} containing
#'i.i.d. realization of a random variable with mean=0 and variance=1.
#'@param rho a vector of rhos, such as one found created from the range returned by \code{\link{findRhoInterval}},
#'used for tuning non-linear kernel. Only used if \code{kernel} is not \code{"linear"}. Default is \code{NA}.
#'Currently not available for use by user-defined kernels.
#'@param kernel a character string indicating which kernel is used. Possible values (currently implemented) are
#'\code{"linear"}, \code{"gaussian"} or \code{"poly"}. Otherwise, this can also be a user defined kernel function.
#'See \code{\link{genericKernelEval}}.
#'@param pca_thres a number between \code{0} and \code{1} giving the threshold to be used for PCA.
#'Default is \code{0.9}. If \code{NULL}, no PCA is performed.
#'@param d if \code{kernel} is \code{"poly"}, the polynomial power. Default is 2 (quadratic kernel).
#'@param get_ptb_pvals a logical flag indicating whether perturbed p-values should be returned
#'as part of the results. Default is \code{FALSE}.
#'@param ... extra parameters to be passed to a user-defined kernel.
#'
#'
#'
#'@return either a \code{vector} of p-values for 5 different models with names:\itemize{
#'\item{\code{"SCR"}:}{ Semi-Comperting Risks}
#'\item{\code{"PFS"}:}{ Progression Free Survival}
#'\item{\code{"CR"}:}{ Competing Risks}
#'\item{\code{"OS"}:}{ Overall Survival}
#'\item{\code{"SCR_alt"}:}{ SCR allowing different tuning params for the two event time processes}
#'} or else if \code{get_ptb_pvals} is \code{TRUE}, a \code{list} with 2 elements:\itemize{
#'\item{\code{"obs_pvals"}:}{ a vector containing the observed p-values for each of the 5 models as described above}
#'\item{\code{"null_pvals_perts"}:}{ a matrix of dimensions \code{num_perts x 5} containing the corresponding
#'perturbed p-values}
#'}
#'
#'@references Neykov M, Hejblum BP, Sinnot JA, Kernel Machine Score Test for Pathway Analysis in the
#'Presence of Semi-Competing Risks, submitted, 2016.
#'@importFrom stats rnorm sd
#'@export
#'
#'@examples
#'
#'## First generate some Data
#'feat_m_fun <- function(X){
#'  sin(X[,1]+X[,2]^2)-1
#'}
#'feat_d_fun <-  function(X){
#'  (X[,4]-X[,5])^2/8
#'}
#'mydata <- sim_SCR_data(data_size = 400, ncol_gene_mat = 20, feat_m = feat_m_fun,
#'                       feat_d = feat_d_fun, mu_cen = 40, cov=0.5)
#'
#' #initial range
#'ind_gene <- c(7:ncol(mydata))
#'my_rho_init <- seq(0.01, 20, length=300)*length(ind_gene)
#'range(my_rho_init)
#'
#'\dontrun{
#'# compute the interval for rho
#'rho_set <- findRhoInterval(tZ=t(mydata[,ind_gene]), rho_init = my_rho_init, kernel="gaussian")
#'rho_set
#'range(my_rho_init) # good to check that the interval produced here is strictly contained in rho_init
#'# otherwise, expand rho.init and rerun
#'
#'rhos <- exp(seq(log(rho_set[1]),log(rho_set[2]), length=50))
#'
#'# run the tests with Gaussian kernel
#'compute_all_tests(data = mydata, num_perts=1000, rho=rhos, kernel="gaussian")
#'# run the tests with linear kernel
#'compute_all_tests(data=mydata, num_perts=1000, kernel="linear")
#'}
#'
compute_all_tests <- function(data, ind_gene=7:ncol(data), num_perts=1000, Ws=NULL, rho = NA,
                              kernel = c("linear", "gaussian", "poly"), d = 2, pca_thres=0.9,
                              get_ptb_pvals=FALSE, ...){
  n <- nrow(data)
  if(is.null(Ws)){
    Ws <- stats::rnorm(n*num_perts)
  }else{
    num_perts <- length(Ws)/n
  }
  p_gene <- length(ind_gene)

  ## data_M sorted by X_M (col #1), data_D sorted by X_D (col #2)
  o_data_M <- order(data[, 1])
  data_M <- data[o_data_M, ]
  o_data_DfromM <- order(data_M[, 2])
  data_D <- data_M[o_data_DfromM, ]

  M_Mc <- M_vec(Inf, data_M[, 1], 1*(data_M[,3]==1), rep(0, 2), matrix(0, n, 2)) ## cause specific hazard for M
  M_Dc <- M_vec(Inf, data_M[, 1], 1*(data_M[,3]==2), rep(0, 2), matrix(0, n, 2)) ## cause specific hazard for D
  M_Dm <- M_vec(Inf, data_D[, 2], 1*(data_D[,4]==1), rep(0, 2), matrix(0, n, 2)) ## marginal hazard for D
  M_Mc <- M_Mc[, o_data_DfromM, drop = FALSE] ## this ensures that all sorting is done according to sorting of data_D (col #2) ##
  M_Dc <- M_Dc[, o_data_DfromM, drop = FALSE] ## this ensures that all sorting is done according to sorting of data_D (col #2) ##

    if(kernel=="linear"){
    K_rho <- linKernelEval(tZ=t(data_D[, ind_gene]))
    K_rho <- matrix(K_rho, nrow=1)
  }else if(kernel=="gaussian"){
    K_rho <- gaussKernelEval_multipleRhos(tZ=t(data_D[, ind_gene]), rho=rho)
  }else if(kernel=="poly"){
    K_rho <- polyKernelEval_multipleRhos(tZ = t(data_D[, ind_gene]), rho=rho, d=d)
  }else{
    K_rho <- genericKernelEval(tZ = t(data_D[, ind_gene]), ...)
    K_rho <- matrix(K_rho, nrow=1)
  }
  sqrt_ncol_K_rho <- sqrt(ncol(K_rho))
  K_rho <- matrix(c(t(K_rho)), ncol = sqrt_ncol_K_rho, byrow = TRUE)

  # do kernel PCA
  tempcomp <- c()
  if(is.null(pca_thres)){
    warning("Not performing PCA: potential loss of power, especially with small samples")
  }else{
    if(kernel=="linear"){
      K_rho_eig <- eigen(K_rho)
      ncomp <- which(cumsum(K_rho_eig$values/sum(K_rho_eig$values)) >= pca_thres)[1]
      K_rho <- tcrossprod(VTM(sqrt(K_rho_eig$values[1:ncomp]), sqrt_ncol_K_rho)*K_rho_eig$vectors[, 1:ncomp])
      tempcomp <- ncomp
    }else{
      K_rho_new_list <- list()
      for (i in 1:length(rho)){
        K_rho_eig <- eigen(K_rho[1:sqrt_ncol_K_rho + (i-1)*sqrt_ncol_K_rho, ])
        ncomp <- which(cumsum(K_rho_eig$values/sum(K_rho_eig$values)) >= pca_thres)[1]
        tempcomp <- c(tempcomp, ncomp)
        K_rho_new_list[[i]] <- tcrossprod(VTM(sqrt(K_rho_eig$values[1:ncomp]), sqrt_ncol_K_rho)*K_rho_eig$vectors[, 1:ncomp])
      }
      K_rho <- do.call(rbind, K_rho_new_list)
      cat(paste("range of eigenvalues considered:", paste(range(tempcomp), collapse=" - ")), "\n")
    }
  }
  t_K_rho <- t(K_rho)

  # returns a to obtaing the K_rhos for a fixed matrix K_rho of the type above
  obtain_K_rhos <- function(t_K_rho, M){
    return(c(M%*%matrix(M%*%t_K_rho, nrow = sqrt_ncol_K_rho)))
  }

  ## including cause specific hazard for M, D as well as marginal for D ##
  stats_all <- cbind("Mc"=obtain_K_rhos(t_K_rho, M_Mc),
                     "Dc"=obtain_K_rhos(t_K_rho, M_Dc),
                     "Dm"=obtain_K_rhos(t_K_rho, M_Dm)
  )

  ## align perturbations with the ordered data
  Ws_mat <- matrix(Ws, n, num_perts)
  Ws_mat_M <- Ws_mat[o_data_M, ]
  Ws_mat_D <- Ws_mat_M[o_data_DfromM, ]

  ### using the new method based on cause-specific hazard model ###
  M_Mc_pert <- M_vec_pert(Ws_mat_M, Inf, data_M[, 1], 1*(data_M[, 3]==1), rep(0, 2), matrix(0, n, 2)) ## cause specific for M
  M_Dc_pert <- M_vec_pert(Ws_mat_M, Inf, data_M[, 1], 1*(data_M[, 3]==2), rep(0, 2), matrix(0, n, 2)) ## cause specific for D
  M_Mc_pert <- M_Mc_pert[o_data_DfromM, ] ## this ensures that sorting is done according to sorting of data_D
  M_Dc_pert <- M_Dc_pert[o_data_DfromM, ] ## this ensures that sorting is done according to sorting of data_D
  M_Dm_pert <- M_vec_pert(Ws_mat_D, Inf, data_D[, 2], 1*(data_D[, 4]==1), rep(0, 2), matrix(0, n, 2)) ## marginal for D

  part_1_Mc <- crossprod(M_Mc_pert, t_K_rho)
  part_1_Dc <- crossprod(M_Dc_pert, t_K_rho)
  part_1_Dm <- crossprod(M_Dm_pert, t_K_rho)

  part_2_Mc <- c(t(M_Mc_pert))*c(part_1_Mc)
  part_2_Dc <- c(t(M_Dc_pert))*c(part_1_Dc)
  part_2_Dm <- c(t(M_Dm_pert))*c(part_1_Dm)

  part_3.list <- list("Mc"=matrix(part_2_Mc, ncol = num_perts, byrow = T), ## cause specific M, n*K_rho rows: (1:n), n+(1:n), 2n+(1:n)
                      "Dc"=matrix(part_2_Dc, ncol = num_perts, byrow = T), ## cause specific D
                      "Dm"=matrix(part_2_Dm, ncol = num_perts, byrow = T)) ## marginal D

  stat_pert_std_list <-  as.list(1:3)
  names(stat_pert_std_list) <- c("Mc", "Dc", "Dm")
  sd_all <- NULL
  for(l in 1:3){
    tmpres = NULL
    for(i in 0:(length(rho) - 1)){ ## num_perts x n.rho matrix
      tmpres <-  cbind(tmpres, colSums(part_3.list[[l]][i*n + 1:n,]))
    }
    tmp_sd <- apply(tmpres, 2, stats::sd)
    stat_pert_std_list[[l]]=as.matrix(tmpres/VTM(tmp_sd, num_perts))
    sd_all <- cbind(sd_all, tmp_sd)
  }
  stats_all_std <- stats_all/sd_all ## standardized statistic for Mc, Dc and Dm
  colnames(stats_all_std) <- c("Mc", "Dc", "Dm")

  ## McDm, allowing different rhos for each process if needed ##
  stats_tune2_std <- sum(apply(stats_all_std[, c(1,3), drop=FALSE], 2, max)) ## maximum statistic across rho and then outcome
  perts_tune2_std <- apply(matrix(unlist(lapply(stat_pert_std_list, function(xx){apply(xx, 1, max)})), ncol=3)[, c(1, 3)], 1, sum)
  pvals_tune2_std <- mean(perts_tune2_std > stats_tune2_std)
  ## sum over cause M and cause D then take max over rho ##
  stats_McDc_std <- max(apply(stats_all_std[,c("Mc","Dc"),drop=F], 1, sum))
  perts_McDc_std <- apply(stat_pert_std_list[["Mc"]] + stat_pert_std_list[["Dc"]], 1, max)
  pvals_McDc_std <- mean(perts_McDc_std > stats_McDc_std)
  ## sum over cause M and cause D then take max over rho (check how much additional info gained by marg D) ##
  stats_McDm_std <- max(apply(stats_all_std[,c("Mc","Dm"),drop=F], 1, sum))
  perts_McDm_std <- apply(stat_pert_std_list[["Mc"]] + stat_pert_std_list[["Dm"]], 1, max)
  pvals_McDm_std <- mean(perts_McDm_std > stats_McDm_std)

  ## take the maximum statistic across rho for each of the outcome, cause M, cause D and marg D
  stats_max_std <- apply(stats_all_std, 2, max) ## maximum statistic across rho
  perts_max_std <- matrix(unlist(lapply(stat_pert_std_list,function(xx){apply(xx, 1, max)})), ncol=3) ## num_perts x 3 matrix
  pvals_max_std <- apply(perts_max_std > VTM(stats_max_std, num_perts), 2, mean)
  #names(pvals_max_std) = c("Mc","Dc","Dm")

  ## using time to min(death,progression) with a single cox model: PFS
  o_data_PFS <- order(data[, 5])
  data_Min <- data[o_data_PFS, ]
  Ws_mat_min <- Ws_mat[o_data_PFS, ]
  M_min <- M_vec(Inf, data_Min[, 5], data_Min[, 6], rep(0, 2), matrix(0, n, 2))

  if(kernel=="linear"){
    K_rho <- linKernelEval(tZ=t(data_Min[, ind_gene]))
    K_rho <- matrix(K_rho,nrow = 1)
  }else if(kernel=="gaussian"){
    K_rho <- gaussKernelEval_multipleRhos(tZ=t(data_Min[, ind_gene]), rho = rho)
  }else if(kernel=="poly"){
    K_rho <- polyKernelEval_multipleRhos(tZ=t(data_Min[, ind_gene]), rho=rho, d=d)
  }else{
    K_rho <- genericKernelEval(tZ = t(data_Min[, ind_gene]), ...)
    K_rho <- matrix(K_rho, nrow=1)
  }
  K_rho <- matrix(c(t(K_rho)), ncol = sqrt_ncol_K_rho, byrow = T)

  # do kernel PCA
  if(kernel=="linear"){
    K_rho_eig <- eigen(K_rho)
    ncomp <- which(cumsum(K_rho_eig$values/sum(K_rho_eig$values)) >= pca_thres)[1]
    K_rho <- tcrossprod(VTM(sqrt(K_rho_eig$values[1:ncomp]), sqrt_ncol_K_rho)*K_rho_eig$vectors[, 1:ncomp])
  }else{
    K_rho_new_list <- list()
    for (i in 1:length(rho)){
      K_rho_eig <- eigen(K_rho[1:sqrt_ncol_K_rho + (i-1)*sqrt_ncol_K_rho, ])
      ncomp <- which(cumsum(K_rho_eig$values/sum(K_rho_eig$values)) >= pca_thres)[1]
      K_rho_new_list[[i]] <- tcrossprod(VTM(sqrt(K_rho_eig$values[1:ncomp]), sqrt_ncol_K_rho)*K_rho_eig$vectors[, 1:ncomp])
    }
    K_rho <- do.call(rbind, K_rho_new_list)
  }
  t_K_rho <- t(K_rho)

  stats_min <- obtain_K_rhos(t_K_rho, M_min)
  M_min_pert <- M_vec_pert(Ws_mat_min, Inf, data_Min[, 5], data_Min[, 6], rep(0, 2), matrix(0, n, 2))

  part_1_min <- crossprod(M_min_pert, t_K_rho)
  part_2_min <- c(t(M_min_pert))*c(part_1_min)
  part_3_min <- matrix(part_2_min, ncol = num_perts, byrow = T)

  res_min_list <- list()
  for(i in 0:(length(rho) - 1)){
    res_min_list[[i+1]] <- colSums(part_3_min[i*n + 1:n, ])
  }
  res_min <- do.call(cbind, res_min_list)

  std_devs_min <- sqrt(colSums((res_min - matrix(rep(colMeans(res_min), num_perts), nrow = nrow(res_min),byrow=T))^2)/(num_perts - 1))
  pval_min <-  sum(apply(t(res_min)/std_devs_min, 2, max) >= max(stats_min/std_devs_min))/num_perts

  res_pvals <- c("SCR"=pvals_McDm_std, "PFS"=pval_min, "CR"=pvals_McDc_std,"OS"=pvals_max_std[3], "SCR_alt"=pvals_tune2_std)

  if(get_ptb_pvals){
    tmpout <- list("obs_pvals" = res_pvals,
                   "pert_pvals" = cbind("SCR"=1-rank(perts_McDm_std)/num_perts, "PFS"=1-rank(apply(t(res_min)/std_devs_min, 2, max))/num_perts,
                                        "CR"=1-rank(perts_McDc_std)/num_perts, "OS"=1-rank(perts_max_std[, 3])/num_perts,
                                        "SCR_alt"=1-rank(perts_tune2_std)/num_perts))
  }else{
    tmpout <- res_pvals
  }


  return(tmpout)
}
