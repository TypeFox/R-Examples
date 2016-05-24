#' Find an interval contraining the rho parameter for a non linear kernel
#'
#'@param tZ a \code{P x N} matrix of genomic covariates (i.e., the usual data array Z transposed)
#'
#'@param rho_init an initial large range of possible rhos, which will be considered to see if they are
#' reasonable tuning parameters for the kernel. Default is \code{seq(0.01, 20, length=300)*P}. See Details.
#'
#'@param kernel character string specifying a nonlinear kernel. Currently supported options are:
#' \code{"gaussian"} or \code{"poly"}
#'
#'@param d if \code{kernel} is \code{"poly"}, the polynomial power (e.g. d=2 for quadratic kernel).
#'Default is \code{NA}.
#'
#'@param rate_range a vector of length 2 indicating the range of alpha in the paper. Default is \code{c(1.5,4)}.
#'
#'@param pca_thres a number between \code{0} and \code{1} giving the threshold to be used for PCA.
#'Default is \code{0.9}. If \code{NULL}, no PCA is performed.
#'
#'@param warning_suppress logical flag. Indicating wether the warnings should be suppress during
#'the linear model fitting step. Default is \code{TRUE}. See details.
#'
#'@details This function will print \code{rho_init} range and the range of valid tuning parameters.
#'If that range butts up against either the upper or lower bound of rho_init, you can rerun this function
#'with a bigger \code{rho_init}.
#'
#'Finding the right tuning parameters includes a step of fitting a linear model which can fail
#'because some tuning parameters yield only one eigenvector. We want to eliminate those tuning parameters,
#'so this is OK. However, in case one want to suppress (numerous) annoying warning messages, use the
#'\code{warning_suppress} argument.
#'
#'@return an upper and lower bound to look for rho
#'@importFrom stats na.omit
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
#'                       feat_d = feat_d_fun, mu_cen = 30, cov=0.5)
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
#'#rhos <- exp(seq(log(rho_set[1]),log(rho_set[2]), length=50))
#'}

findRhoInterval <- function(tZ, rho_init = seq(0.01, 20, length=300)*nrow(tZ), kernel = c("gaussian", "poly"),
                            d = NA, rate_range=c(1.5,4), pca_thres=0.9, warning_suppress=TRUE){

  rho <- rho_init

  if(is.null(rho0)){
    rho0 <- 1
  }
  rho0 <- pca_thres

  if (kernel == "gaussian"){
    G <- gaussKernelEval_multipleRhos(tZ, rho=rho)
  }
  if (kernel == "poly"){
    G <- polyKernelEval_multipleRhos(tZ, rho=rho, d=d)
  }

  slope_eigen <- function(ind, GG = G, tZZ = tZ, rho00 = rho0){
    eigenvals <- eigen(matrix(GG[ind,], ncol(tZZ)), symmetric = TRUE, only.values=TRUE)$values
    log_eig <-  log(eigenvals[1:sum(cumsum(eigenvals)/sum(eigenvals) <= rho00)])
    log_j <-  log(1:length(log_eig))
    slope <- try(-MASS::rlm(log_eig~log_j)$coef[2], silent=TRUE)
    if(inherits(slope, "try-error")){
      slope <- NA
    }
    return(c(slope, length(log_j)))
  }

  if(warning_suppress){
    options(warn = -1)
  }
  slope_vec <-  stats::na.omit(cbind("rho" = rho, "slope" = t(sapply(1:length(rho), slope_eigen))))
  if(warning_suppress){
    options(warn = 0)
  }

  colnames(slope_vec)[2:3] = c("slope", "m_keep")
  rho <- slope_vec[,"rho"]
  return(c(max(c(min(rho), slope_vec[slope_vec[,"slope"] <= rate_range[1], "rho"])),
           min(c(max(rho), slope_vec[slope_vec[, "slope"] >= rate_range[2], "rho"]))
  ))
}
