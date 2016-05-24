#'Data Simulation Function
#'
#'@param data_size an integer giving the simulated sample size \code{N}
#'@param ncol_gene_mat an integer giving the simulated number of genomic covariates \code{P}
#'@param feat_m a function that transforms the genomic features into the signal for the metastasis
#'process. This function should a matrix of dimensions \code{N X P} as its only argument.
#'@param feat_d a function that transforms the genomic features into the signal for the death
#'process. This function should a matrix of dimensions \code{N X P} as its only argument.
#'@param mu_cen mean of the exponential censoring process
#'@param cov the correlation between the genomic covariates
#'@param lam_m baseline hazard constant for metastasis process. Default is \code{1/15}.
#'@param lam_d baseline hazard constant for death process. Default is \code{1/20}.
#'@param norm_vcov vector of length 4 of correlation between errors between the two processes on
#'the normal scale before being complementary-log-log-transformed. Default is \code{c(1,.5,.5,1)}.
#'@return a \code{data.frame} with columns:\itemize{
#'\item{\code{XR}:}{ time to recurrence / death / censoring}
#'\item{\code{XD}:}{ time to death / censoring}
#'\item{\code{DeltaR}:}{ Indicator of censoring (0), recurrence (1), or death (2) for this earliest time XR}
#'\item{\code{DeltaD}:}{ Indicator of censoring (0) or death (1)}
#'\item{\code{XPFS}:}{ time to recurrence / death / censoring (=XR)}
#'\item{\code{DeltaPFS}:}{ Indicator of censoring (0) or recurrence or death, whichever came first (1)}
#'\item{\code{Z_1,...,Z_P}:}{ genomic variables}
#'}
#'
#'@importFrom mvtnorm rmvnorm
#'@importFrom stats pnorm rexp
#'@export
#'
#'@examples
#'feat_m_fun <- function(X){
#'  sin(X[,1]+X[,2]^2)-1
#'}
#'feat_d_fun <-  function(X){
#'  (X[,4]-X[,5])^2/8
#'}
#'mydata <- sim_SCR_data(data_size = 400, ncol_gene_mat = 20, feat_m = feat_m_fun,
#'                       feat_d = feat_d_fun, mu_cen = 30, cov=0.5)
#'head(mydata)
#'## how many experience both events
#'mean(mydata[,"DeltaR"]==1 & mydata[,"DeltaD"]==1)
#'## how many only recur
#'mean(mydata[,"DeltaR"]==1 & mydata[,"DeltaD"]==0)
#'## how many only die
#'mean(mydata[,"DeltaR"]==2 & mydata[,"DeltaD"]==1)
#'## how many are censored
#'mean(mydata[,"DeltaR"]==0 & mydata[,"DeltaD"]==0)
#'
#'
sim_SCR_data <- function(data_size, ncol_gene_mat, feat_m, feat_d, mu_cen, cov,
                         lam_m = 1/15, lam_d = 1/20, norm_vcov = c(1,.5,.5,1)){
  # the genomic covariates
  S <- matrix(cov, ncol_gene_mat, ncol_gene_mat)
  diag(S) <- 1
  Zgen.mat <- mvtnorm::rmvnorm(data_size, sigma = S)

  # the correlated errors
  eps <- mvtnorm::rmvnorm(n=data_size, mean=c(0,0), sigma=matrix(norm_vcov,2,2))
  eps <- log(-log(stats::pnorm(eps)))
  Tm <- exp(-log(lam_m) + feat_m(Zgen.mat) + eps[,1])
  Td <- exp(-log(lam_d) + feat_d(Zgen.mat) + eps[,2])

  # the censoring process
  Tc <- stats::rexp(data_size, 1/mu_cen)

  Tmat1 <- cbind(C=as.vector(Tc), Tm=as.vector(Tm), Td=as.vector(Td))
  Tmin1 <- apply(Tmat1, 1, min) # smallest event time
  delta1 <- apply(Tmat1 == Tmin1, 1, which) - 1 # 0 for censoring, 1 for Mets, 2 for death

  Tmat2 <- cbind(as.vector(Tc), as.vector(Td))
  Tmin2 <- apply(Tmat2, 1, min) # smallest of censoring, death
  delta2 <- apply(Tmat2 == Tmin2, 1, which) - 1 # 0 for censoring, 1 for Death

  colnames(Zgen.mat) <- paste("Z", 1:ncol(Zgen.mat), sep="")
  data <- data.frame(XR=Tmin1, XD=Tmin2, DeltaR=delta1, DeltaD=delta2, XPFS=Tmin1, DeltaPFS=as.numeric(delta1 != 0), Zgen.mat)

  return(data)
}