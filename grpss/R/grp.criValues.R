#' Compute group screening criterion values
#' @description Computes the screening criterion values for each group.
#' @param X A matrix of grouped predictors.
#' @param y A numeric vector of response.
#' @param group A vector of group indices for each predictor. Numeric and consecutive group indices are
#' recommended.
#' @param criterion The group screening criterion. The default is \code{gSIS}.
#' @param family A description of the error distribution and link function to be used
#' in the model. The default is \code{gaussian}.
#' @param scale The type of scaling of the predictors. The default is "\code{standardize}".
#' @param norm The type of norm for "\code{gSIS}" or "\code{gHOLP}" screening criterion.
#' See \code{\link{norm_vec}} for details. The default is \code{L1}.
#'
#' @details In the group screening procedure, we first have to calculate the values which
#' measure the strength of relationship between entire predictors of each group and response. These values can
#' be used to screen out the important grouped variables (equivalently, remove the
#' unimportant grouped variables) so that we can reduce the dimension of data from high or
#' ultra-high to moderate or even small one.
#'
#' In greater details, let \eqn{X = (x_{11},x_{12},...,x_{1p_1},...,x_{j1},x_{j2},...,
#' x_{jp_j},...,x_{J1},x_{J2},...,x_{Jp_J})} be the grouped predictors,
#' where \eqn{J} is the number of groups and
#' \eqn{p_j} is the number of predictors in the \eqn{j}-th group.
#'
#' For the case in which
#' \code{family = "gaussian"}, four approaches are applied to calculate
#' such criterion values.
#'
#' The first criterion is "\code{gSIS}" that is the grouped version of sure
#' independence screening [SIS, Fan and Lv (2008)] and defined as
#' \deqn{\hat{w} = X^{T}y = (w_{11},w_{12},...,w_{1p_1},...,w_{j1},w_{j2},...,
#' w_{jp_j},...,w_{J1},w_{J2},...,w_{Jp_J}).}
#' Then we take the norm of the vector \eqn{(w_{j1},w_{j2},...,
#' w_{jp_j})} from the \eqn{j}-th group divided by its size \eqn{p_j}, defined as \eqn{W_j}
#' and thus we obtain the criterion values for the whole groups defined as
#' \deqn{\hat{W} = (W_1,...,W_J).} The details of \code{norm} type can be seen in
#' \code{\link{norm_vec}}.
#'
#' The second criterion is "\code{gHOLP}" that is a grouped version of High-dimensional
#' Ordinary Least-squares Projector [HOLP, Wang and Leng (2015)] and defined as
#' \deqn{\hat{\beta} = X^{T}(XX^{T})^{-1}y = (\beta_{11},\beta_{12},...,\beta_{1p_1},...,
#' \beta_{j1},\beta_{j2},...,\beta_{jp_j},...,\beta_{J1},\beta_{J2},...,\beta_{Jp_J})}
#' and then we proceed the same way as "\code{gSIS}" to incorporate the group structure.
#'
#' The third criterion is "\code{gAR2}" which is called groupwise adjusted r.squared. The
#' basic idea is that we fit a linear model for each group separately and compute the
#' adjusted r.squared that measures the correlation between each group and response. Note
#' that in order to calculate the adjusted r.squared, the maximum group size
#' \eqn{\max(p_j),j=1,...,J} should not be larger than sample size \eqn{n}.
#'
#' The last criterion is "\code{gDC}" which is called grouped distance correlation.
#' The distance correlation [Szekely, Rizzo and Bakirov (2007)] measures the dependence
#' between two random variables or two random vectors.
#' Thus, similar to the idea of "\code{gAR2}", we compute the distance correlation between
#' each group and response. It is worthwhile pointing out that distance correlation can not only
#' measure the linear relationship, but also nonlinear relationship. However, it may take
#' longer time in computation due to the three steps of calculating distance correlation.
#' The distance correlation has been applied to screen the individual variables
#' as in Li, Zhong and Zhu (2012).
#'
#' For the case in which \code{family = "binomial"} and \code{family = "poisson"}, a different
#' screening criterion is used for computing the relationship between response and
#' predictors in each group. To measure the strength of relationship between predictors and
#' response, the Akaike's Information Criterion (AIC) is utilized and defined as
#' \deqn{AIC = -2*LogLikelihood + 2*npar,} where \eqn{LogLikelihood} is the log-likelihood for
#' a fitted generalized linear model, and \eqn{npar} is the number of parameters in the
#' fitted model. In this case, \eqn{npar} is the number of variables within each group,
#' i.e., \eqn{npar = p_j, j = 1,...,J}.
#'
#' Note that the individual "SIS", "HOLP" can be regarded as a special case of "\code{gSIS}",
#' and "\code{gHOLP}" when each group has only one predictor.
#'
#' @return A numeric matrix with two columns: the first column is the group index, and the
#' second column is the grouped screening criterion values corresponding to the first column.
#' @author Debin Qiu, Jeongyoun Ahn
#' @references
#' Fan, J. and Lv, J. (2008). Sure independence screening for ultrahigh dimensional feature
#' space (with discussion). \emph{Journal of the Royal Statistical Society B}, 70, 849-911.
#'
#' Li, R., Zhong,W., and Zhu, L. (2012). Feature screening via distance correlation learning.
#' \emph{Journal of American Statistical Association}, 107, 1129-1139.
#'
#' Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007), Measuring and Testing Dependence
#' by Correlation of Distances, \emph{Annals of Statistics}, Vol. 35 No. 6, pp. 2769-2794.
#'
#' Wang, X. and Leng, C. (2015). High-dimensional Ordinary Least-squares Projector for
#' screening variables.\emph{Journal of the Royal Statistical Society: Series B}. To appear.
#' @seealso \code{\link{grpss}}
#' @examples
#' library(MASS)
#' n <- 30 # sample size
#' p <- 3  # number of predictors in each group
#' J <- 50 # number of groups
#' group <- rep(1:J,each = 3)  # group indices
#' Sigma <- diag(p*J)  # covariance matrix
#' X <- mvrnorm(n,seq(0,5,length.out = p*J),Sigma)
#' beta <- runif(12,-2,5)  # coefficients
#' y <- X%*%matrix(c(beta,rep(0,p*J-12)),ncol = 1) + rnorm(n)
#'
#' grp.criValues(X,y,group)  # gSIS
#' grp.criValues(X,y,group,"gHOLP")  # gHOLP
#' grp.criValues(X,y,group,"gAR2")   # gAR2
#' grp.criValues(X,y,group,"gDC")    # gDC
#'
#' @importFrom stats aggregate
#' @importFrom MASS ginv
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom stats glm
#' @importFrom stats lm
#' @export
grp.criValues <- function(X,y,group,criterion = c("gSIS","gHOLP","gAR2","gDC"),
                          family = c("gaussian","binomial","poisson"),
                          scale = c("standardize","normalize","none"),
                          norm = c("L1","L2","Linf")) {
  if (length(group) != NCOL(X))
    stop("length of group does not match columns of X")
  criterion <- match.arg(criterion)
  family <- match.arg(family)
  norm <- match.arg(norm)
  type <- match.arg(scale)
  X <- Scale(X,type)
  i <- NULL
  if (family == "gaussian") {
    grp.values <- switch(criterion,
                         gSIS = aggregate(abs(t(X)%*%y),list(group),norm_vec,norm)[,2],
                         gHOLP = aggregate(abs(ginv(X)%*%y),list(group),norm_vec,norm)[,2],
                         gAR2 = foreach(i = 1:max(group), .combine = cbind) %dopar%
                           summary(lm(y ~ as.matrix(X[,group == i]) - 1))$adj.r.squared,
                         gDC = foreach(i = 1:max(group), .combine = cbind) %dopar%
                           distcor(as.matrix(X[,group == i]),y))
  }
  else {
    if (any(criterion == c("gSIS","gHOLP")))
      stop("use 'gAR2' or 'gDC' for family = 'binomial' or 'poisson'")
    grp.values <- switch(criterion,
                         gAR2 = foreach(i = 1:max(group), .combine = cbind) %dopar%
                           -glm(y ~ X[,group ==i], family = family)$aic,
                         gDC = foreach(i = 1:max(group), .combine = cbind) %dopar%
                           distcor(X[,group == i],y))
  }
  result <- matrix(c(1:max(group), as.numeric(grp.values)), ncol = 2)
  colnames(result) <- c("group","grp.values")
  return(result)
}
