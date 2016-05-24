#' Group screening and selection
#' @description Performs the grouped variable screening and selection.
#' @param X A matrix of grouped predictors.
#' @param y A numeric vector of response.
#' @param formula An object of class "\code{\link{formula}}".
#' @param data An optional data frame.
#' @param group A vector of describing the grouping of the predictors.
#' Numeric and consectutive group indices are recommended.
#' @param threshold A threshold meaning how many groups are screened out. The default is \code{NULL}.
#' See details.
#' @param scale The type of scaling of the predictors. The default is "\code{standardize}".
#' @param criterion The screening criterion. The default is "\code{gSIS}".
#' @param family A description of the error distribution and link function to be used
#' in the model. The default is "\code{gaussian}".
#' @param select A logical value indicating whether to perform the grouped variable selection.
#' The default is \code{FALSE}.
#' @param penalty The penalty to be applied to the screened model. The default is
#' "\code{grSCAD}". Only valid when \code{select = TRUE}.
#' @param cross.validation A logical value indicating whether to perform the k-fold
#' cross-validation when conducting the grouped variable selection. Only valid when
#' \code{select = TRUE}. The default is \code{FALSE}.
#' @param norm The type of norm applied to \code{criterion = gSIS} and
#' \code{criterion = gHOLP}. The default is \code{L1} norm.
#' @param q A quantile for calculating the data-driven threshold in the permutation-based
#' grouped screening. The default is \code{1}. (i.e., the maximum absolute value of the
#' permuted estimates). See details for more information.
#' @param perm.seed A seed of the random number generator used for the permutation-based
#' screening to obtain the threshold. See details.
#' @param nfolds The number of folds to perform the cross-validation. The default is \code{10}.
#' @param cv.seed A seed of the random number generator used for the cross-validation.
#' @param parallel A logical value indicating whether to use the parallel computing. The
#' default is \code{FALSE}.
#' @param cl A cluster object as returned by makeCluster, or the number of nodes to be
#' created in the cluster.
#' @param cores The number of cores to use for parallel execution. If not specified, the number
#' of core is set to be 3.
#' @param ... Optional arguments passed to \code{\link[grpreg]{grpreg}}.
#'
#' @details
#' The grouped variable selection will have big challenges or even fail in the presence of
#' ultra-high dimension of groups. To solve these issues, we implement a two-stage procedure.
#' At the first stage, a grouped screening procedure is
#' applied to reduce the dimensions of groups from ultra-high to moderate or even small one,
#' then we can use the grouped variable selection for the screened data without facing the
#' big challenges at the second stage. At the first stage, the sure screening property
#' ensures that the screening procedure can retain all important groups with overwhelming
#' probability.
#'
#' This function is used to accomplish this two-stage procedure. At the first stage,
#' we apply different screening criteria for grouped variables by calculating the grouped
#' screening values that measures the strength of relationship between response and entire
#' predictors of each group. See \code{\link{grp.criValues}} for the details of calculating
#' the grouped screening values.
#' For the \code{family = "gaussian"} case, we select the groups which
#' have the largest \code{threshold} values of screening criterion indices.
#' On the contrary, for the \code{family = "binomial"} or \code{"poisson"}
#' case, we keep the groups which have the smallest \code{threshold} values of screening criterion
#' values.
#'
#' If \code{threshold = NULL}, we use the random permutation strategy to gain the threshold
#' (\code{threshold}), which is called the data-driven threshold. The details can be seen in
#' Fan, Feng and Song (2011). Larger threshold (\code{threshold}) will lead to larger probability
#' of containing the true important groups, but may result in more intense computation in
#' grouped variable selection and larger false positive rate.
#'
#' At the second stage, we use the function \code{\link[grpreg]{grpreg}} in \code{grpreg}
#' package
#' developed by Patrick Breheny to fit the penalized regression model for the grouped
#' variables that are screened out at the first stage. More details of the grouped variable
#' selection can be refered to the details of \code{\link[grpreg]{grpreg}}.
#'
#' Also, we use the parallel computation in this function by importing the
#' \code{doParallel} package to improve the computation efficiency.
#'
#' @note The missing values are removed before the analysis.
#'
#' @return If \code{select = FALSE}, a list with class "\code{grpss}" containing the following
#' components:
#' \item{call}{The function call.}
#' \item{y}{The response.}
#' \item{X}{The screened predictors.}
#' \item{group.screen}{The indices of screened groups.}
#' \item{threshold}{The threshold.}
#' \item{criterion}{The screening criterion.}
#' If \code{select = TRUE}, a list with class "\code{grpreg}" or "\code{cv.grpreg}"
#' (when \code{cross.validation = TRUE}) containing the similar components as in function
#' \code{\link[grpreg]{grpreg}} or \code{\link[grpreg]{cv.grpreg}}, plus the following
#' three elements:
#' \item{call}{Same as above.}
#' \item{group.screen}{Same as above.}
#' \item{criterion}{Same as above.}
#'
#' @author Debin Qiu, Jeongyoun Ahn
#'
#' @references
#' Fan J, Feng Y, Song R (2011). Nonparametric independence screening in sparse ultra-high
#' dimensional additive models. \emph{Journal of the American Statistical Association}.
#' 106:544-557.
#' @seealso \code{\link[grpreg]{grpreg}}, \code{\link[grpreg]{cv.grpreg}},
#' \code{\link{grp.criValues}}
#' @examples
#' library(MASS)
#' set.seed(23)
#' n <- 30 # sample size
#' p <- 3  # number of predictors in each group
#' J <- 50  # group size
#' group <- rep(1:J,each = 3)  # group indices
#' ##autoregressive correlation
#' Sigma <- 0.6^abs(matrix(1:(p*J),p*J,p*J) - t(matrix(1:(p*J),p*J,p*J)))
#' X <- mvrnorm(n,seq(0,5,length.out = p*J),Sigma)
#' betaTrue <- runif(12,-2,5)
#' mu <- X%*%matrix(c(betaTrue,rep(0,p*J-12)),ncol = 1)
#'
#' # normal distribution
#' y <- mu + rnorm(n)
#'
#' # only conduct screening procedure
#' (gss01 <- grpss(X,y,group)) # gSIS
#'
#' # perform both screening and selection procedures
#' ## use grpss.default with cross-validation
#' gss11 <- grpss(X,y,group,select = TRUE,cross.validation = TRUE)
#' summary(gss11)
#' ## without cross-validation
#' gss12 <- grpss(X,y,threshold = 10,group,select = TRUE,criterion = "gHOLP")
#' summary(gss12)
#'
#' ## binomial distribution
#' y1 <- rbinom(n,1,1/(1 + exp(-mu)))
#' (gss21 <- grpss(X,y1,group, criterion = "gAR2")) # use AIC
#' (gss22 <- grpss(X,y1,group, criterion = "gDC"))  # use gDC
#'
#' ## poisson distribution
#' y2 <- rpois(n,lambda = exp(mu))
#' (gss31 <- grpss(X,y2,group, criterion = "gAR2"))
#' (gss22 <- grpss(X,y2,group, criterion = "gDC"))
#'
#' @importFrom stats quantile
#' @importFrom grpreg grpreg
#' @importFrom grpreg cv.grpreg
#' @importFrom MASS mvrnorm
#' @importFrom doParallel registerDoParallel
#' @importFrom stats complete.cases
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response

#' @export

grpss <- function(...) UseMethod("grpss")
#' @rdname grpss
#' @export
grpss.default <- function(X, y, group, threshold = NULL,
                          scale = c("standardize","normalize","none"),
                          criterion = c("gSIS","gHOLP","gAR2","gDC"),
                          family = c("gaussian","binomial","poisson"), select = FALSE,
                          penalty = c("grSCAD","grLasso","grMCP","gel","cMCP"),
                          cross.validation = FALSE,norm = c("L1","L2","Linf"),q = 1,
                          perm.seed = 1,nfolds = 10,cv.seed = NULL,parallel = FALSE,cl = NULL,
                          cores = NULL, ...) {
  if (class(X) != "matrix") {
    tempX <- try(X <- as.matrix(X), silent = TRUE)
    if (class(tempX)[1] == "try-error")
      stop("'X' must be a matrix or can be coerced to a matrix")
  }
  if (!is.numeric(y))
    stop("'y' must be a numeric vector or a matrix")
  if (parallel)
    registerDoParallel(cl = ifelse(is.null(cl),3,cl),cores)
  type <- match.arg(scale)
  criterion <- match.arg(criterion)
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  norm <- match.arg(norm)
  ok <- complete.cases(X,y)
  if (sum(!ok) > 0)
    warning("Missing values exist and have been removed")
  X <- X[ok,]
  y <- y[ok]
  if (length(group) != ncol(X))
    stop("length of group must be equal to ncol(X)")
  if (is.null(colnames(X)))
    colnames(X) <- paste0("X",group)
  X0 <- g0 <- NULL
  if (any(group == 0)) {
    grp0 <- group == 0
    X0 <- X[,grp0]
    X <- X[,!grp0]
    g0 <- group[grp0]
    group <- group[!grp0]
  }
  X <- XX <- X[,order(group)]
  group <- sort(as.numeric(as.factor(group)))
  grp.values <- grp.criValues(X,y,group,criterion,family,type,norm)
  grp.index <- grp.values[order(grp.values[,2],decreasing = TRUE),]
  if (is.null(threshold)) {
    set.seed(perm.seed)
    grp.values0 <- grp.criValues(X[sample(nrow(X)),],y,group,criterion,family,type,norm)[,2]
    set.seed(NULL)
    thres <- grp.index[,2] > as.numeric(quantile(grp.values0,q))
    threshold <- as.integer(sum(thres))
    if (threshold == 0 || threshold == max(group))
      threshold <- as.integer(length(y)/log(length(y)))
  }
  grp.select <- sort(grp.index[1:threshold,1])
  X <- cbind(X0,XX[,group %in% grp.select])
  if (!select) {
    result <- list(call = match.call(),y = y, X = X,group.screen = c(unique(g0),grp.select),
                   threshold = threshold,criterion = criterion)
    class(result) <- "grpss"
  }
  else {
    grp0 <- table(group[group %in% grp.select])
    group <- rep(1:length(grp0),times = grp0)
    if (cross.validation) {
      if (!is.null(cv.seed))
        set.seed(cv.seed)
      grpfit <- cv.grpreg(X,y,group,family = family, penalty = penalty,
                            nfolds = nfolds,...)
    }
    else {
      grpfit <- grpreg(X,y,group,penalty,family,...)
    }
    result <- c(list(call = match.call(),group.screen = c(g0,grp.select),
                     criterion = criterion),grpfit)
    class(result) <- if (cross.validation) "cv.grpreg" else "grpreg"
  }
  return(result)
}

#' @rdname grpss
#' @export
grpss.formula <- function(formula, data, group, ...) {
  mf <- model.frame(formula = formula, data = data)
  X <- model.matrix(attr(mf, "terms"), data = mf)[,-1]
  y <- model.response(mf)
  colnames(X) <- NULL
  result <- grpss.default(X,y,group,...)
  result$call <- match.call()
  result$formula <- formula
  return(result)
}
