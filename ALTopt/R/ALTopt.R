#' Transform the array to the model matrix
#'
#' The internal function to make the model matrix corresponded to linear
#'   predictor model from the array (vector) containg coordinates of stress factors.
#'
#' @keywords internal
#' @importFrom stats model.matrix
ExtendedForm <- function(array, formula, nf) {
  terms <- attr(terms(formula), "term.labels")
  mtx <- as.data.frame(matrix(array, ncol = nf))
  colnames(mtx) <- terms[1:nf]
  out <- model.matrix(formula, mtx)
}

#' Calculates the prediction variance at the particular use condition
#'
#' The internal function to calculate the prediction variance
#'
#' @keywords internal
PreVar <- function(location, formula, nf, infMtxInv) {
  # Calculates the prediction variance in a particular use condition
  use <- ExtendedForm(location, formula, nf)
  as.numeric(use %*% infMtxInv %*% t(use))
}

#' Perform the k-means clustering and make design table
#'
#' The internal function to perform the k-means clustering and make design table.
#'
#' @keywords internal
#' @importFrom stats kmeans
kmeansCls <- function(Mtx, nCls) {
  kmeansOut <- kmeans(Mtx, nCls)
  Tbl <- cbind(kmeansOut$centers, kmeansOut$size)
  colnames(Tbl)[ncol(Tbl)] <- paste("allocation")
  Tbl
}

#' Objective function of D optimal design with right censoring
#'
#' The internal function to calculate the objective function value of
#' D optimal design with right censoring plan.
#'
#' @keywords internal
Dobj.rc <- function(x, formula, coef, nf, tc, alpha) {
  X <- ExtendedForm(x, formula, nf)
  b <- coef
  eta <- X %*% b #linear predictor
  phi <- 1 - exp(- exp(eta) * tc ^ alpha)
  W <- diag(phi[, 1])
  XWX <- t(X) %*% W %*% X
  det(XWX)
}

#' Objective function of U optimal design with right censoring
#'
#' The internal function to calculate the objective function value of
#' U optimal design with right censoring plan.
#'
#' @keywords internal
Uobj.rc <- function(x, formula, coef, nf, tc, alpha, useCond) {
  X <- ExtendedForm(x, formula, nf)
  b <- coef
  eta <- X %*% b #linear predictor
  phi <- 1 - exp(- exp(eta) * tc ^ alpha)
  W <- diag(phi[, 1])
  XWX <- t(X) %*% W %*% X
  c <- try(qr.solve(XWX), silent = TRUE)
  if (class(c) == "try-error")
    return("cannot calculated ; information matrix is near singular")
  else
    PreVar(location = useCond, formula = formula, nf = nf, infMtxInv = c)
}

#' Objective function of I optimal design with right censoring
#'
#' The internal function to calculate the objective function value of
#' I optimal design with right censoring plan.
#'
#' @keywords internal
Iobj.rc <- function(x, formula, coef, nf, tc, alpha, useLower, useUpper) {
  X <- ExtendedForm(x, formula, nf)
  b <- coef
  eta <- X %*% b #linear predictor
  phi <- 1 - exp(- exp(eta) * tc ^ alpha)
  W <- diag(phi[, 1])
  XWX <- t(X) %*% W %*% X
  c <- try(qr.solve(XWX), silent = TRUE)
  if (class(c) == "try-error")
    return("cannot calculated ; information matrix is near singular")
  else {
    # numerical integration
    intgratedPV <- cubature::adaptIntegrate(PreVar, lowerLimit = useLower,
                                            upperLimit = useUpper, formula = formula,
                                            nf = nf, infMtxInv = c)$integral
    volume <- 1
    for (i in 1:nf) volume <- volume * (useUpper[i] - useLower[i])
    intgratedPV / volume
  }
}

#' Design evaluation with right censoring.
#'
#' \code{\link{alteval.rc}} calculates the objective function value
#'   (D, U or I) for a given design with right censoring plan.
#'
#' @param designTable a data frame containing the coordinates and the number of
#'   allocation of each design point. The design created by either
#'   \code{\link{altopt.rc}} or \code{\link{altopt.ic}} or any design matrix
#'   with the same form as those can be provided for this argument.
#' @param optType the choice of \code{"D"}, \code{"U"} and \code{"I"} optimality.
#' @param tc the censoring time.
#' @param nf the number of stress factors.
#' @param alpha the value of the shape parameter of Weibull distribution.
#' @param formula the object of class formula which is the linear predictor model.
#' @param coef the numeric vector containing the coefficients of each term in \code{formula}.
#' @param useCond the numeric vector of use condition.
#'   It should be provided when \code{optType} is \code{"U"}. The length of the vector
#'   should be same as the number of stress factors.
#' @param useLower the numeric vector of lower bound of use region.
#'   It should be provided when \code{optType} is \code{"I"}. The length of the vector
#'   should be same as the number of stress factors.
#' @param useUpper the numeric vector of upper bound of use region.
#'   It should be provided when \code{optType} is \code{"I"}. The length of the vector
#'   should be same as the number of stress factors.
#' @return The objective function value corresponded by \code{optType}
#'   for a given design with right censoring plan.
#' @seealso \code{\link{altopt.rc}}
#' @examples
#' # Evaluation of factorial design for right censoring.
#' x1 <- c(0, 1, 0, 1)
#' x2 <- c(0, 0, 1, 1)
#' allocation <- c(25, 25, 25, 25)
#' facDes <- data.frame(x1, x2, allocation)
#'
#' alteval.rc(facDes, "D", 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01))
#'
#' alteval.rc(facDes, "U", 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useCond = c(1.758, 3.159))
#'
#' alteval.rc(facDes, "I", 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#' @export
alteval.rc <- function(designTable, optType, tc, nf, alpha, formula, coef,
                       useCond, useLower, useUpper) {

  # Transform design to the single column array.
  x <- NULL
  for (col in (1:nf)) {
    for (row in (1:nrow(designTable))) {
      x <- c(x, rep(designTable[row, col], designTable[row, nf + 1]))
    }
  }

  if      (optType == "D")
    value <- Dobj.rc(x, formula, coef, nf, tc, alpha)
  else if (optType == "U")
    value <- Uobj.rc(x, formula, coef, nf, tc, alpha, useCond)
  else if (optType == "I")
    value <- Iobj.rc(x, formula, coef, nf, tc, alpha, useLower, useUpper)
  else stop('Wrong optimization criteria')
  value
}



#' Optimal design with right censoring.
#'
#' \code{\link{altopt.rc}} creates D, U or I optimal design
#' of the accelerated life testing with right censoring plan.
#'
#' @param optType the choice of \code{"D"}, \code{"U"} and \code{"I"} optimality.
#' @param N the number of test units.
#' @param tc the censoring time.
#' @param nf the number of stress factors.
#' @param alpha the value of the shape parameter of Weibull distribution.
#' @param formula the object of class formula which is the linear predictor model.
#' @param coef the numeric vector containing the coefficients of each term in \code{formula}.
#' @param useCond the numeric vector of use condition.
#'   It should be provided when \code{optType} is \code{"U"}. The length of the vector
#'   should be same as the number of stress factors.
#' @param useLower the numeric vector of lower bound of use region.
#'   It should be provided when \code{optType} is \code{"I"}. The length of the vector
#'   should be same as the number of stress factors.
#' @param useUpper the numeric vector of upper bound of use region.
#'   It should be provided when \code{optType} is \code{"I"}. The length of the vector
#'   should be same as the number of stress factors.
#' @param nOpt the number of repetition of optimization process. Default is 1.
#' @param nKM the number of repetition of k-means clustering. Default is 20.
#' @param nCls the number of clusters used for k-means clustering. If not specified,
#'   it is set as the number of parameters in the linear predictor model.
#' @return A list with components
#' \itemize{
#'   \item{call:}{ the matched call.}
#'   \item{opt.design.ori:}{ the original optimal design.}
#'   \item{opt.value.ori:}{ the objective function value of \code{opt.design.ori}.}
#'   \item{opt.design.rounded:}{ the optimal design clustered by rounding in third decimal points.}
#'   \item{opt.value.rounded:}{ the objective function value of \code{opt.design.rounded}.}
#'   \item{opt.design.kmeans:}{ the optimal design clustered by \code{\link[stats]{kmeans}}.}
#'   \item{opt.value.kmeans:}{ the objective function value of \code{opt.design.kmeans}.}
#' }
#' @references
#' {
#' Monroe, E. M., Pan, R., Anderson-Cook, C. M., Montgomery, D. C. and
#' Borror C. M. (2011) A Generalized Linear Model Approach to Designing
#' Accelerated Life Test Experiments, \emph{Quality and Reliability Engineering
#'  International} \bold{27(4)}, 595--607
#'
#' Yang, T., Pan, R. (2013) A Novel Approach to Optimal Accelerated Life Test
#'  Planning With Interval Censoring, \emph{Reliability, IEEE Transactions on}
#'   \bold{62(2)}, 527--536
#' }
#' @seealso \code{\link[stats]{kmeans}}, \code{\link{alteval.rc}}
#' @examples
#' \dontrun{
#' # Generating D optimal design for right censoring.
#' altopt.rc("D", 100, 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01))
#'
#' # Generating U optimal design for right censoring.
#' altopt.rc("D", 100, 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useCond = c(1.758, 3.159))
#'
#' # Generating I optimal design for right censoring.
#' altopt.rc("D", 100, 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useLower = c(1.458, 2.859),
#' useUpper = c(2.058, 3.459))
#' }
#' @importFrom stats aggregate
#' @importFrom stats optim
#' @importFrom stats runif
#' @export
altopt.rc <- function(optType, N, tc, nf, alpha, formula, coef,
                      useCond, useLower, useUpper,
                      nOpt = 1, nKM = 30, nCls = NULL) {

  # function for optimization
  Opt.rc <- function() {
    xInit <- runif(nf * N, min = 0, max = 1)
    lb <- rep(0, nf * N)
    ub <- rep(1, nf * N)
    if      (optType == "D")
      solution <- optim(xInit, Dobj.rc, NULL, formula, coef, nf, tc, alpha,
                        method = "L-BFGS-B", lower = lb, upper = ub,
                        control = list(fnscale = -1) # Maximization
      )
    else if (optType == "U")
      solution <- optim(xInit, Uobj.rc, NULL, formula, coef, nf, tc, alpha,
                        useCond,
                        method = "L-BFGS-B", lower = lb, upper = ub)
    else if (optType == "I")
      solution <- optim(xInit, Iobj.rc, NULL, formula, coef, nf, tc, alpha,
                        useLower, useUpper,
                        method = "L-BFGS-B", lower = lb, upper = ub)
    else stop('Wrong optimization criteria')
    solution
  }

  # Repeat optimization with different initial points
  for (i in 1:nOpt) {
    Curr_sol <- Opt.rc()
    if (i == 1) Best_sol <- Curr_sol
    else if (optType == "D" && Curr_sol$value > Best_sol$value)
      Best_sol <- Curr_sol
    else if ((optType %in% c("U", "I")) && Curr_sol$value < Best_sol$value)
      Best_sol <- Curr_sol
  }

  terms <- attr(terms(formula), "term.labels")
  optDesignMtx <- as.data.frame(matrix(Best_sol$par, ncol = nf))
  colnames(optDesignMtx) <- terms[1:nf]

  # Original result

  columnlist <- apply(optDesignMtx, 2, list)
  columnlist <- lapply(columnlist, unlist)
  optDesignOri <- aggregate(optDesignMtx, columnlist, length)
  while (length(optDesignOri) != (nf + 1))
    optDesignOri <- optDesignOri[, -length(optDesignOri)]
  names(optDesignOri)[ncol(optDesignOri)] <- paste("allocation")
  optValueOri <- alteval.rc(optDesignOri, optType, tc, nf, alpha,
                            formula, coef, useCond, useLower, useUpper)

  # Rounding result
  optDesignMtxRound <- round(optDesignMtx, digits = 3)
  columnlist <- apply(optDesignMtxRound, 2, list)
  columnlist <- lapply(columnlist, unlist)
  optDesignRound <- aggregate(optDesignMtxRound, columnlist, length)
  while (length(optDesignRound) != (nf + 1))
    optDesignRound <- optDesignRound[, -length(optDesignRound)]
  names(optDesignRound)[ncol(optDesignRound)] <- paste("allocation")
  optValueRound <- alteval.rc(optDesignRound, optType, tc, nf, alpha,
                              formula, coef, useCond, useLower, useUpper)


  # k-means result
  if (is.null(nCls)) nCls <- length(coef)
  for (j in 1:nKM) {
    curDesignKmeans <- kmeansCls(optDesignMtx, nCls)
    curValueKmeans <- alteval.rc(curDesignKmeans, optType, tc, nf, alpha,
                                 formula, coef, useCond, useLower, useUpper)
    if (j == 1) {
      optDesignKmeans <- curDesignKmeans
      optValueKmeans <- curValueKmeans
    } else if (optType == "D" && curValueKmeans > optValueKmeans) {
      optDesignKmeans <- curDesignKmeans
      optValueKmeans <- curValueKmeans
    } else if ((optType %in% c("U", "I")) && curValueKmeans < optValueKmeans) {
      optDesignKmeans <- curDesignKmeans
      optValueKmeans <- curValueKmeans
    }
  }

  # Creates output
  out <- list(call = match.call(),
              opt.design.ori = optDesignOri,
              opt.value.ori = as.numeric(optValueOri),
              opt.design.rounded = optDesignRound,
              opt.value.rounded = as.numeric(optValueRound),
              opt.design.kmeans = optDesignKmeans,
              opt.value.kmeans = ifelse (class(optValueKmeans) == "character",
                                     optValueKmeans, as.numeric(optValueKmeans))
  )
  out
}




#' Objective function of D optimal design with interval censoring
#'
#' The internal function to calculate the objective function value of
#' D optimal design with interval censoring plan.
#'
#' @keywords internal
Dobj.ic <- function(x, formula, coef, nf, t, k, alpha) {
  X <- ExtendedForm(x, formula, nf)
  b <- coef
  dt <- t / k
  eta <- X %*% b # linear predictor
  temp1 <- exp(2 * eta)
  temp2 <- diag(temp1[, 1])
  sum <- 0
  for(j in 1:k) {
    c <- ((j - 1) ^ alpha - j ^ alpha) * dt ^ alpha
    temp3 <- c ^ 2 * exp(- exp(eta) * (dt ^ alpha) * (j ^ alpha))
    temp4 <- temp3 / (1 - exp(c * exp(eta)))
    sum = sum + temp4
  }
  W <- sum[, 1] * temp2 # W plays a role of weight matrix
  XWX <- t(X) %*% W %*% X
  det(XWX)
}

#' Objective function of U optimal design with interval censoring
#'
#' The internal function to calculate the objective function value of
#' U optimal design with interval censoring plan.
#'
#' @keywords internal
Uobj.ic <- function(x, formula, coef, nf, t, k, alpha, useCond) {
  X <- ExtendedForm(x, formula, nf)
  b <- coef
  dt <- t / k
  eta <- X %*% b # linear predictor
  temp1 <- exp(2 * eta)
  temp2 <- diag(temp1[, 1])
  sum <- 0
  for(j in 1:k) {
    c <- ((j - 1) ^ alpha - j ^ alpha) * dt ^ alpha
    temp3 <- c ^ 2 * exp(- exp(eta) * (dt ^ alpha) * (j ^ alpha))
    temp4 <- temp3 / (1 - exp(c * exp(eta)))
    sum = sum + temp4
  }
  W <- sum[, 1] * temp2 # W plays a role of weight matrix
  XWX <- t(X) %*% W %*% X
  c <- try(qr.solve(XWX), silent = TRUE)
  if (class(c) == "try-error")
    return("cannot calculated ; information matrix is near singular")
  else
    PreVar(location = useCond, formula = formula, nf = nf, infMtxInv = c)
}

#' Objective function of U optimal design with interval censoring
#'
#' The internal function to calculate the objective function value of
#' U optimal design with interval censoring plan.
#'
#' @keywords internal
Iobj.ic <- function(x, formula, coef, nf, t, k, alpha, useLower, useUpper) {
  X <- ExtendedForm(x, formula, nf)
  b <- coef
  dt <- t / k
  eta <- X %*% b # linear predictor
  temp1 <- exp(2 * eta)
  temp2 <- diag(temp1[, 1])
  sum <- 0
  for(j in 1:k) {
    c <- ((j - 1) ^ alpha - j ^ alpha) * dt ^ alpha
    temp3 <- c ^ 2 * exp(- exp(eta) * (dt ^ alpha) * (j ^ alpha))
    temp4 <- temp3 / (1 - exp(c * exp(eta)))
    sum = sum + temp4
  }
  W <- sum[, 1] * temp2 # W plays a role of weight matrix
  XWX <- t(X) %*% W %*% X
  c <- try(qr.solve(XWX), silent = TRUE)
  if (class(c) == "try-error")
    return("cannot calculated ; information matrix is near singular")
  else {
    # numerical integration
    intgratedPV <- cubature::adaptIntegrate(PreVar, lowerLimit = useLower,
                                            upperLimit = useUpper, formula = formula,
                                            nf = nf, infMtxInv = c)$integral
    volume <- 1
    for (i in 1:nf) volume <- volume * (useUpper[i] - useLower[i])
    intgratedPV / volume
  }
}

#' Design evaluation with interval censoring.
#'
#' \code{\link{alteval.ic}} calculates the objective function value
#'   (D, U or I) for a given design with interval censoring plan.
#'
#' @param designTable a data frame containing the coordinates and the number of
#'   allocation of each design point. The design created by either
#'   \code{\link{altopt.rc}} or \code{\link{altopt.ic}} or any design matrix
#'   with the same form as those can be provided for this argument.
#' @param optType the choice of \code{"D"}, \code{"U"} and \code{"I"} optimality.
#' @param t the total testing time.
#' @param k the number of time intervals.
#' @param nf the number of stress factors.
#' @param alpha the value of the shape parameter of Weibull distribution.
#' @param formula the object of class formula which is the linear predictor model.
#' @param coef the numeric vector containing the coefficients of each term in \code{formula}.
#' @param useCond the numeric vector of use condition.
#'   It should be provided when \code{optType} is \code{"U"}. The length of the vector
#'   should be same as the number of stress factors.
#' @param useLower the numeric vector of lower bound of use region.
#'   It should be provided when \code{optType} is \code{"I"}. The length of the vector
#'   should be same as the number of stress factors.
#' @param useUpper the numeric vector of upper bound of use region.
#'   It should be provided when \code{optType} is \code{"I"}. The length of the vector
#'   should be same as the number of stress factors.
#' @return The objective function value corresponded by \code{optType}
#'   for a given design with interval censoring plan.
#' @seealso \code{\link{altopt.ic}}
#' @examples
#' # Evaluation of factorial design for interval censoring.
#' x1 <- c(0, 1, 0, 1)
#' x2 <- c(0, 0, 1, 1)
#' allocation <- c(25, 25, 25, 25)
#' facDes <- data.frame(x1, x2, allocation)
#'
#' alteval.ic(facDes, "D", 30, 5, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01))
#'
#' alteval.ic(facDes, "U", 30, 5, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useCond = c(1.758, 3.159))
#'
#' alteval.ic(facDes, "I", 30, 5, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#' @export
alteval.ic <- function(designTable, optType, t, k, nf, alpha, formula, coef,
                       useCond, useLower, useUpper) {
  # Transform design to the single column array.
  x <- NULL
  for (col in (1:nf)) {
    for (row in (1:nrow(designTable))) {
      x <- c(x, rep(designTable[row, col], designTable[row, nf + 1]))
    }
  }

  if      (optType == "D")
    value <- Dobj.ic(x, formula, coef, nf, t, k, alpha)
  else if (optType == "U")
    value <- Uobj.ic(x, formula, coef, nf, t, k, alpha, useCond)
  else if (optType == "I")
    value <- Iobj.ic(x, formula, coef, nf, t, k, alpha, useLower, useUpper)
  else stop('Wrong optimization criteria')
  value
}

#' Optimal design with interval censoring.
#'
#' \code{\link{altopt.ic}} creates D, U or I optimal design
#' of the accelerated life testing with interval censoring plan.
#'
#' @param optType the choice of \code{"D"}, \code{"U"} and \code{"I"} optimality.
#' @param N the number of test units.
#' @param t the total testing time.
#' @param k the number of time intervals.
#' @param nf the number of stress factors.
#' @param alpha the value of the shape parameter of Weibull distribution.
#' @param formula the object of class formula which is the linear predictor model.
#' @param coef the numeric vector containing the coefficients of each term in \code{formula}.
#' @param useCond the numeric vector of use condition.
#'   It should be provided when \code{optType} is \code{"U"}. The length of the vector
#'   should be same as the number of stress factors.
#' @param useLower the numeric vector of lower bound of use region.
#'   It should be provided when \code{optType} is \code{"I"}. The length of the vector
#'   should be same as the number of stress factors.
#' @param useUpper the numeric vector of upper bound of use region.
#'   It should be provided when \code{optType} is \code{"I"}. The length of the vector
#'   should be same as the number of stress factors.
#' @param nOpt the number of repetition of optimization process. Default is 1.
#' @param nKM the number of repetition of k-means clustering. Default is 20.
#' @param nCls the number of clusters used for k-means clustering. If not specified,
#'   it is set as the number of parameters in the linear predictor model.
#' @return A list with components
#' \itemize{
#'   \item{call:}{ the matched call.}
#'   \item{opt.design.ori:}{ the original optimal design.}
#'   \item{opt.value.ori:}{ the objective function value of \code{opt.design.ori}.}
#'   \item{opt.design.rounded:}{ the optimal design clustered by rounding in third decimal points.}
#'   \item{opt.value.rounded:}{ the objective function value of \code{opt.design.rounded}.}
#'   \item{opt.design.kmeans:}{ the optimal design clustered by \code{\link[stats]{kmeans}}.}
#'   \item{opt.value.kmeans:}{ the objective function value of \code{opt.design.kmeans}.}
#' }
#' @references
#' {
#' Monroe, E. M., Pan, R., Anderson-Cook, C. M., Montgomery, D. C. and
#' Borror C. M. (2011) A Generalized Linear Model Approach to Designing
#' Accelerated Life Test Experiments, \emph{Quality and Reliability Engineering
#'  International} \bold{27(4)}, 595--607
#'
#' Yang, T., Pan, R. (2013) A Novel Approach to Optimal Accelerated Life Test
#'  Planning With Interval Censoring, \emph{Reliability, IEEE Transactions on}
#'   \bold{62(2)}, 527--536
#' }
#' @seealso \code{\link[stats]{kmeans}}, \code{\link{alteval.ic}}
#' @examples
#' \dontrun{
#' # Generating D optimal design for interval censoring.
#' altopt.ic("D", 100, 30, 5, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01))
#'
#' # Generating U optimal design for interval censoring.
#' altopt.ic("D", 100, 30, 5, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useCond = c(1.758, 3.159))
#'
#' # Generating I optimal design for interval censoring.
#' altopt.ic("D", 100, 30, 5, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useLower = c(1.458, 2.859),
#' useUpper = c(2.058, 3.459))
#' }
#' @importFrom stats aggregate
#' @importFrom stats optim
#' @importFrom stats runif
#' @export
altopt.ic <- function(optType, N, t, k, nf, alpha, formula, coef,
                      useCond, useLower, useUpper,
                      nOpt = 1, nKM = 30, nCls = NULL) {

  # function for optimization
  Opt.ic <- function() {
    xInit <- runif(nf * N, min = 0, max = 1)
    lb <- rep(0, nf * N)
    ub <- rep(1, nf * N)
    if      (optType == "D")
      solution <- optim(xInit, Dobj.ic, NULL, formula, coef, nf, t, k, alpha,
                        method = "L-BFGS-B", lower = lb, upper = ub,
                        control = list(fnscale = -1) # Maximization
      )
    else if (optType == "U")
      solution <- optim(xInit, Uobj.ic, NULL, formula, coef, nf, t, k, alpha,
                        useCond,
                        method = "L-BFGS-B", lower = lb, upper = ub)
    else if (optType == "I")
      solution <- optim(xInit, Iobj.ic, NULL, formula, coef, nf, t, k, alpha,
                        useLower, useUpper,
                        method = "L-BFGS-B", lower = lb, upper = ub)
    else stop('Wrong optimization criteria')
    solution
  }

  # Repeat optimization with different initial points
  for (i in 1:nOpt) {
    Curr_sol <- Opt.ic()
    if (i == 1) Best_sol <- Curr_sol
    else if (optType == "D" && Curr_sol$value > Best_sol$value)
      Best_sol <- Curr_sol
    else if ((optType %in% c("U", "I")) && Curr_sol$value < Best_sol$value)
      Best_sol <- Curr_sol
  }

  terms <- attr(terms(formula), "term.labels")
  optDesignMtx <- as.data.frame(matrix(Best_sol$par, ncol = nf))
  colnames(optDesignMtx) <- terms[1:nf]

  # Original result

  columnlist <- apply(optDesignMtx, 2, list)
  columnlist <- lapply(columnlist, unlist)
  optDesignOri <- aggregate(optDesignMtx, columnlist, length)
  while (length(optDesignOri) != (nf + 1))
    optDesignOri <- optDesignOri[, -length(optDesignOri)]
  names(optDesignOri)[ncol(optDesignOri)] <- paste("allocation")
  optValueOri <- alteval.ic(optDesignOri, optType, t, k, nf, alpha,
                              formula, coef, useCond, useLower, useUpper)

  # Rounding result
  optDesignMtxRound <- round(optDesignMtx, digits = 3)
  columnlist <- apply(optDesignMtxRound, 2, list)
  columnlist <- lapply(columnlist, unlist)
  optDesignRound <- aggregate(optDesignMtxRound, columnlist, length)
  while (length(optDesignRound) != (nf + 1))
    optDesignRound <- optDesignRound[, -length(optDesignRound)]
  names(optDesignRound)[ncol(optDesignRound)] <- paste("allocation")
  optValueRound <- alteval.ic(optDesignRound, optType, t, k, nf, alpha,
                              formula, coef, useCond, useLower, useUpper)


  # k-means result
  if (is.null(nCls)) nCls <- length(coef)
  for (j in 1:nKM) {
    curDesignKmeans <- kmeansCls(optDesignMtx, nCls)
    curValueKmeans <- alteval.ic(curDesignKmeans, optType, t, k, nf, alpha,
                                 formula, coef, useCond, useLower, useUpper)
    if (j == 1) {
      optDesignKmeans <- curDesignKmeans
      optValueKmeans <- curValueKmeans
    } else if (optType == "D" && curValueKmeans > optValueKmeans) {
      optDesignKmeans <- curDesignKmeans
      optValueKmeans <- curValueKmeans
    } else if ((optType %in% c("U", "I")) && curValueKmeans < optValueKmeans) {
      optDesignKmeans <- curDesignKmeans
      optValueKmeans <- curValueKmeans
    }
  }

  # Creates output
  out <- list(call = match.call(),
              opt.design.ori = optDesignOri,
              opt.value.ori = as.numeric(optValueOri),
              opt.design.rounded = optDesignRound,
              opt.value.rounded = as.numeric(optValueRound),
              opt.design.kmeans = optDesignKmeans,
              opt.value.kmeans = ifelse (class(optValueKmeans) == "character",
                                     optValueKmeans, as.numeric(optValueKmeans))
  )
  out
}

#' Design plot.
#'
#' \code{\link{design.plot}} draws design plot as a form of a bubble plot
#' of any two stress factors which are specified by \code{xAxis} and \code{yAxis}.
#' The size of each bubble indicates the relative magnitude of allocation on
#' each design point.
#'
#' @param design the data frame containing the coordinates and the number of
#'   allocation of each design point. The design created by either
#'   \code{\link{altopt.rc}} or \code{\link{altopt.ic}} or any design matrix
#'   with the same form as those can be provided for this argument.
#' @param xAxis the name of the factor to be displayed in x axis.
#' @param yAxis the name of the factor to be displayed in y axis.
#' @return The bubble plot of a design with two stress factors.
#' @examples
#' \dontrun{
#' # Design plot of D optimal design with right censoring.
#' Design1 <- altopt.rc("D", 100, 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01))
#'
#' design.plot(Design1$opt.design.rounded, x1, x2)
#' }
#' @importFrom graphics plot
#' @importFrom graphics rect
#' @importFrom graphics symbols
#' @importFrom graphics text
#' @importFrom stats aggregate
#' @export
design.plot <- function (design, xAxis, yAxis) {
  designName <- deparse(substitute(design))
  xAxisName <- deparse(substitute(xAxis))
  yAxisName <- deparse(substitute(yAxis))

  plot(1, type = "n", main = designName, xlab = xAxisName, ylab = yAxisName,
       xaxp = c(0, 1, 5), yaxp = c(0, 1, 5),
       xlim = c(-0.2, 1.2), ylim = c(-0.2, 1.2), frame = FALSE)
  rect(0, 0, 1, 1)
  agg.des <- aggregate(design, by = list(design[, colnames(design) == xAxisName],
                                         design[, colnames(design) == yAxisName]),
                       FUN = sum)
  symbols(agg.des$Group.1, agg.des$Group.2, circles = agg.des$allocation / 300,
          inches = FALSE, add = TRUE, fg = "blue", bg = "white", lwd = 1.5)
  text(agg.des$Group.1, agg.des$Group.2, agg.des$allocation, cex = .75)
}

#' Contour plot of prediction variance for a design with right censoring.
#'
#' \code{\link{pv.contour.rc}} draws the contour plot of prediction variance
#' for a given design with right censoring plan. Either \code{useCond} or
#' use region (\code{useLower} and \code{useUpper}) should be
#' provided.
#'
#' @param design the data frame containing the coordinates and the number of
#'   allocation of each design point. The design created by either
#'   \code{\link{altopt.rc}} or \code{\link{altopt.ic}} or any design matrix
#'   with the same form as those can be provided for this argument.
#' @param xAxis the name of the factor to be displayed in x axis.
#' @param yAxis the name of the factor to be displayed in y axis.
#' @param tc the censoring time.
#' @param nf the number of stress factors.
#' @param alpha the value of the shape parameter of Weibull distribution.
#' @param formula the object of class formula which is the linear predictor model.
#' @param coef the numeric vector containing the coefficients of each term in \code{formula}.
#' @param useCond the vector of specified use condition. If it is provided,
#'   the contour line will be generated up to this point.
#' @param useLower,useUpper the vector of the use region. If these are
#'   provided, the contour line will be generated up to this region.
#'   Note that either \code{useCond} or both of \code{useLower, useUpper}
#'   should be provided.
#' @return The contour plot of prediction variance for right censoring.
#' @seealso \code{\link{altopt.rc}}
#' @examples
#' \dontrun{
#' # Contour plot of prediction variance of U optimal design with right censoring.
#' Design <- altopt.rc("D", 100, 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useCond = c(1.758, 3.159))
#'
#' pv.contour.rc(Design$opt.design.rounded, x1, x2, 100, 2, 1,
#' formula = ~ x1 + x2 + x1:x2, coef = c(0, -4.086, -1.476, 0.01), useCond = c(1.758, 3.159))
#' }
#' @importFrom graphics contour
#' @importFrom graphics points
#' @importFrom graphics rect
#' @importFrom graphics segments
#' @importFrom graphics symbols
#' @importFrom graphics text
#' @importFrom graphics title
#' @importFrom stats aggregate
#' @export
pv.contour.rc <- function (design, xAxis, yAxis, tc, nf, alpha, formula, coef,
                           useCond = NULL, useLower = NULL, useUpper = NULL) {
  designName <- deparse(substitute(design))
  xAxisName <- deparse(substitute(xAxis))
  yAxisName <- deparse(substitute(yAxis))
  xColNum <- which(colnames(design) == xAxisName)
  yColNum <- which(colnames(design) == yAxisName)

  pv <- function (useCond, ux, uy) {
    useCond[xColNum] <- ux # change use condition of two factors given by user
    useCond[yColNum] <- uy # with maintaining values of the other factors
    as.numeric(alteval.rc(design, optType = "U", tc, nf,
                          alpha, formula, coef, useCond))
  }

  if (!is.null(useCond)) {
    range <- ceiling(max(useCond[xColNum], useCond[yColNum]))
    pv.use <- round(pv(useCond, useCond[xColNum], useCond[yColNum]), digits = 2)
  }
  else if (!is.null(useLower) && !is.null(useUpper)) {
    range <- ceiling(max(useUpper[xColNum], useUpper[yColNum]))
    pv.use <- round(pv(useUpper, useUpper[xColNum], useUpper[yColNum]), digits = 2)
  }
  else stop('Use condition missing, either useCond or
            both of useLower and useUpper should be provided.')

  x <- seq(0, range, 0.1)
  y <- seq(0, range, 0.1)

  pv.grid <- matrix(nrow = length(x), ncol = length(y))
  for (r in 1:length(x)) {
    for (c in 1:length(y)) {
      if (!is.null(useCond))
        pv.grid[r, c] <- pv(useCond, x[r], y[c])
      else if (!is.null(useLower) && !is.null(useUpper))
        pv.grid[r, c] <- pv((useLower + useUpper) / 2, x[r], y[c])
    }
  }

  mylevels <- seq(0, pv.use, 0.1)

  contour(x, y, pv.grid, levels = mylevels, method = "edge", pty = "s")
  title(main = paste("PV contour of ", designName, sep = ""),
        xlab = xAxisName, ylab = yAxisName, cex.main = .75)
  agg.des <- aggregate(design, by = list(design[, colnames(design) == xAxisName],
                                         design[, colnames(design) == yAxisName]),
                       FUN = sum)
  symbols(agg.des$Group.1, agg.des$Group.2, circles = agg.des$allocation / 300,
          inches = FALSE, add = TRUE, bg = "gray")
  text(agg.des$Group.1, agg.des$Group.2, agg.des$allocation, cex = .75,
       adj = c(-.5, 1))

  segments(0, 0, range, 0)
  segments(0, 0, 0, range)
  segments(range, 0, range, range)
  segments(0, range, range, range)
  segments(0, 1, 1, 1)
  segments(1, 0, 1, 1)
  if (!is.null(useCond)) {
    points(useCond[xColNum], useCond[yColNum], pch=22, bg="white")
    text(useCond[xColNum], useCond[yColNum], paste("PV =", pv.use),
         cex = .75, adj = c(-.2, -.2))
    segments(0, useCond[yColNum], useCond[xColNum], useCond[yColNum])
    segments(useCond[xColNum], 0, useCond[xColNum], useCond[yColNum])
  }
  else if (!is.null(useLower) && !is.null(useUpper)) {
    rect(useLower[xColNum], useLower[yColNum],
         useUpper[xColNum], useUpper[yColNum], col="white")
    text((useLower[xColNum] + useUpper[xColNum]) / 2,
         (useLower[yColNum] + useUpper[yColNum]) / 2,
         "Use Region", cex = .5, adj = c(0.5, 0.5))
  }
}

#' Contour plot of prediction variance for a design with interval censoring.
#'
#' \code{\link{pv.contour.ic}} draws the contour plot of prediction variance
#' for a given design with interval censoring plan. Either \code{useCond} or
#' use region (\code{useLower} and \code{useUpper}) should be
#' provided.
#'
#' @param design the data frame containing the coordinates and the number of
#'   allocation of each design point. The design created by either
#'   \code{\link{altopt.rc}} or \code{\link{altopt.ic}} or any design matrix
#'   with the same form as those can be provided for this argument.
#' @param xAxis the name of the factor to be displayed in x axis.
#' @param yAxis the name of the factor to be displayed in y axis.
#' @param t the total testing time.
#' @param k the number of time intervals.
#' @param nf the number of stress factors.
#' @param alpha the value of the shape parameter of Weibull distribution.
#' @param formula the object of class formula which is the linear predictor model.
#' @param coef the numeric vector containing the coefficients of each term in \code{formula}.
#' @param useCond the vector of specified use condition. If it is provided,
#'   the contour line will be generated up to this point.
#' @param useLower,useUpper the vector of the use region. If these are
#'   provided, the contour line will be generated up to this region.
#'   Note that either \code{useCond} or both of \code{useLower, useUpper}
#'   should be provided.
#' @return The contour plot of prediction variance for interval censoring.
#' @seealso \code{\link{altopt.ic}}
#' @examples
#' \dontrun{
#' # Contour plot of prediction variance of U optimal design with interval censoring.
#' Design <- altopt.ic("D", 100, 30, 5, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useCond = c(1.758, 3.159))
#'
#' pv.contour.ic(Design$opt.design.rounded, x1, x2, 30, 5, 2, 1,
#' formula = ~ x1 + x2 + x1:x2, coef = c(0, -4.086, -1.476, 0.01), useCond = c(1.758, 3.159))
#' }
#' @importFrom graphics contour
#' @importFrom graphics points
#' @importFrom graphics rect
#' @importFrom graphics segments
#' @importFrom graphics symbols
#' @importFrom graphics text
#' @importFrom graphics title
#' @importFrom stats aggregate
#' @export
pv.contour.ic <- function (design, xAxis, yAxis, t, k, nf, alpha, formula, coef,
                           useCond = NULL, useLower = NULL, useUpper = NULL) {
  designName <- deparse(substitute(design))
  xAxisName <- deparse(substitute(xAxis))
  yAxisName <- deparse(substitute(yAxis))
  xColNum <- which(colnames(design) == xAxisName)
  yColNum <- which(colnames(design) == yAxisName)

  pv <- function (useCond, ux, uy) {
    useCond[xColNum] <- ux # change use condition of two factors given by user
    useCond[yColNum] <- uy # with maintaining values of the other factors
    as.numeric(alteval.ic(design, optType = "U", t, k, nf,
                          alpha, formula, coef, useCond))
  }

  if (!is.null(useCond)) {
    range <- ceiling(max(useCond[xColNum], useCond[yColNum]))
    pv.use <- round(pv(useCond, useCond[xColNum], useCond[yColNum]), digits = 2)
  }
  else if (!is.null(useLower) && !is.null(useUpper)) {
    range <- ceiling(max(useUpper[xColNum], useUpper[yColNum]))
    pv.use <- round(pv(useUpper, useUpper[xColNum], useUpper[yColNum]), digits = 2)
  }
  else stop('Use condition missing, either useCond or
            both of useLower and useUpper should be provided.')

  x <- seq(0, range, 0.1)
  y <- seq(0, range, 0.1)

  pv.grid <- matrix(nrow = length(x), ncol = length(y))
  for (r in 1:length(x)) {
    for (c in 1:length(y)) {
      if (!is.null(useCond))
        pv.grid[r, c] <- pv(useCond, x[r], y[c])
      else if (!is.null(useLower) && !is.null(useUpper))
        pv.grid[r, c] <- pv((useLower + useUpper) / 2, x[r], y[c])
    }
  }

  mylevels <- seq(0, pv.use, 0.1)

  contour(x, y, pv.grid, levels = mylevels, method = "edge", pty="s")
  title(main = paste("PV contour of ", designName, sep = ""),
        xlab = xAxisName, ylab = yAxisName, cex.main = .75)
  agg.des <- aggregate(design, by = list(design[, colnames(design) == xAxisName],
                                         design[, colnames(design) == yAxisName]),
                       FUN = sum)
  symbols(agg.des$Group.1, agg.des$Group.2, circles = agg.des$allocation / 300,
          inches = FALSE, add = TRUE, bg = "gray")
  text(agg.des$Group.1, agg.des$Group.2, agg.des$allocation, cex = .75,
       adj = c(-.5, 1))

  segments(0, 0, range, 0)
  segments(0, 0, 0, range)
  segments(range, 0, range, range)
  segments(0, range, range, range)
  segments(0, 1, 1, 1)
  segments(1, 0, 1, 1)
  if (!is.null(useCond)) {
    points(useCond[xColNum], useCond[yColNum], pch=22, bg="white")
    text(useCond[xColNum], useCond[yColNum], paste("PV =", pv.use),
         cex = .75, adj = c(-.2, -.2))
    segments(0, useCond[yColNum], useCond[xColNum], useCond[yColNum])
    segments(useCond[xColNum], 0, useCond[xColNum], useCond[yColNum])
  }
  else if (!is.null(useLower) && !is.null(useUpper)) {
    rect(useLower[xColNum], useLower[yColNum],
         useUpper[xColNum], useUpper[yColNum], col="white")
    text((useLower[xColNum] + useUpper[xColNum]) / 2,
         (useLower[yColNum] + useUpper[yColNum]) / 2,
         "Use Region", cex = .5, adj = c(0.5, 0.5))
  }
}

#' FUS (Fraction of Use Space) plot for right censoring.
#'
#' \code{\link{pv.fus.rc}} draws the FUS plot of prediction variance
#' for a given design with right censoring plan. The use region
#' (\code{useLower} and \code{useUpper}) should be
#' provided.
#'
#' @param useLower,useUpper the vectors containing the lower bound and upper
#'   bound for the use region. They should be provided for FUS plot.
#' @return The "trellis" object which includes the FUS plot
#'   for right censoring.
#' @inheritParams design.plot
#' @inheritParams altopt.rc
#' @seealso \code{\link{altopt.rc}}
#' @examples
#' \dontrun{
#' # FUS plot of I optimal design with right censoring.
#' Design <- altopt.rc("I", 100, 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#'
#' pv.fus.rc(Design$opt.design.rounded, 100, 2, 1,
#' formula = ~ x1 + x2 + x1:x2, coef = c(0, -4.086, -1.476, 0.01),
#' useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#' }
#' @export
pv.fus.rc <- function (design, tc, nf, alpha, formula, coef,
                       useLower = NULL, useUpper = NULL) {
  if (is.null(useLower) || is.null(useUpper)) stop('Use condition missing')

  d <- round(5000 ^ (1 / nf))
  pv.grid <- matrix(nrow = d ^ nf, ncol = nf + 1)
  for (p in 1:nf) {
    range <- seq(useLower[p], useUpper[p], length.out = d)
    pv.grid[, p] <- rep(rep(range, each = d ^ (nf - p)), d ^ (p - 1))
  }

  for (r in 1:nrow(pv.grid)) {
    pv.grid[r, ncol(pv.grid)] <- as.numeric(alteval.rc(design, optType = "U",
                                                       tc, nf, alpha, formula, coef, useCond = pv.grid[r, 1:nf]))
  }

  fus <- cbind(sort(pv.grid[, ncol(pv.grid)]), c(1:nrow(pv.grid)) / nrow(pv.grid))
  plot <- lattice::xyplot(fus[, 1] ~ fus[, 2],
                          aspect = 1 / 2,
                          main = paste("FUS of ", deparse(substitute(design)), sep = ""),
                          xlab = paste("Fraction of Use Space"),
                          ylab = paste("Prediction Variance"),
                          type = "a",
                          grid = TRUE,
                          scales = list(x = list(tick.number = 11)))
  plot
}

#' FUS (Fraction of Use Space) plot for interval censoring.
#'
#' \code{\link{pv.fus.ic}} draws the FUS plot of prediction variance
#' for a given design with interval censoring plan. The use region
#' (\code{useLower} and \code{useUpper}) should be
#' provided.
#'
#' @param useLower,useUpper the vectors containing the lower bound and upper
#'   bound for the use region. They should be provided for FUS plot.
#' @return The "trellis" object which includes the FUS plot
#'   for interval censoring.
#' @inheritParams design.plot
#' @inheritParams altopt.ic
#' @seealso \code{\link{altopt.ic}}
#' @examples
#' \dontrun{
#' # FUS plot of I optimal design with interval censoring.
#' Design <- altopt.ic("I", 100, 30, 5, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#'
#' pv.fus.ic(Design$opt.design.rounded, 30, 5, 2, 1,
#' formula = ~ x1 + x2 + x1:x2, coef = c(0, -4.086, -1.476, 0.01),
#' useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#' }
#' @export
pv.fus.ic <- function (design, t, k, nf, alpha, formula, coef,
                       useLower = NULL, useUpper = NULL) {
  if (is.null(useLower) || is.null(useUpper)) stop('Use condition missing')

  d <- round(5000 ^ (1 / nf))
  pv.grid <- matrix(nrow = d ^ nf, ncol = nf + 1)
  for (p in 1:nf) {
    range <- seq(useLower[p], useUpper[p], length.out = d)
    pv.grid[, p] <- rep(rep(range, each = d ^ (nf - p)), d ^ (p - 1))
  }

  for (r in 1:nrow(pv.grid)) {
    pv.grid[r, ncol(pv.grid)] <- as.numeric(alteval.ic(design, optType="U",
                                                       t, k, nf, alpha, formula, coef, useCond = pv.grid[r, 1:nf]))
  }

  fus <- cbind(sort(pv.grid[, ncol(pv.grid)]), c(1:nrow(pv.grid)) / nrow(pv.grid))
  plot <- lattice::xyplot(fus[, 1] ~ fus[, 2],
                          aspect = 1 / 2,
                          main = paste("FUS of ", deparse(substitute(design)), sep = ""),
                          xlab = paste("Fraction of Use Space"),
                          ylab = paste("Prediction Variance"),
                          type = "a",
                          grid = TRUE,
                          scales = list(x = list(tick.number = 11)))
  plot
}

#' Comparing designs using FUS
#'
#' \code{\link{compare.fus}} draws the FUS plots of multiple designs on a
#' single frame.
#'
#' @param ... Objects created by \code{\link{pv.fus.rc}} or
#' \code{\link{pv.fus.ic}}.
#' @return FUS plots of multiple designs.
#' @seealso \code{\link{pv.fus.rc}}, \code{\link{pv.fus.ic}}
#' @examples
#' \dontrun{
#' # Generating D optimal design and FUS plot.
#' Dopt <- altopt.rc("D", 100, 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01))
#'
#' FUS.D <- pv.fus.rc(Dopt$opt.design.rounded, 100, 2, 1,
#' formula = ~ x1 + x2 + x1:x2, coef = c(0, -4.086, -1.476, 0.01),
#' useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#'
#' # Generating U optimal design and FUS plot.
#' Uopt <- altopt.rc("U", 100, 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useCond = c(1.758, 3.159))
#'
#' FUS.U <- pv.fus.rc(Uopt$opt.design.rounded, 100, 2, 1,
#' formula = ~ x1 + x2 + x1:x2, coef = c(0, -4.086, -1.476, 0.01),
#' useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#'
#' # Comparing D and U optimal designs.
#' compare.fus(FUS.D, FUS.U)
#' }
#' @export
compare.fus <- function (...) {
  # Keep the name of arguments as string
  nm <- unlist(strsplit(deparse(substitute(list(...))), split = ","))
  nm <- unlist(strsplit(nm, split = " "))
  nm <- unlist(strsplit(nm, split = "list(", fixed = TRUE))
  nm <- unlist(strsplit(nm, split = ")", fixed = TRUE))

  fus.obj <- list(...)

  data <- data.frame()
  for (i in 1:length(fus.obj)) {
    data <- rbind(data, data.frame(y = fus.obj[[i]]$panel.args[[1]]$y,
                                   x = fus.obj[[i]]$panel.args[[1]]$x,
                                   z = nm[i]))
  }

  plot <- lattice::xyplot(y ~ x, data, groups = data$z, auto.key = list(corner = c(0, 1)),
                          aspect = 1 / 2,
                          main = paste("FUS of ", deparse(substitute(list(...))), sep = ""),
                          xlab = paste("Fraction of Use Space"),
                          ylab = paste("Prediction Variance"),
                          type = "a",
                          grid = TRUE,
                          scales = list(x = list(tick.number = 11)))
  plot
}

#' VDUS (Variance Dispersion of Use Space) plot for right censoring.
#'
#' \code{\link{pv.vdus.rc}} draws the VDUS plot of prediction variance
#' for a given design with right censoring plan. The use region
#' (\code{useLower} and \code{useUpper}) should be
#' provided.
#'
#' @param useLower,useUpper the vectors containing the lower bound and upper
#'   bound for the use region. They should be provided for VDUS plot.
#' @return The "trellis" object which includes the VDUS plot
#'   for right censoring.
#' @inheritParams design.plot
#' @inheritParams altopt.rc
#' @seealso \code{\link{altopt.rc}}
#' @examples
#' \dontrun{
#' # VDUS plot of I optimal design with right censoring.
#' Design <- altopt.rc("I", 100, 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#'
#' pv.vdus.rc(Design$opt.design.rounded, 100, 2, 1,
#' formula = ~ x1 + x2 + x1:x2, coef = c(0, -4.086, -1.476, 0.01),
#' useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#' }
#' @importFrom stats aggregate
#' @export
pv.vdus.rc <- function (design, tc, nf, alpha, formula, coef,
                        useLower = NULL, useUpper = NULL) {
  if (is.null(useLower) || is.null(useUpper)) stop('Use condition missing')

  d <- round(5000 ^ (1 / nf))
  if (d %% 2 == 0) d <- d + 1
  pv.grid <- matrix(nrow = d ^ nf, ncol = nf)
  for (p in 1:nf) {
    range <- seq(useLower[p], useUpper[p], length.out = d)
    pv.grid[, p] <- rep(rep(range, each = d ^ (nf - p)), d ^ (p - 1))
  }

  center <- (useLower + useUpper) / 2
  dx <- (useUpper - center) / ((d - 1) / 2)
  rank <- numeric(length = nrow(pv.grid))
  for (r in 1:nrow(pv.grid)) {
    if (all(pv.grid[r, ] == center)) rank[r] <- 0
    else {
      for (n in 1:(d - 1) / 2) {
        for (p in 1:nf) {
          if ((round(pv.grid[r, p], 3) == round(center[p] + n * dx[p], 3)
               || round(pv.grid[r, p], 3) == round(center[p] - n * dx[p], 3))
              && round(pv.grid[r, -p], 3) >= round(center[-p] - n * dx[-p], 3)
              && round(pv.grid[r, -p], 3) <= round(center[-p] + n * dx[-p], 3))
            rank[r] <- n
        }
      }
    }
  }
  pv.grid <- cbind(pv.grid, rank)
  pv <- numeric(length = nrow(pv.grid))
  for (r in 1:nrow(pv.grid)) {
    pv[r] <- as.numeric(alteval.rc(design, optType = "U", tc, nf, alpha,
                                   formula, coef, useCond = pv.grid[r, 1:nf]))
  }
  pv.grid <- cbind(pv.grid, pv)

  min <- aggregate(pv.grid, by = list(pv.grid[, "rank"]), FUN = min)
  avg <- aggregate(pv.grid, by = list(pv.grid[, "rank"]), FUN = mean)
  max <- aggregate(pv.grid, by = list(pv.grid[, "rank"]), FUN = max)
  vdus <- rbind(cbind(min, gr = "min"), cbind(avg, gr = "avg"),
                cbind(max, gr = "max"))

  plot <- lattice::xyplot(vdus$pv ~ vdus$rank,
                          aspect = 1 / 2,
                          main = paste("VDUS of ", deparse(substitute(design)), sep = ""),
                          xlab = paste("Relative radius from origin"),
                          ylab = paste("Prediction Variance"),
                          type = "smooth",
                          grid = TRUE,
                          group = vdus$gr,
                          auto.key = list(corner = c(0, 1)))
  plot
}

#' VDUS (Variance Dispersion of Use Space) plot for interval censoring.
#'
#' \code{\link{pv.vdus.ic}} draws the VDUS plot of prediction variance
#' for a given design with interval censoring plan. The use region
#' (\code{useLower} and \code{useUpper}) should be
#' provided.
#'
#' @param useLower,useUpper the vectors containing the lower bound and upper
#'   bound for the use region. They should be provided for VDUS plot.
#' @return The "trellis" object which includes the VDUS plot
#'   for interval censoring.
#' @inheritParams design.plot
#' @inheritParams altopt.ic
#' @seealso \code{\link{altopt.ic}}
#' @examples
#' \dontrun{
#' # VDUS plot of I optimal design with interval censoring.
#' Design <- altopt.ic("I", 100, 30, 5, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#'
#' pv.vdus.ic(Design$opt.design.rounded, 30, 5, 2, 1,
#' formula = ~ x1 + x2 + x1:x2, coef = c(0, -4.086, -1.476, 0.01),
#' useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#' }
#' @importFrom stats aggregate
#' @export
pv.vdus.ic <- function (design, t, k, nf, alpha, formula, coef,
                        useLower = NULL, useUpper = NULL) {
  if (is.null(useLower) || is.null(useUpper)) stop('Use condition missing')

  d <- round(5000 ^ (1 / nf))
  if (d %% 2 == 0) d <- d + 1
  pv.grid <- matrix(nrow = d ^ nf, ncol = nf)
  for (p in 1:nf) {
    range <- seq(useLower[p], useUpper[p], length.out = d)
    pv.grid[, p] <- rep(rep(range, each = d ^ (nf - p)), d ^ (p - 1))
  }

  center <- (useLower + useUpper) / 2
  dx <- (useUpper - center)/((d - 1) / 2)
  rank <- numeric(length = nrow(pv.grid))
  for (r in 1:nrow(pv.grid)) {
    if (all(pv.grid[r, ] == center)) rank[r] <- 0
    else {
      for (n in 1:(d - 1) / 2) {
        for (p in 1:nf) {
          if ((round(pv.grid[r, p], 3) == round(center[p] + n * dx[p], 3)
               || round(pv.grid[r, p], 3) == round(center[p] - n * dx[p], 3))
              && round(pv.grid[r, -p], 3) >= round(center[-p] - n * dx[-p], 3)
              && round(pv.grid[r, -p], 3) <= round(center[-p] + n * dx[-p], 3))
            rank[r] <- n
        }
      }
    }
  }
  pv.grid <- cbind(pv.grid, rank)
  pv <- numeric(length = nrow(pv.grid))
  for (r in 1:nrow(pv.grid)) {
    pv[r] <- as.numeric(alteval.ic(design, optType = "U", t, k, nf, alpha,
                                   formula, coef, useCond = pv.grid[r, 1:nf]))
  }
  pv.grid <- cbind(pv.grid, pv)

  min <- aggregate(pv.grid, by = list(pv.grid[, "rank"]), FUN = min)
  avg <- aggregate(pv.grid, by = list(pv.grid[, "rank"]), FUN = mean)
  max <- aggregate(pv.grid, by = list(pv.grid[, "rank"]), FUN = max)
  vdus <- rbind(cbind(min, gr = "min"), cbind(avg, gr = "avg"),
                cbind(max, gr = "max"))

  plot <- lattice::xyplot(vdus$pv ~ vdus$rank,
                          aspect = 1 / 2,
                          main = paste("VDUS of ", deparse(substitute(design)), sep = ""),
                          xlab = paste("Relative radius from origin"),
                          ylab = paste("Prediction Variance"),
                          type = "smooth",
                          grid = TRUE,
                          group = vdus$gr,
                          auto.key = list(corner = c(0, 1)))
  plot
}

#' Comparing designs using VDUS
#'
#' \code{\link{compare.vdus}} draws the VDUS plots of multiple designs on a
#' single frame.
#'
#' @param ... Objects created by \code{\link{pv.vdus.rc}} or
#' \code{\link{pv.vdus.ic}}.
#' @return VDUS plots of multiple designs.
#' @seealso \code{\link{pv.vdus.rc}}, \code{\link{pv.vdus.ic}}
#' @examples
#' \dontrun{
#' # Generating D optimal design and VDUS plot.
#' Dopt <- altopt.rc("D", 100, 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01))
#'
#' VDUS.D <- pv.vdus.rc(Dopt$opt.design.rounded, 100, 2, 1,
#' formula = ~ x1 + x2 + x1:x2, coef = c(0, -4.086, -1.476, 0.01),
#' useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#'
#' # Generating U optimal design and VDUS plot.
#' Uopt <- altopt.rc("U", 100, 100, 2, 1, formula = ~ x1 + x2 + x1:x2,
#' coef = c(0, -4.086, -1.476, 0.01), useCond = c(1.758, 3.159))
#'
#' VDUS.U <- pv.vdus.rc(Uopt$opt.design.rounded, 100, 2, 1,
#' formula = ~ x1 + x2 + x1:x2, coef = c(0, -4.086, -1.476, 0.01),
#' useLower = c(1.458, 2.859), useUpper = c(2.058, 3.459))
#'
#' # Comparing D and U optimal designs.
#' compare.vdus(VDUS.D, VDUS.U)
#' }
#' @export
compare.vdus <- function (...) {
  nm <- unlist(strsplit(deparse(substitute(list(...))), split = ","))
  nm <- unlist(strsplit(nm, split = " "))
  nm <- unlist(strsplit(nm, split = "list(", fixed = TRUE))
  nm <- unlist(strsplit(nm, split = ")", fixed = TRUE))

  vdus.obj <- list(...)

  data <- data.frame()
  for (i in 1:length(vdus.obj)) {
    data <- rbind(data, data.frame(y = vdus.obj[[i]]$panel.args[[1]]$y,
                                   x = vdus.obj[[i]]$panel.args[[1]]$x,
                                   z = paste(nm[i], vdus.obj[[i]]$panel.args.common$groups),
                                   sep = "."))
  }

  plot <- lattice::xyplot(y ~ x, data, groups = data$z, auto.key = list(corner = c(0, 1)),
                          aspect = 1 / 2,
                          main = paste("VDUS of ", deparse(substitute(list(...))), sep = ""),
                          xlab = paste("Variance Dispersion of Use Space"),
                          ylab = paste("Prediction Variance"),
                          type = "a",
                          grid = TRUE,
                          scales = list(x = list(tick.number = 11)))
  plot
}

#' Coding and decoding stress level
#'
#' Convert the stress levels from the actual levels to standardized levels,
#'   and vice versa.
#'
#' @param lowStLv a numeric vector containing the actual lowest stress level
#'   of each stress variable in design region.
#' @param highStLv a numeric vector containing the actual highest stress level
#'   of each stress variable in design region.
#' @param actual a data frame or numeric vector containing the design points
#'   in actual units.
#' @param stand a data frame or numeric vector containing the design points
#'   in standardized units.
#' @return When \code{actual} is provided, the function converts it to the
#'   standardized units and when \code{stand} is provided, the function converts
#'   it to the actual units.
#' @examples
#' \dontrun{
#'   # Generating D optimal design in coded unit.
#'   Design <- altopt.rc(optType = "D", N = 100, tc = 100, nf = 2, alpha = 1,
#'   formula = ~x1 + x2 + x1:x2, coef = c(0, -4.086, -1.476, 0.01))
#'
#'   # Transform the coded unit to actual stress variable's level.
#'   convert.stress.level(lowStLv = c(34.834, 4.094), highStLv = c(30.288, 4.5),
#'   stand = Design$opt.design.rounded)
#'
#'   # Transform the actual stress level to coded units.
#'   use <- c(38.281, 3.219)
#'   convert.stress.level(lowStLv = c(34.834, 4.094), highStLv = c(30.288, 4.5),
#'   actual = use)
#'   }
#' @export
convert.stress.level <- function(lowStLv, highStLv,
                                 actual = NULL, stand = NULL) {
  nf <- length(lowStLv)
  if (is.null(actual) && is.null(stand))
    stop ('Either actual or Stand should be provided.')
  else if (!is.null(actual) && !is.null(stand))
    stop ('Only one of actual or Stand should be provided')
  else if (!is.null(stand)) {
    # Convert from stand to actual
    if (class(stand) == "numeric")
      stand <- as.data.frame(matrix(stand, ncol = nf))
    out <- stand
    for (c in 1:nf) {
      for (r in 1:nrow(stand))
        out[r, c] <- stand[r, c] * (lowStLv[c] - highStLv[c]) + highStLv[c]
    }
  }
  else if (!is.null(actual)) {
    # Convert from actual to stand
    if (class(actual) == "numeric")
      actual <- as.data.frame(matrix(actual, ncol = nf))
    out <- actual
    for (c in 1:nf) {
      for (r in 1:nrow(actual))
        out[r, c] <- (actual[r, c] - highStLv[c]) / (lowStLv[c] - highStLv[c])
    }
  }
  out
}
