#' Performing multivariate dose-response meta-analysis
#' 
#' @description Estimates a dose-response relation from either a single or multiple summarized data, taking into account the correlation among set of log relative risks. 
#' The covariances are approximated according to two different methods, proposed respectively by Greeland S., Longnecker M., and Hamling J.; alternatively the user can provide directly the covariance matrices or the avarage covariances (Easton D.).
#' The study-specific estimates are combined according to principles of multivariate random-effects meta-analysis.
#' 
#' @param formula an object of class "\code{\link{formula}}" offering a symbolic description of the dose-response functional relation. Terms in the formula need to be
#' provided in the \code{data} below.
#' @param id an optional vector to specify the id variable for the studies included in the analysis.
#' @param type a vector (or a string) to specify the study-specific design. The values for case-control, incidence-rate, and cumulative incidence data are \code{cc},
#' \code{ir}, and \code{ci}.
#' @param v a vector to specify the variances of the reported log relative risks. Alternatively the user can provide the standard error in the \code{se} argument, or the confidence interval for the reported relative risks 
#' in the \code{lb} and \code{ub} arguments.
#' @param cases a vector to specify the number of cases for each exposure level.
#' @param n a vector to specify the total number of subjects for each exposure level. For incidence-rate data \code{n} indicates the amount of person-time for each exposure level.
#' @param data a data frame (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the previous arguments.
#' @param intercept a logical value to specify if an intercept term needs to be included in the model. See details.
#' @param center a logical value to specify if the exposures level need to be center at the referent ones. See details.
#' @param se an optional vector to specify the standard error of the reported log relative risks; needed if \code{v} is not provided.
#' @param lb an optional vector to specify the lower bound of the confidence interval for the reported relative risks; needed if \code{v} and \code{se} are not provided.
#' @param ub an optional vector to specify the upper bound of the confidence interval for the reported relative risks; needed if \code{v} and \code{se} are not provided.
#' @param covariance method to approximate the coviariance among set of reported log relative risks, "\code{gl}" for the method proposed by Greenland and Longnecker, "\code{h}" for the method proposed by Hamling (default), "\code{fl}" for absolute floated
#'  risks presented by Easton D., "\code{user}" if provided by the user, or "\code{independent}" for assuming independence.
#' @param method  method used to estimate the pooled dose-response relation: "\code{fixed}" for fixed-effects models, "\code{ml}" or "\code{reml}" for random-effects models fitted through (restricted) maximum likelihood, and "\code{mm}" for random-effects models fitted through 
#' method of moments. The default method is "\code{reml}". See \code{\link{mvmeta}}.
#' @param fcov an optional vector to specify the avarage covariances for the set of reported log relative risks. It is required if \code{covariance = "fl"}.
#' @param ucov an optional list of matrices to specify the covariance matrices for the set of reported log relative risks. It is required if \code{covariance = "user"}.
#' @param alpha a scalar to specify the alpha nominal value used in the published data, by defaul equal to .05.
#' @param \dots other useful agurments related to \code{\link{mvmeta}} model. 
#' 
#' @details The function estimates the dose response-relation specified in the \code{formula} for each study included in the analysis. Typically the model does not have an intercept (\code{intercept = FALSE} by default) term since the log relative risk for the exposure 
#' level (usually zero) is zero (RR = 1). For that reason, the values in the desing matrix need to be centered at the referent values, as described by Qin Liu et al, 2009. This is automatically done by the function when \code{center = TRUE} (default value).
#' The study-specific trends are efficienly estimated taking into account the covariance among relative risks. For a theorical description see Orsini et al, 2006.
#' The study specific trends are then combined according to the principles of multivariate random-effects meta-analysis, and relies on \code{\link{mvmeta}} package.
#' 
#' @note The function requires the packages \code{\link{mvmeta}} and \code{\link{aod}} to be installed and loaded.
#' 
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @references
#' Greenland, S.,  Longnecker, M. P. (1992). Methods for trend estimation from summarized dose-response data, with applications to meta-analysis. American journal of epidemiology, 135(11), 1301-1309.
#' 
#' Orsini, N., Bellocco, R.,  Greenland, S. (2006). Generalized least squares for trend estimation of summarized dose-response data. Stata Journal, 6(1), 40.
#'
#' Liu, Q., Cook, N. R., Bergstrom, A., Hsieh, C. C. (2009). A two-stage hierarchical regression model for meta-analysis of epidemiologic nonlinear dose-response data. Computational Statistics & Data Analysis, 53(12), 4157-4167. 
#' 
#' Gasparrini, A., Armstrong, B.,  Kenward, M. G. (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. Statistics in Medicine, 31(29), 3821-3839.
#' 
#' @seealso \code{\link{dosresmeta}}, \code{\link{grl}}, \code{\link{hamling}}
#' 
#' @return The \code{dosresmeta} function typically returns a list of object of class \code{dosresmeta} which resembles a \code{\link{mvmetaObject}}, 
#' with differences in case of trend estimation for a sigle study.
#' 
#' @export dosresmeta
#' 
#' @examples
#'## FIRST EXAMPLE: Single case-control study
#'## Linear trend estimation
#'## Inspect data
#'data("cc_ex")
#'
#'## Fitting the model
#'mod1 <- dosresmeta(formula = logrr ~ dose, type = "cc", cases = case,
#'                   n = n, lb = lb, ub = ub, data= cc_ex)
#'summary(mod1)
#'## Results
#'predict(mod1, delta = 1)
#'
#'
#'## SECOND EXAMPLE: Multiple studies
#'## Linear and quadratic trend using random-effects meta-analysis
#'## Inspect data
#'data("alcohol_cvd")
#'
#'## Linear trend
#'lin <- dosresmeta(formula = logrr ~ dose, type = type, id = id,
#'                  se = se, cases = cases, n = n, data = alcohol_cvd) 
#'## Summarize the results
#'summary(lin)
#'predict(lin, delta = 1)
#'
#'## Non-linear (quadratic) trend
#'quadr <- dosresmeta(formula = logrr ~ dose + I(dose^2), type = type, id = id,
#'                    se = se, cases = cases, n = n, data = alcohol_cvd) 
#'## Summarize the results
#'summary(quadr)
#'
#'## Graphical results
#'with(predict(quadr), {
#'  plot(dose, pred, log = "y", type = "l",
#'       xlim = c(0, 45), ylim = c(.4, 2))
#'  lines(dose,  ci.lb, lty = 2)
#'  lines(dose, ci.ub, lty = 2)
#'  rug(dose, quiet = TRUE)
#'})
#'
dosresmeta <- function(formula, id, type, v, cases, n,
                   data, intercept = F, center = T, se, lb, ub,
                   covariance = "gl", method = "reml", fcov, ucov,
                   alpha = .05, ...)
{
  if (missing(data)) data <- NULL
  if (is.null(data)) {
    data <- sys.frame(sys.parent())
  } else {
    if (!is.data.frame(data)) {
      data <- data.frame(data)
    }
  }
  mf <-  match.call(expand.dots = FALSE)
  mf.id <- mf[[match("id", names(mf))]]
  mf.type <- mf[[match("type", names(mf))]]
  mf.v <- mf[[match("v", names(mf))]]
  mf.cases <- mf[[match("cases", names(mf))]]
  mf.n <- mf[[match("n", names(mf))]]
  mf.se <- mf[[match("se", names(mf))]]
  mf.lb <- mf[[match("lb", names(mf))]]
  mf.ub <- mf[[match("ub", names(mf))]]
  type <- eval(mf.type, data, enclos = sys.frame(sys.parent()))
  if (is.null(type)) type <- as.vector(mf.type)
  id <- as.factor(eval(mf.id, data, enclos = sys.frame(sys.parent())))
  if (is.null(mf.v)){
    se <- eval(mf.se, data, enclos = sys.frame(sys.parent()))
    if (is.null(mf.se)){
      lb <- eval(mf.lb, data, enclos = sys.frame(sys.parent()))
      ub <- eval(mf.ub, data, enclos = sys.frame(sys.parent()))
      v <- ((log(ub) - log(lb)) / (2 * qnorm(1 - alpha/2)))^2
    } else {
      v <- se^2
    }
  } else {
    v <- eval(mf.v, data, enclos = sys.frame(sys.parent()))
  }
  v[is.na(v)] <- 0
  if (is.null(mf.id)) id <- as.factor(rep(1, length(v)))
  cases <- eval(mf.cases, data, enclos = sys.frame(sys.parent()))
  n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
  if (intercept == F) formula <- update(formula, . ~ . + 0)
  mfm <- model.frame(formula, data)
  logrr <- as.matrix(model.response(mfm, "numeric"))
  x <- model.matrix(attr(mfm, "terms"), data = mfm)
   if (center == T){
      for (i in levels(id)){
         x[id == i, ] <- scale(x[id == i, , drop = FALSE], 
                                             t(x[id == i & v== 0, ]), 
                                             scale = FALSE)
      }
   }
  xref <- min(unlist(tapply(x[, 1], id, head, n = 1)))
  fit <- list()
  if (covariance == "gl") {
    cases1 <- grl(logrr, v, cases, n, type, id)$cases1
    n1 <- grl(logrr, v, cases, n, type, id,)$n1
  }
  if (covariance == "h") {
    cases1 <- hamling(logrr, v, cases, n, type, id)$cases1
    n1 <- hamling(logrr, v, cases, n, type, id)$n1
  }
  Ccov <- list()
  vb <- list()
  b <- list()
  for (j in unique(id)) {
    if (covariance %in% c("h", "gl")) {
      if (type[id == j][1] == "cc" | type[id == j][1] == 1) {
        s0 <- 1/cases1[id == j & v == 0] + 1/n1[id == j & v == 0]
        si <- s0 + 1/cases1[id == j & v != 0] + 1/n1[id == j & v != 0]
      }
      if (type[id == j][1] == "ir" | type[id == j][1] == 2) {
        s0 <- 1/cases1[id == j & v == 0]
        si <- s0 + 1/cases1[id == j & v != 0]
      }
      if (type[id == j][1] == "ci" | type[id == j][1] == 3) {
        s0 <- 1/cases1[id == j & v == 0] - 1/n1[id == j & v == 0]
        si <- s0 + 1/cases1[id == j & v != 0] - 1/n1[id == j & v != 0]
      }
      rcorr <- s0/tcrossprod(si)^.5
      diag(rcorr) <- 1
      ccov <- tcrossprod(v[id == j & v != 0])^.5 * rcorr
    }
    if (covariance == "fl") {
      ccov <- matrix(fcov[which(unique(id) == j)], ncol = length(v[id == j & v != 0]),
                     nrow = length(v[id == j & v != 0]))
      diag(ccov) <- v[id == j & v != 0]
    }
    if (covariance == "independent") {
      ccov <- diag(v[id == j & v != 0])
    }
    if (covariance == "user") {
      ccov <- ucov[[which(unique(id) == j)]]
    }
    Ccov <- c(Ccov, list(ccov))
    L <- chol(solve(ccov))
    vbi <- solve(crossprod(L %*% x[id == j & v != 0, ]))
    vb <- c(vb, list(vbi))
    bi <- vbi %*% crossprod(L %*% x[id == j & v != 0, ], L %*% logrr[id == j & v != 0])
    b <- c(b, list(bi))
  }
  tmfm <- solve(t(chol(as.matrix(bdiag(Ccov))))) %*% cbind(logrr, x)[v != 0, ]
  tmod <- lm(tmfm[, 1] ~ 0 + tmfm[, -1])
  colnames(tmfm) <- paste0("t", colnames(as.matrix(mfm)))
  if (length(unique(id)) == 1){
    fit$coefficients <- t(as.matrix(bi))
    fit$vcov <- as.matrix(vbi)
    fit$method <- ""
    fit$dim <- list(k = ncol(x), m = length(unique(id)), p = 1,
                    j = nrow(x))
  }
  if (length(unique(id)) > 1){
    yi <- matrix(unlist(b), ncol = ncol(x), byrow = T)
    colnames(yi) <-  colnames(as.matrix(mfm))[-1]
    fit <- mvmeta(yi, vb, method = method, ...)
  }
  fit$call <- mf
  fit$covariance <- covariance
  fit$design <- model.matrix(attr(mfm, "terms"), data = mfm)
  fit$response <- as.matrix(model.response(mfm, "numeric"))
  fit$ccov <- Ccov
  fit$tdata <- tmfm
  fit$R2 <- summary(tmod)$r.squared
  fit$R2adj <- summary(tmod)$adj.r.squared
  fit$termsi <- attr(mfm, "terms")
  fit$formula <- formula
  fit$coefficients <- matrix(fit$coefficients, ncol = ncol(x), byrow = T)
  colnames(fit$coefficients) <- colnames(as.matrix(mfm))[-1]
  rownames(fit$coefficients) <- ""
  dimnames(fit$vcov) <- rep(list(colnames(fit$coefficients)), 2)
  fit$vb <- vb
  fit$xref <- xref
  class(fit) <- c("dosresmeta", "mvmeta")
  fit
}