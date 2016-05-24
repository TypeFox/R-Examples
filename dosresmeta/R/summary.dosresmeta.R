#' summarizing dosresmeta Models
#' 
#' @description Print and summary method functions for dose-response models represented in objects of class "\code{dosresmeta}".
#' 
#' @param object an object of class \code{dosresmeta} produced by \code{\link{dosresmeta}}.
#' @param x an object of class \code{dosresmeta} or \code{summary.dosresmeta} produced by \code{\link{dosresmeta}} or \code{summary.dosresmeta}, respectively.
#' @param ci.level the confidence level used for defining the confidence intervals for the estimates of the (fixed-effects) coefficients.
#' @param digits an integer specifying the number of digits to which printed results must be rounded.
#' @param \dots further arguments passed to or from other methods.
#'
#' @details the \code{print} method for class \code{dosresmeta} only returns basic information of the fitted model, namely the call, 
#' estimated (fixed-effects) coefficients, and dimensions).
#' 
#' The \code{summary} method function computes additional statistics and tests, and produces a list object of class \code{summary.dosresmeta}. 
#' The \code{print} method function for this class, depending on the number of studies included in the analysis, shows additional information, 
#' such as tables reporting the estimates for the fixed and random-effects parts of the model, Chi-square test for model significance, 
#' Cochran Q test for heterogeneity and I-square. 
#' 
#' @return The \code{summary} method function for \code{dosresmeta} objects produces a list of class "\code{summary.dosresmeta}" which resembles
#' a list object of class \code{\link{summary.mvmeta}}.
#' 
#'
#' As usual, the \code{print} method functions for classes "\code{dosresmeta}" and "\code{summary.dosresmeta}" do not return any value.
#'
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' @seealso \code{\link{dosresmeta}}, \code{\link{summary}}
#' 
#' @examples
#' ## Load data and run the model
#'data("alcohol_cvd")
#'model <- dosresmeta(formula = logrr ~ dose + I(dose^2), type = type, id = id,
#'                    se = se, cases = cases, n = n, data = alcohol_cvd) 
#'## Defult print
#'model
#'## Specify digits
#'print(model, digit = 2)
#'## summary with 90th confidence intervals
#'summary(model, ci.level = .8)
#'
#' @rdname summary.dosresmeta
#' @method print dosresmeta
#' @export print.dosresmeta
#' @S3method print dosresmeta
#' 

print.dosresmeta <- function (x, digits = 4, ...) 
{
  cat("Call:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Fixed-effects coefficients:", "\n", sep = "")
  table <- formatC(x$coefficients, digits = digits, format = "f")
  print(table, quote = FALSE, right = TRUE, print.gap = 2)
  cat("\n")
  if (x$dim$m > 1){
    cat(x$dim$m, " studies, ", x$df$nall, " values, ", x$df$fixed, " fixed and ", 
        x$df$random, " random-effects parameters", sep = "")
  }
  cat("\n")
}

#' @rdname summary.dosresmeta
#' @method summary dosresmeta
#' @export summary.dosresmeta
#' @S3method summary dosresmeta
#' 
summary.dosresmeta <- function (object, ci.level = 0.95, ...) 
{
  if (ci.level <= 0 || ci.level >= 1)  stop("'ci.level' must be within 0 and 1")
  coef <- object$coefficients
  vcov <- object$vcov
  dim <- object$dim
  Psi <- if (dim$m > 1){ 
    object$Psi
  } else NULL 
  lab <- object$lab
  coef <- as.numeric(coef)
  coef.se <- sqrt(diag(vcov))
  zval <- coef/coef.se
  zvalci <- qnorm((1 - ci.level)/2, lower.tail = FALSE)
  pvalue <- 2 * (1 - pnorm(abs(zval)))
  ci.lb <- coef - zvalci * coef.se
  ci.ub <- coef + zvalci * coef.se
  cilab <- paste(signif(ci.level, 2) * 100, "%ci.", c("lb", "ub"), sep = "")
  test <- wald.test(b = coef(object), Sigma = vcov(object), 
                    Terms = 1:ncol(object$coefficients))[["result"]][["chi2"]]
  tabfixed <- cbind(coef, coef.se, zval, pvalue, ci.lb, ci.ub)
  dimnames(tabfixed) <- list(colnames(object$coefficients) , 
                             c("Estimate", "Std. Error", "z", "Pr(>|z|)", cilab))
  corFixed <- vcov/outer(coef.se, coef.se)
  if (dim$m > 1){
    corRandom <- if (object$method != "fixed") {
      ran.sd <- sqrt(diag(Psi))
      Psi/outer(ran.sd, ran.sd)
    } else NULL
    qstat <- unclass(qtest(object))
    keep <- match(c("vcov", "Psi", "df.res", "rank", "logLik", 
                    "converged", "niter", "negeigen", "dim", "df", "lab", 
                    "na.action", "call", "terms", "method", "covariance"), names(object), 0L)
    out <- c(list(coefficients = tabfixed), object[keep], list(AIC = AIC(object), 
                                                               BIC = BIC(object), corFixed = corFixed, corRandom = corRandom, 
                                                               qstat = qstat, ci.level = ci.level, test = test))
  }
  if (dim$m == 1){
    Q <- tcrossprod(crossprod((object$response - tcrossprod(object$design, 
                                                            object$coefficients))[-1], t(chol(solve(object$ccov[[1]])))))
    keep <- match(c("vcov", "Psi", "df.res", "dim", "df", "lab", 
                    "call", "method", "covariance"), names(object), 0L)
    out <- c(list(coefficients = tabfixed), object[keep], Q = Q, test = list(test))
  }
  class(out) <- "summary.dosresmeta"
  return(out)
}

#' @rdname summary.dosresmeta
#' @method print summary.dosresmeta
#' @export print.summary.dosresmeta
#' @S3method print summary.dosresmeta
#' 
print.summary.dosresmeta <- function (x, digits = 4, ...) 
{
  methodname <- c("reml", "ml", "fixed", "mm", "vc")
  methodlabel <- c("REML", "ML", "Fixed", "Method of moments", "Variance components")
  covariancename <- c("h", "gl", "fl", "user")
  covariancelabel <- c("Hamling", "Greenland & Longnecker", "Floated Absolute Risks", 
                       "User defined")
  int <- 1
  cat("Call:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (x$dim$m > 1){ 
    cat(if (x$dim$k == 1L) "Uni" else "Multi", "variate ", ifelse(x$method == "fixed", 
                                                                  "fixed", "random"), "-effects meta-", ifelse(x$dim$p - int > 0, 
                                                                                                               "regression", "analysis"), "\n", sep = "")
    cat("Dimension: ", x$dim$k, "\n", sep = "")
    if (x$method != "fixed") {
      cat("Estimation method: ", methodlabel[which(x$method == 
                                                     methodname)], "\n", sep = "")
      cat("Variance-covariance matrix Psi: ", "unstructured", "\n", sep = "")
    }
  }
  if (x$dim$m == 1){ 
    cat("Trend estimation", "\n", sep = "")
    cat("Dimension: ", x$dim$k, " study", "\n", sep = "")
  }
  cat("Approximate covariance method: ", covariancelabel[which(x$covariance == 
                                                                 covariancename )], "\n", sep = "")
  cat("\n")
  cat("Fixed-effects coefficients", "\n", sep = "")
  signif <- symnum(x$coefficients[, "Pr(>|z|)"], corr = FALSE, 
                   na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                   symbols = c("***", "**", "*", ".", " "))
  tabletot <- formatC(x$coefficients, digits = digits, format = "f")
  tabletot <- cbind(tabletot, signif)
  colnames(tabletot)[7] <- ""
  ## rownames(tabletot) <- colnames(x$coefficients)
  print(tabletot, quote = FALSE, right = TRUE, print.gap = 2)
  cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")
  chi2 <- formatC(x$test[1], digits = digits, format = "f")
  pchi2 <- formatC(x$test[3], digits = digits, format = "f")
  cat("Chi2 model: X2 = ", chi2, " (df = ", x$test[2], "), p-value = ", 
      pchi2, "\n", sep = "")
  if (x$dim$m == 1){
    Q <- formatC(x$Q, digits = digits, format = "f")
    pQ <- 1 - pchisq(x$Q, x$dim$j - 2 )
    pQ <- formatC(pQ, digits = digits, format = "f")
    cat("Goodness-of-fit (chi2): X2 = ", Q, " (df = ", x$dim$j - 2, "), p-value = ", 
        pQ, "\n", sep = "")
  }
  if(x$dim$m > 1){
    Q <- formatC(x$qstat$Q, digits = digits, format = "f")
    pvalue <- formatC(x$qstat$pvalue, digits = digits, format = "f")
    i2 <- formatC(pmax((x$qstat$Q - x$qstat$df)/x$qstat$Q * 100, 0), digits = 1, format = "f")
    cat(if (x$qstat$k == 1) "Uni" else "Multi", "variate ", "Cochran Q-test for ", 
        if (x$qstat$residual) "residual ", "heterogeneity:", "\n", sep = "")
    cat("Q = ", Q[1], " (df = ", x$qstat$df[1], "), p-value = ", pvalue[1], "\n", sep = "")
    cat("I-square statistic = ", i2[1], "%", "\n\n", sep = "")
    cat(x$dim$m, " studies, ", x$df$nall, " observations, ", 
        x$df$fixed, " fixed and ", x$df$random, " random-effects parameters", 
        sep = "")
  }
  cat("\n")
}