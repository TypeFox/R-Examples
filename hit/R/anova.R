#' @title Fast ANOVA
#' 
#' @description A fast sequential analysis of variance (ANOVA). Mainly 
#' developed for internal use.
#' 
#' @param x Design matrix of dimension \code{n * p}.
#' @param y Response vector of observations of length \code{n}.
#' @param assign Integer vector assigning columns to terms can be also given as 
#' \code{x} attribute in which case the argument is ignored. If an intercept 
#' exist it is expected to be the first column in \code{x} and it has to 
#' be specified by a '0' in tis vector. For details about assign see 
#' \code{\link[stats]{model.matrix}}.
#' @param family A description of the error distribution and link function to 
#' be used in the model. For glm this can be a character string naming a family 
#' function or the result of a call to a family function. (See 
#' \code{\link[stats]{family}} for details of family functions.)
#' @param test The name of the test either "LRT" (default) for likelihood ratio 
#' test or "F" for F test.
#' 
#' @seealso \code{\link[stats]{lm}}, \code{\link[stats]{anova}}, and 
#' \code{\link[stats]{aov}}.
#' 
#' @examples
#' y <- rnorm(n = 100)
#' x <- matrix(data = rnorm(1000), nrow = 100)
#' a <- 1:10
#' fast.anova(x = x, y = y, assign = a)
#' 
#' @importFrom stats binomial gaussian Gamma inverse.gaussian poisson quasi 
#' quasibinomial quasipoisson
#' @export 
fast.anova <- function(x, y, assign = NULL, family = gaussian(), 
                       test = c("LRT", "F")) {
  ##### check validity of args
  if (is.null(assign) && !is.null(attr(x, "assign"))) {
    assign <- attr(x, "assign")
  }
  if (is.null(assign)) {
    stop("either a 'x' attribute 'assign' or the 'assign' argument must be specified")
  }
  stopifnot(ncol(x) == length(assign))
  stopifnot(nrow(x) == length(y))
  if (is.character(family)) {
    family <- eval(call(family))
  } else if (is.function(family)) {
    family <- family()
  }
  test <- match.arg(test, c("LRT", "F"))
  # fit anova model
  if (family$family == "gaussian" && test == "F") 
    p <- fast.lmanova(x, y, assign)
  else
    p <- fast.glmanova(x, y, assign, family, test)
  p
}


#' @title Fast LM F test ANOVA
#' 
#' @description A fast sequential analysis of variance (ANOVA). Mainly 
#' developed for internal use.
#' 
#' @param x Design matrix of dimension \code{n * p}.
#' @param y Response vector of observations of length \code{n}.
#' @param assign Integer vector assigning columns to terms.
#' 
#' @importFrom stats lm.fit pf 
#' @keywords internal
fast.lmanova <- function(x, y, assign) {
  # LM fit by pivoted QR decomposition
  fit <- lm.fit(x, y)
  if (assign[1L] == 0L) { 
    ##### with intercept
    full.rank <- 1L:(fit$rank - 1L)
    assign.pivot <- assign[fit$qr$pivot[full.rank + 1L]]
    var <- fit$effects[-1L]^2
  } else { 
    ##### without intercept
    full.rank <- 1L:fit$rank
    assign.pivot <- assign[fit$qr$pivot[full.rank]]
    var <- fit$effects^2
  }
  # Treatment: Sum Sq | DF | Mean Sq
  ss.treat <- tapply(var[full.rank], assign.pivot, "sum")
  df.treat <- table(assign.pivot)
  ms.treat <- ss.treat / df.treat
  # Residuals: Sum Sq | DF | Mean Sq
  ss.res <- sum(var[-full.rank])
  df.res <- fit$df.residual
  ms.res <- ss.res / df.res
  # F value
  f <- ms.treat / ms.res 
  # p value
  p <- rep(1, max(assign))
  p[unique(assign.pivot)] <- pf(f, df.treat, df.res, lower.tail = FALSE)
  p
}


#' @title Fast GLM F test of LR test ANOVA
#' 
#' @description A fast sequential analysis of variance (ANOVA). Mainly 
#' developed for internal use.
#' 
#' @param x Design matrix of dimension \code{n * p}.
#' @param y Response vector of observations of length \code{n}.
#' @param assign Integer vector assigning columns to terms.
#' @param family A description of the error distribution and link function to 
#' be used in the model. For glm this must be the result of a call to a family 
#' function.
#' @param test The name of the test either "LRT" (default) for likelihood ratio 
#' test or "F" for F test.
#' 
#' @importFrom speedglm speedglm.wfit
#' @importFrom stats pf pchisq
#' @keywords internal
fast.glmanova <- function(x, y, assign, family, test) {
  inter <- ifelse(assign[1] == 0L, TRUE, FALSE)
  # Full model fit by pivoted Colesky decomposition
  suppressWarnings(
    full.fit <- speedglm.wfit(X = x, y = y, 
                              intercept = inter, family = family, 
                              method = "Cholesky", tol.solve = 1e-30)
  )
  # Reduced model fits by pivoted Colesky decomposition
  red.dev <- red.df <- NULL 
  if (max(assign) > 1L) {
    for (i in seq_len( max(assign) - 1L)) {
      suppressWarnings(
        red.fit <- speedglm.wfit(X = x[, assign <= i, drop = FALSE], y = y, 
                                 intercept = inter, family = family, 
                                 method = "Cholesky", tol.solve = 1e-30)
      )
      red.dev <- c(red.dev, red.fit$deviance)
      red.df <- c(red.df, red.fit$df)
    }
  }
  # Treatment: Variance | DF
  df.treat <- -diff(c(full.fit$nulldf, red.df, full.fit$df))
  var.treat <- pmax(-diff(c(full.fit$nulldev, red.dev, full.fit$deviance)))
  ##### check for singularity
  nonsing.covar <- which(df.treat != 0)
  df.treat <- df.treat[nonsing.covar]
  var.treat <- var.treat[nonsing.covar]
  ##### Dispersion factor
  disp <- full.fit$dispersion
  # estimated  test statistic
  p <- rep(1, max(assign))
  if (test == "F") {
    #  residuals: DF
    df.res <- ifelse(disp == 1, Inf, full.fit$df)
    # F value
    f <- var.treat / (disp * df.treat)
    # p value
    p[nonsing.covar] <- pf(f, df.treat, df.res, lower.tail = FALSE)
  } else {
    # Chi value
    chi <- var.treat / disp
    # p value
    p[nonsing.covar] <- pchisq(chi, df.treat, lower.tail = FALSE)
  }
  p
}
