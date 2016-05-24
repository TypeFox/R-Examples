#' @title flexmix model driver for adaptive lasso (elastic-net) with GLMs
#' @author F. Mortier (fmortier@cirad.fr) and N. Picard (nicolas.picard@cirad.fr)
#' @param formula    A symbolic description of the model to be fit. 
#' The general form is y~x|g where y is the response, x the set of predictors and g an 
#' optional grouping factor for repeated measurements.
#' @param family a description of the error distribution and link function to be used in the model.
#' "gausian", "poisson" and "binomial" are allowed. 
#' @param adaptive boolean indicating if algorithm should perform adaptive lasso or not
#' @param select boolean vector indicating which covariates will be included in the selection process.
#'  Others will be included in the model.
#' @details Some care is needed to ensure convergence of the
#' algorithm, which is computationally more challenging than a standard EM. 
#' In the proposed method, not only are cluster allocations identified
#' and component parameters estimated as commonly done in mixture models,
#' but there is also variable selection via penalized regression using
#' $k$-fold cross-validation to choose the penalty parameter. 
#' For the algorithm to converge, it is necessary that the same cross-validation
#' partitioning be used across the EM iterations, i.e.,
#' the subsamples for cross-validation must be defined at the beginning
#' This is accomplished using the {\tt foldid} option
#'  as an additional parameter to be passed to  \code{\link{cv.glmnet}} (see \link{glmnet} package documentation).
 
FLXMRglmnet <-
    function(formula = .~., family = c("gaussian", "binomial", "poisson"), adaptive = TRUE, select = TRUE, offset = NULL, ...) {
        family <- match.arg(family)
        z <- FLXMRglm(formula = formula, family = family)
        z@preproc.x <- function(x) {
            if (!isTRUE(all.equal(x[, 1], rep(1, nrow(x)), check.attributes = FALSE)))
                stop("The model needs to include an intercept in the first column.")
            x
        }
        z@fit <- function(x, y, w) {
            if (all(!select)) {
                coef <- if (family == "gaussian")
                            lm.wfit(x, y, w = w)$coef
                        else if (family == "binomial")
                            glm.fit(x, y, family = binomial(), weights = w)$coef
                        else if (family == "poisson")
                            glm.fit(x, y, family=poisson(), weights = w)$coef
            } else {
                if (adaptive) {
                    coef <- if (family == "gaussian")
                                lm.wfit(x, y, w = w)$coef[-1]
                            else if(family == "binomial")
                                glm.fit(x, y, family = binomial(), weights = w)$coef[-1]
                            else if (family == "poisson")
                                glm.fit(x, y, family = poisson(), weights = w)$coef[-1]
                    penalty <- mean(w) / abs(coef)
                } else
                    penalty <- rep(1, ncol(x) - 1)
                if (any(!select)){
                    select <- which(!select)
                    penalty[select] <- 0
                }      
                m <-  glmnet::cv.glmnet(x[, -1, drop = FALSE], y, family = family, weights = w,
                                        penalty.factor = penalty, ...)
                coef <- as.vector(coef(m, s = "lambda.min"))
            }
            df <- sum(coef != 0)
            sigma <- if (family == "gaussian") sqrt(sum(w * (y - x %*% coef)^2/mean(w))/(nrow(x) - df)) else NULL
            with(list(coef = coef, sigma = sigma, df = df + ifelse(family == "gaussian", 1, 0)),
                 eval(z@defineComponent))
        }
        z
    }
    


