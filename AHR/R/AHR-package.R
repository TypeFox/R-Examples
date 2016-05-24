#' AHR
#'
#' This package provides methods for estimation of multivariate average hazard ratios as defined by Kalbfleisch and Prentice.
#' The underlying survival functions of the event of interest in each group can be estimated using either the (weighted) Kaplan-Meier estimator or
#' the Aalen-Johansen estimator for the transition probabilities in Markov multi-state models. Right-censored and left-truncated data is supported.
#' Moreover, the difference in restricted mean survival can be estimated. Currently variance estimation for the average hazard ratio based on the
#' Aalen-Johansen estimator is only supported for competing risks models, i.e. for estimation of the average sub-distribution hazard ratio
#' (Average cause-specific hazard ratios can be estimated by using the Kaplan-Meier estimator with competing risks data).
#'
#' Furthermore estimation of quantiles, ratios and differences of quantiles and corresponding p-values and confidence intervals of
#' survival times based on the (weighted) Kaplan-Meier estimator and the Aalen-Johansen estimator is also supported.
#' 
#' @author Matthias Brueckner \email{mwbr@@math.uni-bremen.de}
#' @references
#' J.~D. Kalbfleisch and R.~L. Prentice. Estimation of the average hazard ratio. \emph{Biometrika}, 68(1):105--112, Apr. 1981.
#' 
#' S.~Murray and A.~A. Tsiatis. Nonparametric survival estimation using prognostic longitudinal covariates. \emph{Biometrics}, 52(1):137--151, Mar. 1996.
#' 
#' C.~A. Struthers and J.~D. Kalbfleisch. Misspecified proportional hazard models. \emph{Biometrika}, 73(2):363--369, Aug. 1986.
#' 
#' @examples
#' T <- c(rexp(100, 1), rexp(100, 2))
#' C <- c(rexp(100, 1), rexp(100, 2))
#' Y <- pmin(T, C)
#' D <- T <= C
#' Z <- rep(c(0,1), c(100, 100))
#'
#' ## uses Kaplan-Meier estimator by default
#' fit <- avgHR(2, data.frame(Y=Y, D=D, Z=Z), formula=Surv(Y, D) ~ Z)
#' fit
#' 
#' ## same as
#' \dontrun{fit <- avgWKM(2, data.frame(Y=Y, D=D, Z=Z), formula=Surv(Y, D) ~ Z)}
#'
#' ## use bootstrap to estimate covariance matrix
#' \dontrun{fit <- avgWKM(2, data.frame(Y=Y, D=D, Z=Z), formula=Surv(Y, D) ~ Z, cov=FALSE,
#'                        bootstrap=10000)}
#'
#' ## calculate restricted mean difference
#' rdm <- rmeanDiff.ahr(fit)
#' rdm
#'
#' ## ventilation status in intensive care unit patients dataset from etm package
#' library(etm)
#' data(sir.cont)
#' df <- sir.cont
#' df$Trt <- factor(rep(0, nrow(df)), levels=c(0, 1))
#' ids <- unique(df$id)
#' df$Trt[df$id %in% sample(ids, floor(length(ids)/2), FALSE)] <- 1
#' 
#' # transition matrix
#' tra <- matrix(FALSE, nrow=3, ncol=3)
#' tra[1, 2:3] <- TRUE
#' tra[2, c(1, 3)] <- TRUE
#'
#' # NOTE: variance estimation not yet supported for Aalen-Johansen based avg. HR
#' sc.fit <- avgHR(2, method="aj", data=df, target="0 2", states=c("0", "1", "2"), transitions=tra,
#'                 censoring="cens", cov=FALSE)
#' sc.fit
#' 
#' @name AHR
#' @docType package
#' @import Rcpp survival stats
#' @importFrom MASS ginv
#' @importFrom etm etm
#' @useDynLib AHR
NULL
