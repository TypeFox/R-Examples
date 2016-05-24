#' Diagnostic tools for hierarchical (multilevel) linear models
#' 
#' HLMdiag provides a suite of diagnostic tools for hierarchical 
#' (multilevel) linear models fit using \code{\link{lmer}}. 
#' These tools are grouped below by purpose. 
#' See the help documentation for additional information
#' about each function.
#' 
#' \bold{Residual analysis}
#' 
#' HLMdiag's \code{\link{HLMresid}} function provides a convenient
#' wrapper to obtain residuals at each level of a hierarchical 
#' linear model. In addition to being a wrapper function for functions 
#' implemented in the \code{lme4} package, HLMresid provides access
#' to the marginal and least squares residuals (through \code{\link{LSresids}}) 
#' that were not previously implemented.
#' 
#' \bold{Influence analysis}
#' 
#' Influence on fitted values
#' 
#' HLMdiag provides \code{\link{leverage}} that calculates the influence
#' that observations/groups have on the fitted values (leverage). 
#' For mixed/hierarchical models leverage can be decomposed into two parts: the 
#' fixed part and the random part. We refer the user to the references
#' cited in the help documentation for additional explanation.
#' 
#' Influence on fixed effects estimates
#' 
#' HLMdiag provides \code{\link{cooks.distance}} and \code{\link{mdffits}}
#' to assess the influence of subsets of observations on the fixed effects.
#' 
#' Influence on precision of fixed effects
#' 
#' HLMdiag provides \code{\link{covratio}} and \code{\link{covtrace}}
#' to assess the influence of subsets of observations on the precision of
#' the fixed effects.
#' 
#' Influence on variance components
#' 
#' HLMdiag's \code{\link{rvc}} calculates the relative variance change to
#' assess the influence of subsets of observations on the variance
#' components.
#' 
#' \bold{Graphics}
#' 
#' HLMdiag also strives to make graphical assessment easier in the 
#' \code{ggplot2} framework by providing dotplots for influence diagnostics
#' (\code{\link{dotplot_diag}}), grouped Q-Q plots (\code{\link{group_qqnorm}}),
#' and Q-Q plots that combine the functionality of \code{\link{qqnorm}} and
#' \code{\link{qqline}} (\code{\link{ggplot_qqnorm}}).
#' 
#' @useDynLib HLMdiag
#' @import lme4
#' @import Matrix
#' # @import methods
#' @import ply
#' @import reshape2
#' @importFrom stats4 coef, confint, plot
#' @importFrom stats cooks.distance, covratio
#' @docType package
#' @name HLMdiag
#' @aliases HLMdiag package-HLMdiag
#' @keywords package
NULL