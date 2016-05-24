#' Generate and plot decision curves.
#'
#'Decision curves are a useful tool to evaluate the population impact of adopting a risk prediction instrument into clinical practice. Given one or more instruments (risk models) that predict the probability of a binary outcome, this package calculates and plots decision curves, which display estimates of the standardized net benefit by the probability threshold used to categorize observations as 'high risk.' Curves can be estimated using data from an observational cohort, or from case-control studies when an estimate of the population outcome prevalence is available. Confidence intervals calculated using the bootstrap can be displayed and a wrapper function to calculate cross-validated curves using k-fold cross-validation is also provided.
#'
#' Functions in this package include:
#'  \code{\link{decision_curve}},
#'  \code{\link{summary.decision_curve}},
#'  \code{\link{plot_decision_curve}},
#'  \code{\link{plot_clinical_impact}},
#'  \code{\link{plot_roc_components}}, and
#'  \code{\link{Add_CostBenefit_Axis}}.
#'
#'
"_PACKAGE"
#> [1] "_PACKAGE"
