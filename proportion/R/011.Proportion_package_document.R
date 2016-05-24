#' Proportion: Let x denote the number of successes in n independent Bernoulli trials with
#' X ~ Binomial (n, p) then phat  = x/n denotes the sample proportion.
#' @section Introduction:
#' Objective of this package is to present interval estimation procedures for 'p'
#' outlined above in a more comprehensive way.Quality assessment procedures such as statistic
#' based on coverage probability, Expected length, Error, p-confidence and p-bias are also included.
#' Also, an array of Bayesian computations (Bayes factor, Empirical Bayesian,
#' Posterior predictive computation, and posterior probability) with conjugate prior are made available.
#' The proportion package provides three categories of important functions:
#' \strong{Confidence Intervals}, \strong{metrics on confidence intrvals} (coverage probability, length, p-confidence and
#' p-bias, error and long term power) and \strong{other methods} (hypothesis
#' testing and general/simulation methods).
#' @section Proportion methods grouping:
#' For finding confidence interval for p we have included
#' \itemize{
#' \item  Methods based on the asymptotic normality of the sample proportion and
#' estimating standard error
#' \item  Exact methods based on inverting equal-tailed binomial tests of H0 : p = p0,
#' \item  Methods based on likelihood ratios
#' \item  Bayesian approaches with beta priors or other suitable priors.}
#' @section Proportion function naming convention:
#' The general guideline for finding functions are given below:
#' \itemize{
#' \item Short names for concepts: ci - Confidence Interval, covp - Coverage Probability,
#' expl - Expected length (simulation), length - Sum of length, pCOpBI - p-Confidence and p-Bias,
#' err - Error and long term power
#' \item Short names for methods: AS - ArcSine, LR - Likelihood Ratio, LT - Logit Wald,
#' SC - Score (also know as Wilson), TW - Wald-T, WD - Wald, BA - Bayesian and
#' EX - Exact in general form that includes Mid-P and Clopper-Pearson.
#' \item For adjusted methods "A" is added to the function name while "C" will be added if it is
#' continuity corrected.
#' \item For generic functions BAF - Bayesian Factor, SIM - Simulation, GEN - Generic, PRE - Predicted,
#' POS - Posterior
#' \item Combining the above you should be able to identify the function. For example,
#' function for coverage probability (covp) using ArcSine (AS) method will be covpAS(). If we need
#' the adjusted coverage probability (covp) using ArcSine (AS) method, then it will be covpAAS().
#' \item	Wherever possible, results are consolidated for all
#' \code{x (0,1...n)} and specific \code{x}
#' (function name succeeds with \code{x}). For example, if we run \code{ciAS(n=5, alp=0.05)}
#' the output of \code{x=5} will be the same as
#' \code{ciASx(x=5, n=5,alp=0.05)}.
#' In the first case the output is printed for all the values of \code{x} till \code{x=n}.
#' \item All refers to six approximate methods (Wald, Score, Likelihood Ratio, ArcSine,
#' Logit Wald and Wald-T) - AAll (Adjusted All) refers to six  methods adjusted with
#' adding factor \code{h} (Wald, Score, Likelihood Ratio, ArcSine, Logit Wald and Wald-T)
#' \item CAll (Continuity corrected All) refers to five methods (Wald, Score, ArcSine,
#' Logit Wald and Wald-T) with continuity correction \code{c}
#' \item Grouping functions for plots end with "g" (PlotciAllxg is the same as PlotciAllx,
#' except the results are grouped by x)
#' \item For almost all the functions, corrosponding plot function is implemented,
#' which plots the output in an apporiate graph. For example, the function
#' \code{covpAll()} will give the numeric output for the
#' coverage probability of the six approximate methods (see explanation of All above).
#' Prefixing this with Plot makes it
#' \code{PlotcovpAll()}
#' and will display the plot for the same six approximate methods.}
#' @section Reproducibility of reference papers:
#' To help the researcher reporduce results in existing papers we have taken six key papers
#' (see references below) [3], [8], [9], [10], [11], [12] and reproduced the results and
#' suggested further
#' items to try. Details are in the vignette.
#' @docType package
#' @references
#' [1] 1993 Vollset SE.
#' Confidence intervals for a binomial proportion.
#' Statistics in Medicine: 12; 809 - 824.
#'
#' [2] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [3] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [4] 2001 Brown LD, Cai TT and DasGupta A.
#' Interval estimation for a binomial proportion.
#' Statistical Science: 16; 101 - 133.
#'
#' [5] 2002 Pan W.
#' Approximate confidence intervals for one proportion and difference of two proportions
#' Computational Statistics and Data Analysis 40, 128, 143-157.
#'
#' [6] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#'
#' [7] 2014 Martin Andres, A. and Alvarez Hernandez. M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#'
#' [8] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#'
#' [9] 2005 Joseph L and Reinhold C.
#' Statistical Inference for Proportions
#' American Journal of Radiologists 184; 1057 - 1064
#'
#' [10] 2008 Zhou, X. H., Li, C.M. and Yang, Z.
#' Improving interval estimation of binomial proportions.
#' Phil. Trans. R. Soc. A, 366, 2405-2418
#'
#' [11] 2012 Wei Yu, Xu Guo and Wangli Xua.
#' An improved score interval with a modified midpoint for a binomial proportion,
#' Journal of Statistical Computation and Simulation, 84, 5, 1-17
#'
#' [12] 2008 Tuyl F, Gerlach R and Mengersen K .
#' A comparison of Bayes-Laplace, Jeffreys, and Other Priors: The case of zero events.
#' The American Statistician: 62; 40 - 44.
#' @import stats
#' @name A package for binomial proportion
NULL
#> NULL
