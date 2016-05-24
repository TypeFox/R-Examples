#' Longitudinal Targeted Maximum Likelihood Estimation
#' 
#' \code{ltmle} is Targeted Maximum Likelihood Estimation (TMLE) of
#' treatment/censoring specific mean outcome for point-treatment and
#' longitudinal data. \code{ltmleMSM} adds Marginal Structural Models. Both
#' always provide Inverse Probability of Treatment/Censoring Weighted estimate
#' (IPTW) as well. Maximum likelihood based G-computation estimate (G-comp) can
#' be obtained instead of TMLE. \code{ltmle} can be used to calculate additive
#' treatment effect, risk ratio, and odds ratio.
#' 
#' The estimates returned by \code{ltmle} are of a treatment specific mean,
#' \eqn{E[Y_{\bar{a}}]}, the mean of the final treatment node, where all
#' treatment nodes, \eqn{A}, are set to \eqn{\bar{a}} (\code{abar}) and all
#' censoring nodes \eqn{C} are set to 1 (uncensored). The estimates returned by
#' \code{ltmleMSM} are similar but are the parameters in a working marginal
#' structural model.
#' 
#' \code{data} should be a data frame where the order of the columns
#' corresponds to the time-ordering of the model.  \itemize{ \item in censoring
#' columns (Cnodes): factor with two levels: "censored" and "uncensored". The
#' helper function \code{CensoringToBinary} can be used to create these
#' factors.  \item in treatment columns (Anodes): 1 = treated, 0 = untreated
#' (must be binary) \item in event columns (Ynodes): If \code{survivalOutcome}
#' is \code{TRUE}, then Y nodes are treated as indicators of a one-time event.
#' See details for \code{survivalOutocme}. If \code{survivalOutcome} is
#' \code{FALSE}, Y nodes are treated as binary if all values are 0 or 1, and
#' are treated as continuous otherwise. If Y nodes are continuous, they may be
#' automatically scaled. See details for \code{Yrange}.  \item time-dependent
#' covariate columns (Lnodes): can be any numeric data \item Data in
#' \code{Cnodes}, \code{Anodes}, \code{Lnodes} and \code{Ynodes} are not used
#' after (to the right of) censoring (or an event when
#' \code{survivalOutcome==TRUE}) and may be coded as \code{NA} or any other
#' value.  \item Columns in \code{data} that are before (to the left of) the
#' first of \code{Cnodes} or \code{Anodes} are treated as baseline variables,
#' even if they are specified as \code{Lnodes}.  \item After the first of
#' \code{Cnodes}, \code{Anodes}, \code{Ynodes}, or \code{Lnodes}, every column
#' must be in one of \code{Cnodes}, \code{Anodes}, \code{Ynodes}, or
#' \code{Lnodes}.  }
#' 
#' If \code{survivalOutcome} is \code{TRUE}, all Y values are indicators of an
#' event (e.g. death) at or before the current time, where 1 = event and 0 = no
#' event. The events in Ynodes must be of the form where once Y jumps to 1, Y
#' remains 1 at subsequent nodes.
#' 
#' For continuous outcomes, (\code{survivalOutcome==FALSE} and some Y nodes are
#' not 0 or 1,) Y values are truncated at the minimum and maximum of
#' \code{Yrange} if specified, and then transformed and scaled to be in [0,1].
#' That is, transformed to \code{(Y-min(Yrange))/(max(Yrange)-min(Yrange))}. If
#' \code{Yrange} is \code{NULL}, it is set to the range of all Y nodes. In that
#' case, Y nodes are only scaled if any values fall outside of [0,1]. For
#' intervention specific means (\code{ltmle}), parameter estimates are
#' transformed back based \code{Yrange}.
#' 
#' \code{Qform} should be \code{NULL}, in which case all parent nodes of each L
#' and Y node will be used as regressors, or a named character vector that can
#' be coerced to class "\code{\link{formula}}". The length of \code{Qform} must
#' be equal to \code{length(Lnodes) + length(Ynodes)}** and the names and order
#' of the formulas must be the same as the names and order of the L and Y nodes
#' in \code{data}. The left hand side of each formula should be
#' "\code{Q.kplus1}". If \code{SL.library} is \code{NULL}, \code{glm} will be
#' called using the elements of \code{Qform}. If \code{SL.library} is
#' specified, \code{\link[SuperLearner:SuperLearner]{SuperLearner}} will be
#' called and all variables appearing on the right hand side of a formula will
#' be passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. The
#' actual functional form of the formula is unimportant if
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}} is called.
#' 
#' ** If there is a "block" of L and Y nodes not separated by A or C nodes,
#' only one regression is required at the first L/Y node in a block. You can
#' pass regression formulas for the other L/Y nodes, but they will be ignored
#' (with a message). See example 5.
#' 
#' \code{gform} should be \code{NULL}, in which case all parent nodes of each L
#' and Y node will be used as regressors, or a character vector that can be
#' coerced to class "\code{\link{formula}}", or a matrix/array of Prob(A=1). If
#' \code{gform} is a character vector, the length of \code{gform} must be equal
#' to \code{length(Anodes) + length(Cnodes)} and the order of the formulas must
#' be the same as the order the A and C nodes appear in \code{data}. The left
#' hand side of each formula should be the name of the Anode or Cnode. If
#' \code{SL.library} is \code{NULL}, \code{glm} will be called using the
#' elements of \code{gform}. If \code{SL.library} is specified,
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}} will be called and all
#' variables appearing on the right hand side of a formula will be passed to
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. The actual functional
#' form of the formula is unimportant if
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}} is called.
#' 
#' In \code{ltmle}, \code{gform} can also be a n x numACnodes matrix where
#' entry (i, j) is the probability that the ith observation of the jth A/C node
#' is 1 (if an Anode) or uncensored (if a Cnode), conditional on following abar
#' up to that node. In \code{ltmleMSM}, \code{gform} can similarly be a n x
#' numACnodes x numRegimes array, where entry (i, j, k) is the probability that
#' the ith observation of the jth A/C node is 1 (if an Anode) or uncensored (if
#' a Cnode), conditional on following regime k up to that node. If \code{gform}
#' is a matrix/array, \code{deterministic.g.function} will not be used and
#' should be \code{NULL}.
#' 
#' \code{abar} specifies the counterfactual values of the Anodes, using the
#' order they appear in \code{data} and should have the same length (if abar is
#' a vector) or number of columns (if abar is a matrix) as \code{Anodes}.
#' 
#' \code{rule} can be used to specify a dynamic treatment rule. \code{rule} is
#' a function applied to each row of \code{data} which returns the a numeric
#' vector of the same length as \code{Anodes}.
#' 
#' \code{abar} and \code{rule} cannot both be specified. If one of them if a
#' list of length 2, additive treatment effect, risk ratio, and odds ratio can
#' be computed using \code{\link{summary.ltmleEffectMeasures}}.
#' 
#' \code{regimes} can be a binary array: n x numAnodes x numRegimes of
#' counterfactual treatment or a list of 'rule' functions as described above
#' for the \code{rule} parameter for the \code{ltmle} function
#' 
#' \code{deterministic.g.function} can be a function used to specify model
#' knowledge about value of Anodes and/or Cnodes that are set
#' deterministically. For example, it may be the case that once a patient
#' starts treatment, they always stay on treatment. For details on the form of
#' the function and examples, see
#' \code{\link{deterministic.g.function_template}}
#' 
#' \code{deterministic.Q.function} can be a function used to specify model
#' knowledge about the final event state. For example, it may be the case that
#' a patient can complete the study at some intermediate time point, in which
#' case the probability of death is 0 (assuming they have not died already).
#' For details on the form of the function and examples, see
#' \code{\link{deterministic.Q.function_template}}
#' 
#' \code{SL.library} may be a character vector of libraries (or \code{NULL} or
#' '\code{default}'), in which case these libraries are used to estimate both
#' \eqn{Q} and \eqn{g} OR a list with two components, \code{Q} and \code{g},
#' where each is a character vector of libraries (or \code{NULL} or
#' '\code{default}').  \code{NULL} indicates \link{glm} should be called
#' instead of \code{\link[SuperLearner:SuperLearner]{SuperLearner}} If
#' \code{SL.library} is the string '\code{default}', \code{SL.library} is set
#' to \code{list("SL.glm", "SL.stepAIC", "SL.bayesglm", c("SL.glm",
#' "screen.corP"), c("SL.step", "screen.corP"), c("SL.step.forward",
#' "screen.corP"), c("SL.stepAIC", "screen.corP"), c("SL.step.interaction",
#' "screen.corP"), c("SL.bayesglm", "screen.corP")}.  Note that the default set
#' of libraries consists of main terms models. It may be advisable to include
#' squared terms, interaction terms, etc in \code{data} or include libraries
#' that consider non-linear terms.
#' 
#' The print method for \code{ltmle} objects only prints the tmle estimates.
#' 
#' @aliases ltmle ltmleMSM
#' @param data data frame following the time-ordering of the nodes. See
#' 'Details'.
#' @param Anodes column names or indicies in \code{data} of treatment nodes
#' @param Cnodes column names or indicies in \code{data} of censoring nodes
#' @param Lnodes column names or indicies in \code{data} of time-dependent
#' covariate nodes
#' @param Ynodes column names or indicies in \code{data} of outcome nodes
#' @param survivalOutcome If \code{TRUE}, then Y nodes are indicators of an
#' event, and if Y at some time point is 1, then all following should be 1.
#' Required to be \code{TRUE} or \code{FALSE} if outcomes are binary and there
#' are multiple Ynodes.
#' @param Qform character vector of regression formulas for \eqn{Q}. See
#' 'Details'.
#' @param gform character vector of regression formulas for \eqn{g} or a
#' matrix/array of prob(A=1). See 'Details'.
#' @param abar binary vector (numAnodes x 1) or matrix (n x numAnodes) of
#' counterfactual treatment or a list of length 2. See 'Details'.
#' @param rule a function to be applied to each row (a named vector) of
#' \code{data} that returns a numeric vector of length numAnodes or a list of
#' length 2. See 'Details'.
#' @param gbounds lower and upper bounds on estimated cumulative probabilities
#' for g-factors. Vector of length 2, order unimportant.
#' @param Yrange NULL or a numerical vector where the min and max of
#' \code{Yrange} specify the range of all Y nodes. See 'Details'.
#' @param deterministic.g.function optional information on A and C nodes that
#' are given deterministically. See 'Details'. Default \code{NULL} indicates no
#' deterministic links.
#' @param stratify if \code{TRUE} stratify on following \code{abar} when
#' estimating Q and g. If \code{FALSE}, pool over \code{abar}.
#' @param SL.library optional character vector of libraries to pass to
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. \code{NULL} indicates
#' \link{glm} should be called instead of
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. '\code{default}'
#' indicates a standard set of libraries. May be separately specified for
#' \eqn{Q} and \eqn{g}. See 'Details'.
#' @param estimate.time if \code{TRUE}, run an initial estimate using only 50
#' observations and use this to print a very rough estimate of the total time
#' to completion. No action if there are fewer than 50 observations.
#' @param gcomp if \code{TRUE}, run the maximum likelihood based G-computation
#' estimate \emph{instead} of TMLE
#' @param regimes binary array: n x numAnodes x numRegimes of counterfactual
#' treatment or a list of 'rule' functions
#' @param working.msm character formula for the working marginal structural
#' model
#' @param summary.measures array: num.regimes x num.summary.measures x
#' num.final.Ynodes - measures summarizing the regimes that will be used on the
#' right hand side of \code{working.msm} (baseline covariates may also be used
#' in the right hand side of \code{working.msm} and do not need to be included
#' in \code{summary.measures})
#' @param final.Ynodes vector subset of Ynodes - used in MSM to pool over a set
#' of outcome nodes
#' @param msm.weights projection weights for the working MSM. If "empirical",
#' weight by empirical proportions of rows matching each regime for each
#' final.Ynode, with duplicate regimes given zero weight. If \code{NULL}, no
#' weights. Or an array of user-supplied weights with dimensions c(n,
#' num.regimes, num.final.Ynodes) or c(num.regimes, num.final.Ynodes).
#' @param iptw.only by default (\code{iptw.only = FALSE}), both TMLE and IPTW
#' are run in \code{ltmle} and \code{ltmleMSM}. If \code{iptw.only = TRUE},
#' only IPTW is run, which is faster.
#' @param deterministic.Q.function optional information on Q given
#' deterministically. See 'Details'. Default \code{NULL} indicates no
#' deterministic links.
#' @param memoize If \code{TRUE}, glm regressions will be memoized. It is
#' recommended to leave this as \code{TRUE} (the default), especially if there
#' are multiple \code{final.Ynodes}, because the code is not written as
#' efficiently as it should be and will end up repeating the same glm call.
#' Will be fixed in a future release.
#' @param observation.weights observation (sampling) weights. Vector of length
#' n. If \code{NULL}, assumed to be all 1.
#' @param IC.variance.only If \code{FALSE}, compute both the robust variance
#' estimate using TMLE and the influence curve based variance estimate (use the
#' larger of the two). If \code{TRUE}, only compute the influence curve based
#' variance estimate, which is faster, but may be substantially
#' anti-conservative if there are positivity violations or rare outcomes.  
#' IC.variance.only=FALSE is not yet available with non-binary outcomes, 
#' gcomp=TRUE, stratify=TRUE, deterministic.Q.function, or numeric gform.
#' 
#' @return \code{ltmle} returns an object of class "\code{ltmle}" (unless
#' \code{abar} or \code{rule} is a list, in which case it returns an object of
#' class \code{ltmleSummaryMeasures}, which has the same components as
#' \code{ltmleMSM}.) The function \code{\link{summary}} (i.e.
#' \code{\link{summary.ltmle}}) can be used to obtain or print a summary of the
#' results. An object of class "\code{ltmle}" is a list containing the
#' following components: \item{estimates}{a named vector of length 4 with
#' elements, each an estimate of \eqn{E[Y_{bar{a}}]}: \itemize{ \item
#' \code{tmle} - Targeted Maximum Likelihood Estimate [NULL if \code{gcomp} is
#' \code{TRUE}] \item \code{iptw} - Inverse Probability of Treatment/Censoring
#' Weighted estimate \item \code{gcomp} - maximum likelihood based
#' G-computation estimate [NULL if \code{gcomp} is \code{FALSE}] } }
#' \item{IC}{a list with the following components of Influence Curve values}
#' \itemize{ \item \code{tmle} - vector of influence curve values for Targeted
#' Maximum Likelihood Estimate [NULL if \code{gcomp} is \code{TRUE}] \item
#' \code{iptw} - vector of influence curve values for Inverse Probability of
#' Treatment/Censoring Weighted estimate \item \code{gcomp} - vector of
#' influence curve values for Targeted Maximum Likelihood Estimate without
#' updating [NULL if \code{gcomp} is \code{FALSE}] } \item{cum.g}{cumulative g,
#' after bounding: for ltmle, n x numACnodes, for ltmleMSM, n x numACnodes x
#' num.regimes} \item{cum.g.unbounded}{cumulative g, before bounding: for
#' ltmle, n x numACnodes, for ltmleMSM, n x numACnodes x num.regimes}
#' \item{call}{the matched call} \item{gcomp}{the \code{gcomp} input}
#' \item{formulas}{a \code{list} with elements \code{Qform} and \code{gform}}
#' \item{fit}{a list with the following components} \itemize{ \item \code{g} -
#' list of length numACnodes - \code{glm} or \code{SuperLearner} return objects
#' from fitting g regressions \item \code{Q} - list of length numLYnodes -
#' \code{glm} or \code{SuperLearner} return objects from fitting Q regressions
#' \item \code{Qstar} - list of length numLYnodes - \code{glm} (or numerical
#' optimization if \code{glm} fails to solve the score equation) return objects
#' from updating the Q fit }
#' 
#' \code{ltmleMSM} returns an object of class "\code{ltmleMSM}" The function
#' \code{\link{summary}} (i.e. \code{\link{summary.ltmleMSM}}) can be used to
#' obtain or print a summary of the results. An object of class
#' "\code{ltmleMSM}" is a list containing the following components:
#' \item{beta}{parameter estimates for working.msm using TMLE (GCOMP if
#' \code{gcomp} input is \code{TRUE})} \item{beta.iptw}{parameter estimates for
#' working.msm using IPTW} \item{IC}{matrix, n x numBetas - influence curve
#' values for TMLE (without updating if \code{gcomp} input is \code{TRUE})}
#' \item{IC.iptw}{matrix, n x numBetas - influence curve values for IPTW}
#' \item{msm}{object of class glm - the result of fitting the working.msm}
#' \item{cum.g}{array, n x numACnodes x numRegimes - cumulative g, after
#' bounding} \item{cum.g.unbounded}{array, n x numACnodes x numRegimes -
#' cumulative g, before bounding} \item{call}{the matched call}
#' \item{gcomp}{the \code{gcomp} input} \item{formulas}{a \code{list} with
#' elements \code{Qform} and \code{gform}} \item{fit}{a list with the following
#' components} \itemize{ \item \code{g} - list of length numRegimes of list of
#' length numACnodes - \code{glm} or \code{SuperLearner} return objects from
#' fitting g regressions \item \code{Q} - list of length numLYnodes -
#' \code{glm} or \code{SuperLearner} return objects from fitting Q regressions
#' \item \code{Qstar} - list of length numLYnodes - \code{glm} (or numerical
#' optimization if \code{glm} fails to solve the score equation) return objects
#' from updating the Q fit }
#' @author Joshua Schwab \email{joshuaschwab@@yahoo.com}, Samuel Lendle, Maya
#' Petersen, and Mark van der Laan
#' @seealso \code{\link{summary.ltmle}}, \code{\link{summary.ltmleMSM}},
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}},
#' \code{\link{deterministic.g.function_template}},
#' \code{\link{deterministic.Q.function_template}}
#' @references Lendle, Samuel, Schwab, Joshua, Petersen, Maya and and van der
#' Laan, Mark J "ltmle: An R Package Implementing Targeted Minimum Loss-based
#' Estimation for Longitudinal Data", Forthcoming
#' 
#' Petersen, Maya, Schwab, Joshua and van der Laan, Mark J, "Targeted Maximum
#' Likelihood Estimation of Marginal Structural Working Models for Dynamic
#' Treatments Time-Dependent Outcomes", Forthcoming
#' 
#' van der Laan, Mark J. and Gruber, Susan, "Targeted Minimum Loss Based
#' Estimation of an Intervention Specific Mean Outcome" (August 2011). U.C.
#' Berkeley Division of Biostatistics Working Paper Series. Working Paper 290.
#' \url{http://biostats.bepress.com/ucbbiostat/paper290}
#' 
#' van der Laan, Mark J. and Rose, Sherri, "Targeted Learning: Causal Inference
#' for Observational and Experimental Data" New York: Springer, 2011.
#' @examples
#' 
#' rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
#'                  
#' # Example 1: Single time point example. True value of E[Y_1] (expected value of
#' #   Y setting A to 1) is approximately 0.5939.
#' set.seed(2)
#' n <- 1000
#' W1 <- rnorm(n)
#' W2 <- rbinom(n, size=1, prob=0.3)   
#' W3 <- rnorm(n)
#' A <- rexpit(-1 + 2 * W1^2)
#' Y <- rexpit(-0.5 + 2 * W1^2 + 0.5 * W2 - 0.5 * A + 0.2 * W3 * A 
#'        - 1.1 * W3 + 0.2 * rnorm(n))
#' 
#' data <- data.frame(W1, W2, W3, A, Y)
#' 
#' 
#' \donttest{ #This takes about 4 seconds to run
#' library(SuperLearner)
#' 
#' #SuperLearner semiparametric estimation using all parents as regressors 
#' result1 <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", abar=1, 
#'   SL.library=c("SL.glm", "SL.step", "SL.mean"))
#' summary(result1)
#' summary(result1, estimator="iptw")
#' 
#' #SuperLearner semiparametric estimation using (incorrectly) specified regressors
#' #note: The functional form for Qform and gform is unimportant if 
#' # using SuperLearner - see 'Details'
#' result1a <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
#'  Qform=c(Y="Q.kplus1 ~ W1 + W3 + A"), gform="A ~ W1", abar=1, 
#'  SL.library=c("SL.glm", "SL.step", "SL.mean"))
#' summary(result1a)
#' }
#' 
#' #glm using correctly specified Qform and gform
#' result.abar1 <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
#'  Qform=c(Y="Q.kplus1 ~ I(W1^2) + W2 + W3*A"), gform="A ~ I(W1^2)", 
#'  abar=1, SL.library=NULL)
#' 
#' \donttest{ #This takes about 18 seconds to run
#' #Get summary measures (additive treatment effect, odds ratio, relative risk) 
#' #  for abar=1 vs abar=0
#' result.compare <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
#'                       Qform=c(Y="Q.kplus1 ~ I(W1^2) + W2 + W3*A"), gform="A ~ I(W1^2)", 
#'                       abar=list(1, 0), SL.library=NULL)
#' summary(result.compare)
#' 
#' 
#' # Example 2: Longitudinal example. Includes informative censoring and treatment. 
#' # Time ordering of data is W, C1, L1, A1, Y1, C2, L2, A2, Y2
#' # True value of E[Y_(1,1,1,1)] (expected value of Y setting C1, A1, C2, A2 all to 1)
#' #  is approximately 0.413.
#' # A1 is known to always be 1 if L1 < -2, and is 1 with probability 0.1 if L1 > 2 
#' # A2 is known to always be 1 if A1 is 1 
#' # We incorporate this knowledge using deterministic.g.function
#' 
#' # Generate data:
#' set.seed(2)
#' ua <- rep(TRUE, n)   #ua = uncensored and alive
#' L1 <- A1 <- Y1 <- C2.binary <- L2 <- A2 <- Y2 <- as.numeric(rep(NA, n))
#' W <- rnorm(n)
#' C1 <- BinaryToCensoring(is.uncensored=rexpit(2 + W))
#' ua <- ua & C1 == "uncensored"
#' L1[ua] <- rnorm(n)[ua] + W[ua]
#' A1[ua] <- rexpit(L1[ua])
#' A1[ua & L1 < -2] <- 1
#' A1[ua & L1 >  2] <- rbinom(n, size=1, prob=0.1)[ua & L1 >  2]
#' Y1[ua] <- rexpit((W + L1 - A1)[ua])
#' ua <- ua & !Y1
#' C2.binary[ua] <- rexpit((1 + 0.7 * L1 - A1)[ua])
#' C2 <- BinaryToCensoring(is.uncensored=C2.binary)
#' ua <- ua & C2 == "uncensored"
#' L2[ua] <- (0.5 * L1 - 0.9 * A1 + rnorm(n))[ua]
#' A2[ua] <- rexpit((0.5 * L1 + 0.8 * L2)[ua]) | A1[ua]
#' Y2[ua] <- rexpit((0.7 * L1 + L2 - 0.8 * A1 - A2)[ua])
#' Y2[Y1 == 1] <- 1  # if a patient dies at time 1, record death at time 2 as well
#' data <- data.frame(W, C1, L1, A1, Y1, C2, L2, A2, Y2)
#' 
#' deterministic.g.function <- function(data, current.node, nodes) {
#'   if (names(data)[current.node] == "A1") {
#'     det <- (data$L1 < -2 | data$L1 > 2) & !is.na(data$L1)
#'     prob1 <- ((data$L1 < -2) * 1 + (data$L1 > 2) * 0.1)[det]
#'   } else if (names(data)[current.node] == "A2") {
#'     det <- data$A1 == 1 & !is.na(data$A1)
#'     prob1 <- 1
#'   } else if (names(data[current.node]) %in% c("C1", "C2")){
#'     return(NULL)  #this returns the default of no deterministic links 
#'     #note that it is not necessary to specify that prior censoring indicates future censoring
#'   } else {
#'     stop("unexpected current.node")
#'   }
#'   return(list(is.deterministic=det, prob1=prob1))  
#' }
#' 
#' result2 <- ltmle(data, Anodes=c("A1","A2"), Cnodes=c("C1", "C2"), 
#'                 Lnodes=c("L1", "L2"), Ynodes=c("Y1", "Y2"), abar=c(1, 1), 
#'                 deterministic.g.function=deterministic.g.function, survivalOutcome=TRUE)
#' summary(result2) 
#'  
#' # Example 3: Dynamic treatment, observation weights
#' # W -> A1 -> L -> A2 -> Y
#' # Treatment regime of interest is: Always treat at time 1 (A1 = 1), 
#' #   treat at at time 2 (A2 = 1), iff L > 0
#' # Weight by pmax(W + 2, 0)
#' 
#' set.seed(2)
#' n <- 1000
#' W <- rnorm(n)
#' A1 <- rexpit(W)
#' L <- 0.3 * W + 0.2 * A1 + rnorm(n)
#' A2 <- rexpit(W + A1 + L)
#' Y <- rexpit(W - 0.6 * A1 + L - 0.8 * A2)
#' data <- data.frame(W, A1, L, A2, Y)
#' 
#' abar <- matrix(nrow=n, ncol=2)
#' abar[, 1] <- 1
#' abar[, 2] <- L > 0
#' 
#' result3 <- ltmle(data, Anodes=c("A1", "A2"), Lnodes="L", Ynodes="Y", 
#'   survivalOutcome=TRUE, abar=abar, observation.weights = pmax(W + 2, 0))
#' summary(result3)
#' 
#' # Example 3.1: The regime can also be specified as a rule function
#' 
#' rule <- function(row) c(1, row["L"] > 0)
#' 
#' result.rule <- ltmle(data, Anodes=c("A1", "A2"), Lnodes="L", Ynodes="Y", 
#'   survivalOutcome=TRUE, rule=rule, observation.weights = pmax(W + 2, 0))
#' # This should be the same as the above result
#' summary(result.rule)
#' 
#' # Example 4: Deterministic Q function
#' # W -> A1 -> Y1 -> L2 -> A2 -> Y2
#' set.seed(2)
#' n <- 200
#' L2 <- A2 <- Y2 <- as.numeric(rep(NA, n))
#' W <- rnorm(n)
#' A1 <- rexpit(W)
#' Y1 <- rexpit(W - A1)
#' alive <- Y1 == 0
#' L2[alive] <- (0.5 * W - 0.9 * A1 + rnorm(n))[alive]
#' completed.study <- alive & L2 > 0
#' 
#' #Specify that Q is deterministically 0 when L2 is in the history of the 
#' # current Q regression and L2 > 0
#' #Note 1: det.Q.fun doesn't condition on called.from.estimate.g so g will also be set 
#' #        deterministically after L2 > 0 
#' #Note 2: It is not necessary to specify that Q is deterministically 1 if Y1 is 1; this is automatic
#' det.Q.fun.4a <- function(data, current.node, nodes, called.from.estimate.g) {
#'   L2.index <- which(names(data) == "L2")
#'   stopifnot(length(L2.index) == 1)
#'   L2.in.history <- L2.index < current.node
#'   if (! L2.in.history) return(NULL)
#'   
#'   is.deterministic <- data$L2 > 0 & !is.na(data$L2)
#'   return(list(is.deterministic=is.deterministic, Q.value=0))
#' }
#' 
#' #patients don't change treatment after leaving study; leave their A2 as NA
#' A2[alive & !completed.study] <- rexpit((0.5 * W + 0.8 * L2)[alive & !completed.study])
#' 
#' Y2[alive & !completed.study] <- rexpit((L2 - 0.8 * A1 - A2)[alive & !completed.study])
#' Y2[alive & completed.study] <- 0
#' Y2[!alive] <- 1  # if a patient dies at time 1, record death at time 2 as well
#' data <- data.frame(W, A1, Y1, L2, A2, Y2)
#' 
#' result4a <- ltmle(data, Anodes=c("A1","A2"), Lnodes="L2", Ynodes=c("Y1", "Y2"), abar=c(1, 1), 
#'   SL.library=NULL, estimate.time=FALSE, deterministic.Q.function=det.Q.fun.4a, survivalOutcome=TRUE)
#' #note: You will get the same result if you pass Lnodes=NULL (see next example)
#' summary(result4a)
#' 
#' #In this variation, suppose that treatment can still change after a patient leaves the study
#' 
#' det.Q.fun.4b <- function(data, current.node, nodes, called.from.estimate.g) {
#'   #there is no deterministic information when calculating g - treatment may still change
#'   if (called.from.estimate.g) return(NULL)  
#'   
#'   L2.index <- which(names(data) == "L2")
#'   stopifnot(length(L2.index) == 1)
#'   L2.in.history <- L2.index < current.node
#'   if (! L2.in.history) return(NULL)
#'   
#'   is.deterministic <- data$L2 > 0 & !is.na(data$L2)
#'   return(list(is.deterministic=is.deterministic, Q.value=0))
#' }
#' 
#' A2[alive] <- rexpit((0.5 * W + 0.8 * L2)[alive])  #patients can change treatment after leaving study
#' Y2[alive & !completed.study] <- rexpit((L2 - 0.8 * A1 - A2)[alive & !completed.study])
#' Y2[alive & completed.study] <- 0
#' Y2[!alive] <- 1  # if a patient dies at time 1, record death at time 2 as well
#' data <- data.frame(W, A1, Y1, L2, A2, Y2)
#' 
#' result4b <- ltmle(data, Anodes=c("A1","A2"), Lnodes="L2", Ynodes=c("Y1", "Y2"), abar=c(1, 1), 
#'  SL.library=NULL, estimate.time=FALSE, deterministic.Q.function=det.Q.fun.4b, survivalOutcome=TRUE)
#' summary(result4b)
#' 
#' # Example 5: Multiple time-dependent covariates and treatments at each time point, 
#' #            continuous Y values
#' # age -> gender -> A1 -> L1a -> L1b -> Y1 -> A2 -> L2a -> L2b -> Y2
#' set.seed(2)
#' n <- 100
#' age <- rbinom(n, 1, 0.5)
#' gender <- rbinom(n, 1, 0.5)
#' A1 <- rexpit(age + gender)
#' L1a <- 2*age - 3*gender + 2*A1 + rnorm(n)
#' L1b <- rexpit(age + 1.5*gender - A1)
#' Y1 <- plogis(age - gender + L1a + 0.7*L1b + A1 + rnorm(n))
#' A2 <- rexpit(age + gender + A1 - L1a - L1b)
#' L2a <- 2*age - 3*gender + 2*A1 + A2 + rnorm(n)
#' L2b <- rexpit(age + 1.5*gender - A1 - A2)
#' Y2 <- plogis(age - gender + L1a + L1b + A1 + 1.8*A2 + rnorm(n))
#' data <- data.frame(age, gender, A1, L1a, L1b, Y1, A2, L2a, L2b, Y2)
#' 
#' #Note that gform is not correctly specified in these examples.
#' 
#' #Also show some different ways of specifying the nodes:
#' 
#' result5a <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), 
#'  Ynodes=grep("^Y", names(data)), abar=c(1, 0), SL.library=NULL, estimate.time=FALSE, 
#'  survivalOutcome=FALSE, gform=c("A1 ~ gender", "A2 ~ age")) 
#' summary(result5a)
#' 
#' #Usually you would specify a Qform for all of the Lnodes and Ynodes but in this case 
#' # L1a, L1b, Y1 is a "block" of L/Y nodes not separated by Anodes or Cnodes (the same is true for 
#' # L2a, L2b, Y2). Only one regression is required at the first L/Y node in a block. You can pass 
#' # regression formulas for the other L/Y nodes, but they'll be ignored.
#' result5b <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), 
#'  Ynodes=grep("^Y", names(data)), abar=c(1, 0), estimate.time=FALSE, survivalOutcome=FALSE, 
#'  gform=c("A1 ~ gender", "A2 ~ age"), Qform=c(L1a="Q.kplus1 ~ 1", L2a="Q.kplus1 ~ 1"))
#' summary(result5b)
#' 
#' 
#' #Gives the same result but prints a message saying some regression formulas will be dropped:
#' result5c <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), 
#'  Ynodes=grep("^Y", names(data)), abar=c(1, 0), estimate.time=FALSE, survivalOutcome=FALSE, 
#'  gform=c("A1 ~ gender", "A2 ~ age"), Qform=c(L1a="Q.kplus1 ~ 1", L1b="Q.klus1~A1", 
#'  Y1="Q.kplus1~L1a", L2a="Q.kplus1 ~ 1", L2b="Q.klus1~A1", Y2="Q.kplus1~A2 + gender"))
#' 
#' summary(result5c)
#' 
#' 
#' #If there were a Anode or Cnode between L1b and Y1, Y1 would also need a Q regression formula
#' 
#' 
#' # Example 6: MSM
#' # Given data over 3 time points where A switches to 1 once and then stays 1. We want to know
#' # how death varies as a function of gender, time and an indicator of whether a patient's 
#' # intended regime was to switch before time.
#' # Note that working.msm includes time and switch.time, which are columns of 
#' # summary.measures; working.msm also includes male, which is ok because it is a baseline
#' # covariate (it comes before any A/C/L/Y nodes).
#' data(sampleDataForLtmleMSM)
#' Anodes <- grep("^A", names(sampleDataForLtmleMSM$data))
#' Lnodes <- c("CD4_1", "CD4_2")
#' Ynodes <- grep("^Y", names(sampleDataForLtmleMSM$data))
#' msm.weights <- matrix(1:12, nrow=4, ncol=3) #just an example (can also use a 200x3x4 array), 
#'                                             #or NULL (for no weights), or "empirical" (the default)
#' 
#' result6 <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, 
#'                    survivalOutcome=TRUE,
#'                    regimes=sampleDataForLtmleMSM$regimes, 
#'                    summary.measures=sampleDataForLtmleMSM$summary.measures, final.Ynodes=Ynodes, 
#'                    working.msm="Y ~ male + time + I(pmax(time - switch.time, 0))", 
#'                    msm.weights=msm.weights, estimate.time=FALSE)
#' print(summary(result6))
#' 
#' 
#' # Example 6.1: regimes can also be specified as a list of rule functions
#' 
#' regimesList <- list(function(row) c(1,1,1),
#'                      function(row) c(0,1,1),
#'                      function(row) c(0,0,1),
#'                      function(row) c(0,0,0))
#' result.regList <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, 
#'                    survivalOutcome=TRUE, regimes=regimesList, 
#'                    summary.measures=sampleDataForLtmleMSM$summary.measures, final.Ynodes=Ynodes, 
#'                    working.msm="Y ~ male + time + I(pmax(time - switch.time, 0))", 
#'                    msm.weights=msm.weights, estimate.time=FALSE)
#' # This should be the same as the above result
#' print(summary(result.regList))         
#' 
#' 
#' # Example 7: variance estimation
#' # A simple point treatment problem W, A, Y. But there is a positivity problem - 
#' # for small values of W, Prob(A = 1) is very small.
#' # The true parameter value, E[Y_1] is approximately 0.697
#' # The true TMLE standard deviation is approximately 0.064, 
#' # the true IPTW standard deviation is approximately 0.058.
#' set.seed(2)
#' n <- 1000
#' W <- rnorm(n)
#' A <- rexpit(8 * W)
#' Y <- rexpit(W + A)
#' r1 <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE)
#' r2 <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE, 
#'  IC.variance.only=TRUE)
#' print(summary(r1))
#' print(summary(r2))
#' print(summary(r1, "iptw"))
#' print(summary(r2, "iptw")) #the same - IC.variance.only only affects TMLE
#' }
#' 
#' @export ltmle
ltmle <- function(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, survivalOutcome=NULL, Qform=NULL, gform=NULL, 
                  abar, rule=NULL, gbounds=c(0.01, 1), Yrange=NULL, deterministic.g.function=NULL, stratify=FALSE, 
                  SL.library=NULL, estimate.time=TRUE, gcomp=FALSE, 
                  iptw.only=FALSE, deterministic.Q.function=NULL, IC.variance.only=FALSE, observation.weights=NULL) {
  msm.inputs <- GetMSMInputsForLtmle(data, abar, rule, gform)
  inputs <- CreateInputs(data=data, Anodes=Anodes, Cnodes=Cnodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=survivalOutcome, Qform=Qform, gform=msm.inputs$gform, Yrange=Yrange, gbounds=gbounds, deterministic.g.function=deterministic.g.function, SL.library=SL.library, regimes=msm.inputs$regimes, working.msm=msm.inputs$working.msm, summary.measures=msm.inputs$summary.measures, final.Ynodes=msm.inputs$final.Ynodes, stratify=stratify, msm.weights=msm.inputs$msm.weights, estimate.time=estimate.time, gcomp=gcomp, iptw.only=iptw.only, deterministic.Q.function=deterministic.Q.function, IC.variance.only=IC.variance.only, observation.weights=observation.weights) 
  result <- LtmleFromInputs(inputs)
  result$call <- match.call()
  return(result)
}
# General code flow:
#  ltmle -> CreateInputs -> LtmleFromInputs -> LtmleMSMFromInputs -> ...
#  ltmleMSM -> CreateInputs -> LtmleMSMFromInputs -> ...

RegimesFromAbar <- function(data, abar, rule) {
  if (!is.null(rule)) {
    if (!(missing(abar) || is.null(abar))) stop("'abar' should not be specified when using a 'rule' function")
    abar <- t(apply(data, 1, rule))
  }
  if (is.vector(abar)) {
    abar <- matrix(rep(abar, each=nrow(data)), nrow=nrow(data))
  } else if (is.null(abar)) {
    abar <- matrix(nrow=nrow(data), ncol=0)
  }
  regimes <- abar
  dim(regimes) <- c(nrow(regimes), ncol(regimes), 1)
  return(regimes)
}

# ltmle is a special case of ltmleMSM - get the arguments used by ltmleMSM for the special case
GetMSMInputsForLtmle <- function(data, abar, rule, gform) {
  if ((!missing(abar) && is.list(abar)) || is.list(rule)) {
    if (is.list(rule)) {
      if (length(rule) != 2) stop("If rule is a list, it must be of length 2")
      regimes1 <- RegimesFromAbar(data, abar, rule[[1]])
      regimes0 <- RegimesFromAbar(data, abar, rule[[2]])
    } else {
      if (length(abar) != 2) stop("If abar is a list, it must be of length 2")
      regimes1 <- RegimesFromAbar(data, abar[[1]], rule)
      regimes0 <- RegimesFromAbar(data, abar[[2]], rule)
    }
    if (ncol(regimes1) != ncol(regimes0)) stop("If abar or rule is a list, both elements must give a matrix with the same number of columns")
    if (nrow(regimes1) != nrow(regimes0)) stop("If abar or rule is a list, both elements must give a matrix with the same number of rows")
    regimes <- c(regimes1, regimes0)
    dim(regimes) <- c(nrow(regimes1), ncol(regimes1), 2)
    summary.measures <- array(1:0, dim=c(2, 1, 1))
    colnames(summary.measures) <- "A"
    working.msm <- "Y ~ A"
    msm.weights <- matrix(1, nrow=2, ncol=1)
  } else {
    regimes <- RegimesFromAbar(data, abar, rule)
    working.msm <- "Y ~ 1"
    msm.weights <- matrix(1, nrow=1, ncol=1)
    summary.measures <- array(dim=c(1, 0, 1))
  }
  if (is.numeric(gform)) {
    stopifnot(is.matrix(gform))
    dim(gform) <- c(nrow(gform), ncol(gform), 1)
  }
  msm.inputs <- list(regimes=regimes, working.msm=working.msm, summary.measures=summary.measures, gform=gform, final.Ynodes=NULL, msm.weights=msm.weights)
  return(msm.inputs)
}

# run ltmle from the ltmleInputs object
LtmleFromInputs <- function(inputs) {
  msm.result <- LtmleMSMFromInputs(inputs)
  num.regimes <- dim(inputs$regimes)[3]
  stopifnot(num.regimes %in% 1:2)
  if (num.regimes == 2) {
    class(msm.result) <- "ltmleEffectMeasures"
    return(msm.result)
  }
  names(msm.result$beta.iptw) <- names(msm.result$beta) <- NULL
  iptw <- plogis(msm.result$beta.iptw)
  iptw.list <- list(iptw.estimate=iptw, iptw.IC=iptw*(1-iptw)*msm.result$IC.iptw[, 1])
  
  r <- list()
  if (inputs$iptw.only) {
    tmle <- NA
    tmle.IC <- rep(NA, nrow(inputs$data))
  } else {
    tmle <- plogis(msm.result$beta)
    tmle.IC <- msm.result$IC[, 1] #only one regime
  }
  r$estimates <- c(tmle=tmle, iptw=iptw.list$iptw.estimate)
  r$IC <- list(tmle=tmle.IC * tmle * (1 - tmle), iptw=iptw.list$iptw.IC)
  if (!is.null(msm.result$variance.estimate)) {
    stopifnot(length(msm.result$variance.estimate) == 1)
    r$variance.estimate <- msm.result$variance.estimate[1] * (tmle * (1 - tmle))^2 
  }
  
  if (inputs$gcomp) {
    names(r$estimates)[1] <- names(r$IC)[1] <- "gcomp"
  }
  
  r$cum.g <- AsMatrix(msm.result$cum.g[, , 1]) #only one regime
  r$cum.g.unbounded <- AsMatrix(msm.result$cum.g.unbounded[, , 1]) #only one regime
  r$gcomp <- inputs$gcomp
  r$fit <- msm.result$fit
  r$fit$g <- r$fit$g[[1]]  #only one regime
  r$fit$Q <- r$fit$Q[[1]]  #only one regime 
  r$Qstar <- msm.result$Qstar[, 1, 1] #1 regime, 1 final.Ynode
  
  r$formulas <- msm.result$formulas
  r$binaryOutcome <- msm.result$binaryOutcome
  r$transformOutcome <- msm.result$transformOutcome==TRUE #Want to store transformOutcome flag without attributes
  
  if (msm.result$transformOutcome) {
    Yrange <- attr(msm.result$transformOutcome, "Yrange")
    #back transform estimate and IC
    r$estimates <- r$estimates*diff(Yrange) + min(Yrange)  
    r$IC <- lapply(r$IC, function (IC) IC * diff(Yrange))
    r$variance.estimate <- r$variance.estimate * (diff(Yrange))^2 
  }
  class(r) <- "ltmle"
  return(r)
}

#' @describeIn ltmle Longitudinal Targeted Maximum Likelihood Estimation for a Marginal Structural Model
#' @export 
ltmleMSM <- function(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, survivalOutcome=NULL, Qform=NULL, gform=NULL, gbounds=c(0.01, 1), Yrange=NULL, deterministic.g.function=NULL, SL.library=NULL, regimes, working.msm, summary.measures, final.Ynodes=NULL, stratify=FALSE, msm.weights="empirical", estimate.time=TRUE, gcomp=FALSE, iptw.only=FALSE, deterministic.Q.function=NULL, memoize=TRUE, IC.variance.only=FALSE, observation.weights=NULL) {
  if (memoize && requireNamespace("memoise")) {
    glm.ltmle.memoized <- memoise::memoize(glm.ltmle)
  }
  
  inputs <- CreateInputs(data, Anodes, Cnodes, Lnodes, Ynodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, estimate.time, gcomp, iptw.only, deterministic.Q.function, IC.variance.only, observation.weights)
  result <- LtmleMSMFromInputs(inputs)
  result$call <- match.call()
  class(result) <- "ltmleMSM"
  return(result) 
}

# run ltmleMSM from ltmleInputs object
LtmleMSMFromInputs <- function(inputs) {  
  if (inputs$estimate.time) EstimateTime(inputs)
  result <- MainCalcs(inputs)
  result$gcomp <- inputs$gcomp
  result$formulas <- list(Qform=inputs$Qform, gform=inputs$gform)
  result$binaryOutcome <- inputs$binaryOutcome
  result$transformOutcome <- inputs$transformOutcome
  result$survivalOutcome <- inputs$survivalOutcome
  return(result)
}

# create the ltmleInputs object used by many other functions - fills in defaults and does error checking
CreateInputs <- function(data, Anodes, Cnodes, Lnodes, Ynodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, estimate.time, gcomp, iptw.only, deterministic.Q.function, IC.variance.only, observation.weights) {
  if (is.list(regimes)) {
    if (!all(do.call(c, lapply(regimes, is.function)))) stop("If 'regimes' is a list, then all elements should be functions.")
    regimes <- aperm(simplify2array(lapply(regimes, function(rule) apply(data, 1, rule)), higher=TRUE), c(2, 1, 3)) 
  }
  if (!(is.null(regimes) || dim(regimes) != 3)) {
    stop("regimes must be an array with 3 dimensions (unless Anodes is NULL, in which case regimes can be NULL)")
  }
  if (is.null(regimes) || dim(regimes)[3]==0) {
    if (length(Anodes) != 0) {
      stop("regimes must not be NULL (or have dim(regimes)[3]==0) unless Anodes is also NULL")
    }
    regimes <- array(dim=c(nrow(data), 0, 1))
  }
  if (is.logical(regimes)) {
    regimes <- regimes * 1
    message("abar or regimes was passed as logical and was converted to numeric")
  }
  nodes <- CreateNodes(data, Anodes, Cnodes, Lnodes, Ynodes)
  Qform <- CreateLYNodes(data, nodes, check.Qform=TRUE, Qform=Qform)$Qform
  data <- ConvertCensoringNodes(data, Cnodes, has.deterministic.functions=!is.null(deterministic.g.function) && is.null(deterministic.Q.function))
  if (is.null(final.Ynodes)) {
    final.Ynodes <- max(nodes$Y)
  } else {
    final.Ynodes <- NodeToIndex(data, final.Ynodes)
  }
  
  #Using get to avoid the "no visible binding for global variable" note in R CMD check
  if (identical(SL.library, 'default')) SL.library <- get("Default.SL.Library")
  SL.library.Q <- GetLibrary(SL.library, "Q")
  SL.library.g <- GetLibrary(SL.library, "g")
    
  if (is.null(summary.measures)) {
    summary.measures <- matrix(nrow=dim(regimes)[3], ncol=0)
  }
  if (length(dim(summary.measures)) == 2) {
    num.final.Ynodes <- length(final.Ynodes)
    summary.measures <- array(repmat(summary.measures, m=1, n=num.final.Ynodes), dim=c(nrow(summary.measures), ncol(summary.measures), num.final.Ynodes), dimnames=list(rownames(summary.measures), colnames(summary.measures), NULL))
  }
  if (is.null(observation.weights)) observation.weights <- rep(1, nrow(data))
  
  #error checking (also get value for survivalOutcome)
  check.results <- CheckInputs(data, nodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, deterministic.Q.function, observation.weights, gcomp) 
  survivalOutcome <- check.results$survivalOutcome
  
  
  if (!isTRUE(attr(data, "called.from.estimate.variance", exact=TRUE))) { 
    data <- CleanData(data, nodes, deterministic.Q.function, survivalOutcome)
  }
  untransformed.data <- data
  transform.list <- TransformOutcomes(data, nodes, Yrange)
  data <- transform.list$data
  transformOutcome <- transform.list$transformOutcome
  binaryOutcome <- check.results$binaryOutcome
  
  if (is.null(Qform)) Qform <- GetDefaultForm(data, nodes, is.Qform=TRUE, stratify, survivalOutcome, showMessage=TRUE)
  if (is.null(gform)) gform <- GetDefaultForm(data, nodes, is.Qform=FALSE, stratify, survivalOutcome, showMessage=TRUE)
  
  inputs <- list(data=data, untransformed.data=untransformed.data, nodes=nodes, survivalOutcome=survivalOutcome, Qform=Qform, gform=gform, gbounds=gbounds, Yrange=Yrange, deterministic.g.function=deterministic.g.function, SL.library.Q=SL.library.Q, SL.library.g=SL.library.g, regimes=regimes, working.msm=working.msm, summary.measures=summary.measures, final.Ynodes=final.Ynodes, stratify=stratify, msm.weights=msm.weights, estimate.time=estimate.time, gcomp=gcomp, iptw.only=iptw.only, deterministic.Q.function=deterministic.Q.function, binaryOutcome=binaryOutcome, transformOutcome=transformOutcome, IC.variance.only=IC.variance.only, observation.weights=observation.weights)
  class(inputs) <- "ltmleInputs"
  if (!inputs$IC.variance.only && !is.null(VarianceAvailableWarning(inputs))) inputs$IC.variance.only <- TRUE
  return(inputs)
}

# Prevent misspelled argument after $ from returning NULL
`$.ltmleInputs` <- function(x, name) {
  if (! (name %in% names(x))) stop(paste(name, "is not an element of x"))
  value <- x[[name]]
  if (identical(value, "ltmleInputs-NULL")) {
    return(NULL)
  } else {
    return(value)
  }
}

# Prevent misspelled argument after $<- from adding new element
`$<-.ltmleInputs` <- function(x, name, value) {
  if (! (name %in% names(x))) stop(paste(name, "is not an element of x"))
  if (is.null(value)) {
    value <- "ltmleInputs-NULL"
  }
  x[[name]] <- value
  return(x)
}

# Loop over final Ynodes, run main calculations
MainCalcs <- function(inputs) {
  #if (inputs$iptw.only) { #could use something like this to save time, but we want to pass all final ynodes to CalcIPTW below
  #  inputs$final.Ynodes <- inputs$final.Ynodes[length(inputs$final.Ynodes)]
  #}
  # Several functions in the pooled version are only written to accept main terms MSM
  # Ex: If working.msm is "Y ~ X1*X2", convert to "Y ~ -1 + S1 + S1 + S3 + S4" where 
  # S1 is 1 (intercept), S2 is X1, S3 is X2, S4 is X1:X2
  main.terms <- ConvertToMainTerms(inputs$data, inputs$working.msm, inputs$summary.measures, inputs$nodes)
  inputs$working.msm <- main.terms$msm
  combined.summary.measures <- main.terms$summary.measures   
  baseline.column.names <- main.terms$baseline.column.names
  num.final.Ynodes <- length(inputs$final.Ynodes)
  num.betas <- length(main.terms$beta.names)
  #combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes       note: num.measures is summary measures and baseline covariates, converted to main terms
  n <- nrow(inputs$data)
  num.regimes <- dim(inputs$regimes)[3]
  Qstar <- g.ratio <- array(dim=c(n, num.regimes, num.final.Ynodes))
  msm.weights <- GetMsmWeights(inputs) #n x num.regimes x num.final.Ynodes
  new.var.y <- array(dim=c(num.betas, num.betas, num.final.Ynodes))
  IC <- matrix(0, n, num.betas)
  #store IC for each final Ynode, compare var(IC) to sum(var(IC.ynode))
  IC.y <- array(dim=c(n, num.betas, num.final.Ynodes))
  for (j in 1:num.final.Ynodes) {
    fixed.tmle <- FixedTimeTMLE(SubsetInputs(inputs, final.Ynode=inputs$final.Ynodes[j]), drop3(msm.weights[, , j, drop=FALSE]), dropn(combined.summary.measures[, , , j, drop=FALSE], n=4), baseline.column.names)
    IC <- IC + fixed.tmle$IC
    IC.y[, , j] <- fixed.tmle$IC
    Qstar[, , j] <- fixed.tmle$Qstar # n x num.regimes
    new.var.y[, , j] <- fixed.tmle$est.var 
    g.ratio[, , j] <- fixed.tmle$g.ratio
  }
  iptw <- CalcIPTW(inputs$data, inputs$nodes, inputs$working.msm, inputs$regimes, combined.summary.measures, inputs$final.Ynodes, fixed.tmle$cum.g, msm.weights, inputs$observation.weights)
  names(iptw$beta) <- main.terms$beta.names

  if (inputs$iptw.only) {
    beta <- rep(NA, length(iptw$beta))
    fitted.msm <- NULL
    variance.estimate <- NULL
  } else {
    fitted.msm <- FitPooledMSM(inputs$working.msm, Qstar, combined.summary.measures, msm.weights * inputs$observation.weights)
    if (isTRUE(attr(inputs$data, "called.from.estimate.variance", exact=TRUE))) { 
      #variance estimate is not needed, this avoids some warnings
      IC <- matrix(NA, n, num.betas)
      C.old <- matrix(NA, num.betas, num.betas)
    } else {
      CheckForVarianceWarning(inputs, g.ratio)
      IC <- FinalizeIC(IC, combined.summary.measures, Qstar, fitted.msm$m.beta, msm.weights, g.ratio, inputs$observation.weights) #n x num.betas
      C.old <- NormalizeIC(IC, combined.summary.measures, fitted.msm$m.beta, msm.weights, g.ratio = array(1, dim=c(n, num.regimes, num.final.Ynodes)), inputs$observation.weights) #C without using g.ratio (setting g.ratio to 1)
    }
 
    if (inputs$IC.variance.only) {   
      variance.estimate <- NULL
    } else {
      new.var <- matrix(NA, num.betas, num.betas)
      for (i in 1:num.betas) {
        for (j in 1:num.betas) {
          if (num.final.Ynodes > 1) {
            cov.IC <- cov(IC.y[, i, ], IC.y[, j, ])
            diag(cov.IC) <- new.var.y[i, j, ]
            new.var[i, j] <- sum(cov.IC)
          } else {
            new.var[i, j] <- new.var.y[i, j, 1]
          }
        }
      }
      C <- NormalizeIC(IC, combined.summary.measures, fitted.msm$m.beta, msm.weights, g.ratio, inputs$observation.weights)
      variance.estimate <- safe.solve(C) %*% new.var %*% t(safe.solve(C))
    }
    
    IC <- t(safe.solve(C.old, t(IC))) #IC %*% solve(C) 
    beta <- coef(fitted.msm$m)
    names(beta) <- main.terms$beta.names
  }
  return(list(IC=IC, msm=fitted.msm$m, beta=beta, cum.g=fixed.tmle$cum.g, cum.g.unbounded=fixed.tmle$cum.g.unbounded, fit=fixed.tmle$fit, variance.estimate=variance.estimate, beta.iptw=iptw$beta, IC.iptw=iptw$IC, Qstar=Qstar, g.ratio=g.ratio)) #note: only returns cum.g and fit for the last final.Ynode
}


VarianceAvailableWarning <- function(inputs) {
  if (!inputs$binaryOutcome) return("Robust variance estimate is not currently available with non binary outcomes")
  if (!is.null(inputs$deterministic.Q.function)) return("Robust variance estimate is not currently available with deterministic.Q.function")
  if (inputs$gcomp) return("Robust variance estimate is not currently available with gcomp")
  if (inputs$stratify) return("Robust variance estimate is not currently available with stratify=TRUE")
  if (is.numeric(inputs$gform)) return("Robust variance estimate is not currently available with numeric gform")
  return(NULL)
}

CheckForVarianceWarning <- function(inputs, g.ratio) {
  if (inputs$IC.variance.only) {
    positivity <- mean(g.ratio < 1, na.rm=TRUE) > 0.01
    rare.events <- inputs$binaryOutcome && (colMeans(inputs$data[, inputs$final.Ynodes, drop=FALSE], na.rm=TRUE) < 0.03) 
    if (positivity || rare.events) {
      variance.available.warning <- VarianceAvailableWarning(inputs)
      warning.msg <- "Variance estimate is based on influence curve only, which may be significantly anticonservative because your data appears to contain"
      if (positivity) warning.msg <- paste(warning.msg, "positivity violations")
      if (positivity && rare.events) warning.msg <- paste(warning.msg, "and")
      if (rare.events) warning.msg <- paste(warning.msg, "rare events")
      if (is.null(variance.available.warning)) {
        warning.msg <- paste0(warning.msg, ". It is recommended to use IC.variance.only=FALSE to obtain a more robust variance estimate (but run time may be significantly longer).")
      } else {
        warning.msg <- paste0(warning.msg, ". ", variance.available.warning, " but this will be addressed in a future release.")
      }
      warning(warning.msg)
    }
  }
  invisible(NULL)
}

CalcIPTW <- function(data, nodes, working.msm, regimes, combined.summary.measures, final.Ynodes, cum.g, msm.weights, observation.weights) {
  if (isTRUE(attr(data, "called.from.estimate.variance", exact=TRUE))) { 
    return(list(beta=NA, IC=matrix(NA, 1, 1)))
  }
  n <- nrow(data)
  num.regimes <- dim(regimes)[3]
  num.final.Ynodes <- length(final.Ynodes)
  Y.vec <- X.mat <- weight.vec <- NULL
  save.xy <- list()
  for (j in 1:num.final.Ynodes) {
    final.Ynode <- final.Ynodes[j]
    for (i in 1:num.regimes) {
      abar <- drop3(regimes[, , i, drop=F])
      abar <- abar[, nodes$A < final.Ynode, drop=FALSE]
      uncensored <- IsUncensored(data, nodes$C, final.Ynode)
      intervention.match <- InterventionMatch(data, abar, nodes$A, final.Ynode)  
      index <- uncensored & intervention.match
      col.index <- which.max(nodes$AC[nodes$AC < final.Ynode]) 
      
      Y <- data[index, final.Ynode]
      g <- cum.g[index, col.index, i] 
      X <- combined.summary.measures[index, , i, j]
      if (is.vector(X)) { #if only one summary.measure or sum(index)==1, X is dropped to vector
        X <- matrix(X, nrow=sum(index), ncol=ncol(combined.summary.measures))
      }
      weight <- msm.weights[index, i, j] * observation.weights[index] / g
      weight[msm.weights[index, i, j] == 0 | observation.weights[index] == 0] <- 0 #avoid problems where weight and g are both 0
      
      save.xy[[length(save.xy) + 1]] <- list(X=X, Y=Y, weight=weight, index=index)
      Y.vec <- c(Y.vec, Y)
      X.mat <- rbind(X.mat, X)
      weight.vec <- c(weight.vec, weight) 
    }
  }
  colnames(X.mat) <- colnames(combined.summary.measures)
  
  if (nrow(X.mat) == 0) {
    #this happens if there are no rows uncensored and intervention.match
    warning("no rows uncensored and matching regimes/abar - IPTW returns NA")
    num.beta <- ncol(combined.summary.measures)
    return(list(beta=rep(NA, num.beta), IC=matrix(nrow=n, ncol=num.beta)))
  }
  m.glm <- glm(formula(working.msm), family="quasibinomial", data=data.frame(Y=Y.vec, X.mat, weight.vec), weights=scale(weight.vec, center=FALSE)) #note: scale weights because there were rare problems where large weights caused convergence problems
  beta <- coef(m.glm)
  IC <- matrix(0, nrow=n, ncol=length(beta))  #n x num.betas
  m.beta <- array(dim=c(n, num.regimes, num.final.Ynodes)) 
  cnt <- 0
  for (j in 1:num.final.Ynodes) {
    final.Ynode <- final.Ynodes[j]
    for (i in 1:num.regimes) {
      newdata <- data.frame(combined.summary.measures[, , i, j])
      colnames(newdata) <- colnames(combined.summary.measures) #needed if only one summary measure
      SuppressGivenWarnings(m.beta[, i, j] <- predict(m.glm, newdata=newdata, type="response"), "prediction from a rank-deficient fit may be misleading")
      
      cnt <- cnt + 1
      XY.list <- save.xy[[cnt]]
      IC[XY.list$index, ] <- IC[XY.list$index, ] + XY.list$weight * XY.list$X * (XY.list$Y - m.beta[XY.list$index, i, j]) #recycles weight, Y, m.beta
    }
  }

  C <- NormalizeIC(IC, combined.summary.measures, m.beta, msm.weights, g.ratio=array(1, dim=c(n, num.regimes, num.final.Ynodes)), observation.weights=observation.weights) 
  normalized.IC <- t(safe.solve(C, t(IC)))  
  return(list(beta=beta, IC=normalized.IC))
}

# ltmleMSM for a single final.Ynode
FixedTimeTMLE <- function(inputs, msm.weights, combined.summary.measures, baseline.column.names) {
  inputs$summary.measures <- NULL #just to make sure it isn't used - should only use combined.summary.measures 
  #combined.summary.measures: n x num.measures x num.regimes   (num.measures=num.summary.measures + num.baseline.covariates)
  
  data <- inputs$data
  nodes <- inputs$nodes
  
  num.regimes <- dim(inputs$regimes)[3]
  n <- nrow(data)
  num.betas <- ncol(combined.summary.measures)
  tmle <- rep(NA, num.regimes)
  IC <- matrix(0, nrow=n, ncol=num.betas)
  cum.g <- cum.g.unbounded <- prob.A.is.1 <- array(0, dim=c(n, length(nodes$AC), num.regimes))
  cum.g.meanL <- cum.g.meanL.unbounded <- array(0, dim=c(n, length(nodes$AC), num.regimes, length(nodes$LY)-1))
  fit.g <- vector("list", num.regimes)
  for (i in 1:num.regimes) {
    if (all(msm.weights[, i] == 0)) {
      g.list <- list(fit=list("no g fit because regime weight is 0"))
    } else {
      # estimate each g factor, and cumulative probabilities
      g.list <- EstimateG(inputs, regime.index=i)
      cum.g[, , i] <- g.list$cum.g
      cum.g.unbounded[, , i] <- g.list$cum.g.unbounded
      cum.g.meanL[, , i, ] <- g.list$cum.g.meanL
      cum.g.meanL.unbounded[, , i, ] <- g.list$cum.g.meanL.unbounded
      prob.A.is.1[, , i] <- g.list$prob.A.is.1
    } 
    fit.g[[i]] <- g.list$fit
  }
  if (inputs$iptw.only) return(list(cum.g=cum.g, cum.g.unbounded=cum.g.unbounded, IC=NA, Qstar=NA, est.var=NA, g.ratio=NA))
  
  est.var <- matrix(0, num.betas, num.betas)  
  logitQ <- matrix(nrow=n, ncol=num.regimes)
  regimes.with.positive.weight <- which(apply(msm.weights > 0, 2, any))
  if (length(regimes.with.positive.weight) == 0) stop("All regimes have weight 0 (one possible reason is that msm.weights=NULL and no data rows match any of the regimes and are uncensored)")
  fit.Qstar <- vector("list", length(nodes$LY))
  names(fit.Qstar) <- names(data)[nodes$LY]
  fit.Q <- vector("list", length(regimes.with.positive.weight))
  for (i in regimes.with.positive.weight) fit.Q[[i]] <- fit.Qstar
  Qstar.kplus1 <- matrix(data[, max(nodes$Y)], nrow=n, ncol=num.regimes)
  
  for (LYnode.index in length(nodes$LY):1) {
    cur.node <- nodes$LY[LYnode.index]
    deterministic.list.origdata <- IsDeterministic(data, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=FALSE, inputs$survivalOutcome)
    uncensored <- IsUncensored(data, nodes$C, cur.node)
    intervention.match <- subs <- matrix(nrow=n, ncol=num.regimes)
    for (i in regimes.with.positive.weight) {
      abar <- GetABar(inputs$regimes, i)
      intervention.match[, i] <- InterventionMatch(data, abar=abar, nodes$A, cur.node)  
      newdata <- SetA(data, abar=abar, nodes, cur.node)
      deterministic.list.newdata <- IsDeterministic(newdata, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=FALSE, inputs$survivalOutcome)
      if (inputs$stratify) {
        subs[, i] <- uncensored & intervention.match[, i] & !deterministic.list.origdata$is.deterministic
      } else {
        subs[, i] <- uncensored & !deterministic.list.origdata$is.deterministic
      }
      if (any(subs[, i])) {
        Q.est <- Estimate(inputs$Qform[LYnode.index], data=data.frame(data, Q.kplus1=Qstar.kplus1[, i]), family="quasibinomial", newdata=newdata, subs=subs[, i], SL.library=inputs$SL.library.Q, type="link", nodes=nodes, observation.weights=inputs$observation.weights)
        logitQ[, i] <- Q.est$predicted.values
      } else {
        if (! all(deterministic.list.newdata$is.deterministic)) {
          msg <- paste0("ltmle failed trying to estimate ", inputs$Qform[LYnode.index], " because there are no observations that are\nuncensored", ifelse(inputs$stratify, ", follow abar,", ""), " and are not set deterministically due to death or deterministic.Q.function\n")
          stop(msg)
        }
        Q.est <- list(fit="no estimation of Q at this node because all rows are set deterministically")
      }
      fit.Q[[i]][[LYnode.index]] <- Q.est$fit 
    }
    if (all(deterministic.list.newdata$is.deterministic)) {
      #no updating needed if all rows are set deterministically
      Qstar <- matrix(deterministic.list.newdata$Q, nrow=n, ncol=num.regimes)
      if (max(abs(Qstar.kplus1 - Qstar)) > 1e-8) {
        #if Qstar.kplus1 != Qstar when all deterministic score equation will not be solved
        stop("inconsistency in deterministic data - all rows are set deterministically but the deterministically set values are not equal to Qstar.kplus1") 
      }
      Qstar.est <- list(fit="no updating at this node because all rows are set deterministically")
    } else {
      ACnode.index  <- which.max(nodes$AC[nodes$AC < cur.node])
      update.list <- UpdateQ(Qstar.kplus1, logitQ, combined.summary.measures, subs, cum.g[, ACnode.index, ], inputs$working.msm, uncensored, intervention.match, msm.weights, inputs$gcomp, inputs$observation.weights)
      Qstar <- update.list$Qstar
      Qstar[deterministic.list.newdata$is.deterministic, ] <- deterministic.list.newdata$Q
      curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
      curIC.relative.error <- abs(colSums(curIC)) / apply(abs(combined.summary.measures), 2, mean)
      if (any(curIC.relative.error > 0.001) && !inputs$gcomp) {
        SetSeedIfRegressionTesting()
        fix.score.list <- FixScoreEquation(Qstar.kplus1, update.list$h.g.ratio, uncensored, intervention.match, deterministic.list.newdata, update.list$off, update.list$X, regimes.with.positive.weight)
        Qstar <- fix.score.list$Qstar
        curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
        update.list$fit <- fix.score.list$fit      
      }
      est.var <- est.var + EstimateVariance(inputs, combined.summary.measures, regimes.with.positive.weight, uncensored, deterministic.list.newdata, Qstar, Qstar.kplus1, cur.node, msm.weights, LYnode.index, ACnode.index, cum.g, prob.A.is.1, baseline.column.names, cum.g.meanL, cum.g.unbounded, cum.g.meanL.unbounded, inputs$observation.weights, is.last.LYnode=(LYnode.index==length(nodes$LY)))
    }
    IC <- IC + curIC 
    Qstar.kplus1 <- Qstar
    fit.Qstar[[LYnode.index]] <- update.list$fit
  }
  g.ratio <- CalcGUnboundedToBoundedRatio(inputs, cum.g, cum.g.meanL, cum.g.unbounded, cum.g.meanL.unbounded)
  return(list(IC=IC, Qstar=Qstar, cum.g=cum.g, cum.g.unbounded=cum.g.unbounded, g.ratio=g.ratio, est.var=est.var, fit=list(g=fit.g, Q=fit.Q, Qstar=fit.Qstar))) 
}

EstimateVariance <- function(inputs, combined.summary.measures, regimes.with.positive.weight, uncensored, deterministic.list.newdata, Qstar, Qstar.kplus1, cur.node, msm.weights, LYnode.index, ACnode.index, cum.g, prob.A.is.1, baseline.column.names, cum.g.meanL, cum.g.unbounded, cum.g.meanL.unbounded, observation.weights, is.last.LYnode) {
  if (inputs$IC.variance.only) return(NA)
  TmleOfVariance <- function(Z, Z.meanL) {
    if (all(is.na(Z))) stop("all Z are NA in EstimateVariance")
    #length(Z.meanL) == 0 --> point treatment case, return mean(Z); all Z can be zero when summary.measure is 0
    if (length(Z.meanL) == 0 || all(Z==0 | is.na(Z))) {
      Qstar <- Scale(Z, 0, 1)
      return(list(EZd1 = mean(Z, na.rm=T), Qstar = Qstar))
    }
    sparsity.data <- inputs$data[, 1:cur.node]
    sparsity.data[, cur.node] <- Scale(Z, 0, 1)
    temp.nodes <- lapply(nodes, function (x) x[x <= cur.node])
    if (cur.node %in% temp.nodes$L) {
      #if the current node is an L node, make current node a Y node (last node has to be a Y node)
      temp.nodes$L <- setdiff(temp.nodes$L, cur.node)
      temp.nodes$Y <- c(temp.nodes$Y, cur.node)
    }
    stratify <- FALSE
    Qform <- paste(GetDefaultForm(sparsity.data[, 1:cur.node], nodes=temp.nodes, is.Qform=TRUE, stratify=stratify, survivalOutcome=FALSE, showMessage=FALSE), paste0("+ sparityAdj_Z.meanL_", 1:length(temp.nodes$LY)))
    Qform[length(Qform)] <- "IDENTITY"
    
    Z.meanL <- apply(AsMatrix(Z.meanL), 2, LogitScale)
    sparsity.data <- cbind(Z.meanL, sparsity.data)
    names(sparsity.data)[sseq(1, ncol(Z.meanL))] <- paste0("sparityAdj_Z.meanL_", sseq(1, ncol(Z.meanL)))
    temp.nodes <- lapply(temp.nodes, function (x) x + ncol(Z.meanL))
    
    names(Qform) <- names(sparsity.data)[temp.nodes$LY]
    attr(sparsity.data, "called.from.estimate.variance") <- TRUE
    var.tmle <- ltmle(sparsity.data, Anodes=temp.nodes$A, Cnodes=temp.nodes$C, Lnodes=temp.nodes$L, Ynodes=temp.nodes$Y, survivalOutcome=FALSE, Qform=Qform, gform=drop3(prob.A.is.1[, 1:ACnode.index, d1, drop=FALSE]), abar=drop3(inputs$regimes[, nodes$A <= cur.node, d1, drop=FALSE]), gbounds=inputs$gbounds, stratify=stratify, estimate.time=FALSE, deterministic.Q.function=det.q.function, IC.variance.only=TRUE, observation.weights=observation.weights) 
  
    EZd1 <- var.tmle$estimates["tmle"] * diff(range(Z, na.rm=T)) + min(Z, na.rm=T)
    return(list(EZd1 = EZd1, Qstar=var.tmle$Qstar))
  }
  
  EqualRegimesIndex <- function(dd1, dd2) {
    #index of each observation where regime d1 matches regime d2 up to cur.node
    return(apply(drop3(inputs$regimes[, inputs$nodes$A <= cur.node, dd1, drop=F]) == drop3(inputs$regimes[, inputs$nodes$A <= cur.node, dd2, drop=F]), 1, all)) 
  }
  
  IsStaticTreatment <- function() {
    #static = for all observations, regime d1 matches regime d2 only when d1=d2
    for (dd1 in regimes.with.positive.weight) {  
      for (dd2 in regimes.with.positive.weight[regimes.with.positive.weight > dd1]) {
        if (any(EqualRegimesIndex(dd1, dd2))) return(FALSE)
      }
    }
    return(TRUE)
  }
  
  nodes <- inputs$nodes  
  num.regimes <- dim(inputs$regimes)[3]
  num.betas <- ncol(combined.summary.measures)
  n <- nrow(inputs$data)
  
  #used in ltmle call below
  if (inputs$survivalOutcome) {
    det.q.function <- function(data, current.node, nodes, called.from.estimate.g) {
      if (!any(nodes$Y < current.node)) return(NULL)
      prev.Y <- data[, nodes$Y[nodes$Y < current.node], drop=F]
      prev.Y[is.na(prev.Y)] <- 0
      is.deterministic <- apply(prev.Y == 1, 1, any)
      Q.value <- data[is.deterministic, max(nodes$Y)] #this is 0 before scaling but may be nonzero after scaling
      return(list(is.deterministic=is.deterministic, Q.value=Q.value))   
    }
  } else {
    det.q.function <- NULL
  }
  static.treatment <- IsStaticTreatment()
  variance.estimate <- matrix(0, num.betas, num.betas)
  alive <- !deterministic.list.newdata$is.deterministic
  Sigma <- array(dim=c(n, num.regimes, num.regimes))
  #fixme - Sigma is not exactly symmetric due to SetA on d1, but could probably save time with approximate symmetry
  for (d1 in regimes.with.positive.weight) {
    if (static.treatment) {
      d2.regimes <- d1 #only need diagonal elements
    } else {
      d2.regimes <- regimes.with.positive.weight #need all elements
    }
    for (d2 in d2.regimes) { 
      if (is.last.LYnode) {
        Sigma[, d1, d2] <- Qstar[, d1] * (1 - Qstar[, d1]) 
      } else {
        Q.data <- inputs$data[alive, 1:cur.node, drop=F]
        resid.sq <- (Qstar.kplus1[alive, d1] - Qstar[alive, d1]) * (Qstar.kplus1[alive, d2] - Qstar[alive, d2]) 
        resid.sq.range <- range(resid.sq, na.rm=T)
        if (diff(resid.sq.range) > 0) {
          Q.data[, cur.node] <- (resid.sq - resid.sq.range[1]) / diff(resid.sq.range)
          names(Q.data)[cur.node]  <- "Q.kplus1" #ugly - to match with Qform
          m <- glm(formula = inputs$Qform[LYnode.index], family = "quasibinomial", data = Q.data, control=glm.control(trace=FALSE, maxit=1000)) 
          Q.newdata <- SetA(data = Q.data, abar = GetABar(regimes = inputs$regimes, d1)[alive, , drop=F], nodes = nodes, cur.node = cur.node)
          SuppressGivenWarnings(Q.resid.sq.pred <- predict(m, newdata = Q.newdata, type = "response"), "prediction from a rank-deficient fit may be misleading")
          Sigma[alive, d1, d2] <- Q.resid.sq.pred * diff(resid.sq.range) + resid.sq.range[1]
        } else {
          resid.sq.value <- min(resid.sq, na.rm = T) #all values are the same, just get one non-NA
          Sigma[alive, d1, d2] <- resid.sq.value
        }
        Sigma[!alive, d1, d2] <- 0
      }
    }
  }
  

  if (static.treatment) {
    for (d1 in regimes.with.positive.weight) {   
      #Z.without.h1h1 <- Sigma[, d1, d1] / cum.g[, ACnode.index, d1] #without h1*h1'
      #Z.without.h1h1.meanL <- 1 / cum.g.meanL[, ACnode.index, d1, ]
      Z.without.sum.meas <- Sigma[, d1, d1] / cum.g[, ACnode.index, d1] * cum.g.unbounded[, ACnode.index, d1] / cum.g[, ACnode.index, d1] * msm.weights[, d1]^2 * observation.weights^2
      Z.without.sum.meas.meanL <- 1 / cum.g.meanL[, ACnode.index, d1, ] * cum.g.meanL.unbounded[, ACnode.index, d1, ] / cum.g.meanL[, ACnode.index, d1, ] * msm.weights[, d1]^2 * observation.weights^2
      var.tmle <- TmleOfVariance(Z.without.sum.meas, Z.without.sum.meas.meanL)
      no.V <- all(combined.summary.measures[1, , d1] == combined.summary.measures[, , d1]) 
      if (no.V) {
        #h1 <- combined.summary.measures[1, , d1] * msm.weights[1, d1]  #num.betas x 1
        #variance.estimate <- variance.estimate + (h1 %*% t(h1)) * var.tmle$EZd1
        variance.estimate <- variance.estimate + (combined.summary.measures[1, , d1] %*% t(combined.summary.measures[1, , d1])) * var.tmle$EZd1
      } else {
        #has V (so combined.summary.measures varies)
        baseline.msm <- "Qstar ~ 1"
        if (length(baseline.column.names) > 0) {
          baseline.msm <- paste(baseline.msm, "+", paste(baseline.column.names, collapse=" + "), "+", paste0("I(", baseline.column.names, "^2)", collapse=" + "))
        }            
        m <- glm(baseline.msm, family = "quasibinomial", data=data.frame(Qstar=var.tmle$Qstar, inputs$data[, baseline.column.names, drop=FALSE]), control=glm.control(trace=FALSE, maxit=1000)) 
        SuppressGivenWarnings(pred.Qstar <- predict(m, type = "response") * diff(range(Z.without.sum.meas, na.rm=T)) + min(Z.without.sum.meas, na.rm=T), "prediction from a rank-deficient fit may be misleading")  #n x 1
        variance.estimate.sum <- matrix(0, num.betas, num.betas)
        for (i in 1:n) {
          #h1 <- combined.summary.measures[i, , d1] * msm.weights[i, d1]
          #variance.estimate.sum <- variance.estimate.sum + (h1 %*% t(h1)) * pred.Qstar[i]
          variance.estimate.sum <- variance.estimate.sum + (combined.summary.measures[i, , d1] %*% t(combined.summary.measures[i, , d1])) * pred.Qstar[i]
        }
        variance.estimate <- variance.estimate + variance.estimate.sum / n
      }
    }
  } else {
    for (beta.index2 in 1:num.betas) {
      for (d1 in regimes.with.positive.weight) {          
        Z.base <- rep(0, n)  #Z without h1(d1, V, beta.index1)
        Z.base.meanL <- matrix(0, n, dim(cum.g.meanL)[4])
        for (d2 in regimes.with.positive.weight) {
          equal.regimes.index <- EqualRegimesIndex(d1, d2) #index of each observation where regime d1 matches regime d2 
          h1 <- combined.summary.measures[, beta.index2, d2] * msm.weights[, d2]
          Z.base[equal.regimes.index] <- Z.base[equal.regimes.index] + h1[equal.regimes.index] * Sigma[equal.regimes.index, d1, d2] / cum.g[equal.regimes.index, ACnode.index, d1] * observation.weights[equal.regimes.index] #this is equivalent to using cum.g.unbounded in the denominator and multiplying by phi=cum.g.unbounded/cum.g.bounded  
          Z.base.meanL[equal.regimes.index, ] <- Z.base.meanL[equal.regimes.index, ] + h1[equal.regimes.index] * 1 / cum.g.meanL[equal.regimes.index, ACnode.index, d1, ] * observation.weights[equal.regimes.index] #recycles
        }        
        for (beta.index1 in 1:num.betas) {  
          if (beta.index1 >= beta.index2) {
            Z <- combined.summary.measures[, beta.index1, d1] * msm.weights[, d1] * cum.g.unbounded[, ACnode.index, d1] / cum.g[, ACnode.index, d1] * observation.weights * Z.base 
            Z.meanL <- combined.summary.measures[, beta.index1, d1] * msm.weights[, d1] * cum.g.meanL.unbounded[, ACnode.index, d1, ] / cum.g.meanL[, ACnode.index, d1, ] * observation.weights * Z.base.meanL
            var.tmle <- TmleOfVariance(Z, Z.meanL)
            variance.estimate[beta.index1, beta.index2] <- variance.estimate[beta.index1, beta.index2] + var.tmle$EZd1
          } else {
            variance.estimate[beta.index1, beta.index2] <- variance.estimate[beta.index2, beta.index1] #use symmetry
          }
        }
      } 
    }
  }
  if (max(abs(variance.estimate - t(variance.estimate))) > 1e-5) stop("not symmetric")  
  if (any(eigen(variance.estimate, only.values = TRUE)$values < -1e-8)) {
    #this happens very rarely
    orig.variance.estimate <- variance.estimate
    try.result <- try({
      near.pd <- Matrix::nearPD(variance.estimate) #may cause an error if variance.estimate is negative definite
      variance.estimate <- as.matrix(near.pd$mat)
    })
    if (inherits(try.result, "try-error") || !near.pd$converged || any((abs(orig.variance.estimate - variance.estimate) / orig.variance.estimate) > 0.1)) {
      warning("Covariance matrix from EstimateVariance not positive definite, unable to compute standard errors. You may want to try IC.variance.only=TRUE.")
      variance.estimate <- matrix(nrow=num.betas, ncol=num.betas)
    }
  }
  return(variance.estimate)
}

CalcGUnboundedToBoundedRatio <- function(inputs, cum.g, cum.g.meanL, cum.g.unbounded, cum.g.meanL.unbounded) {
  n <- dim(cum.g)[1]
  num.AC.nodes <- dim(cum.g)[2]
  num.regimes <- dim(cum.g)[3]
  if (! any(is.na(cum.g))) return(AsMatrix(cum.g.unbounded[, num.AC.nodes, ] / cum.g[, num.AC.nodes, ]))
  #cum.g is NA after censoring - for censored observations use cum.g.meanL  
  #If censored at node j, set all nodes > j to meanl. 
  #[,,k] is prob.A.is.1 with all L and Y nodes after and including LYnodes[k] set to mean of L (na.rm=T)
  g.ratio <- matrix(NA, n, num.regimes)
  for (i in 1:num.regimes) {
    g.ratio.temp <- cbind(cum.g.meanL.unbounded[, num.AC.nodes, i, ] / cum.g.meanL[, num.AC.nodes, i, ], cum.g.unbounded[, num.AC.nodes, i] / cum.g[, num.AC.nodes, i])
    index <- max.col(!is.na(g.ratio.temp), "last")
    g.ratio[, i] <- g.ratio.temp[sub2ind(1:n, col = index, num.rows = n)]
  }
  return(g.ratio)
}

# remove any information in ltmleInputs after final.Ynode
SubsetInputs <- function(inputs, final.Ynode) {
  if (is.numeric(inputs$gform)) {
    stopifnot(length(dim(inputs$gform)) == 3)
    inputs$gform <- inputs$gform[, inputs$nodes$AC < final.Ynode, , drop=FALSE]
  } else {
    inputs$gform <- inputs$gform[inputs$nodes$AC < final.Ynode]
  }
  inputs$Qform <- inputs$Qform[inputs$nodes$LY <= final.Ynode]
  inputs$data <- inputs$data[, 1:final.Ynode, drop=FALSE]
  inputs$untransformed.data <- inputs$untransformed.data[, 1:final.Ynode, drop=FALSE]
  inputs$regimes <- inputs$regimes[, inputs$nodes$A <= final.Ynode, , drop=FALSE]
  inputs$summary.measures <- drop3(inputs$summary.measures[, , inputs$final.Ynodes == final.Ynode, drop=FALSE])
  inputs$final.Ynodes <- inputs$final.Ynodes[inputs$final.Ynodes <= final.Ynode]
  inputs$nodes <- lapply(inputs$nodes, function (x) x[x <= final.Ynode])
  
  return(inputs)
}

# Fit the MSM
FitPooledMSM <- function(working.msm, Qstar, combined.summary.measures, msm.weights) {
  #Qstar: n x num.regimes x num.final.Ynodes
  #combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes   (num.measures=num.summary.measures + num.baseline.covariates)
  #msm.weights: n x num.regimes x num.final.Ynodes

  n <- dim(Qstar)[1]
  num.regimes <- dim(Qstar)[2]
  num.final.Ynodes <- dim(Qstar)[3]
  num.summary.measures <- dim(combined.summary.measures)[2]
  
  X <- apply(combined.summary.measures, 2, rbind) 
  Y <- as.vector(Qstar)
  weight.vec <- as.vector(msm.weights)
  
  m <- glm(as.formula(working.msm), data=data.frame(Y, X), family="quasibinomial", weights=scale(weight.vec, center = FALSE), na.action=na.exclude, control=glm.control(maxit=1000)) 
  SuppressGivenWarnings(m.beta <- predict(m, type="response"), "prediction from a rank-deficient fit may be misleading")
  dim(m.beta) <- dim(Qstar)
  return(list(m=m, m.beta=m.beta))
}

#final step in calculating TMLE influence curve
FinalizeIC <- function(IC, combined.summary.measures, Qstar, m.beta, msm.weights, g.ratio, observation.weights) {
  #mBeta, Qstar: n x num.regimes x num.final.Ynodes
  #combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes   (num.measures=num.summary.measures + num.baseline.covariates)
  
  #summary.measures: num.regimes x num.summary.measures x num.final.Ynodes
  #msm.weights: n x num.regimes x num.final.Ynodes
  
  num.betas <- ncol(IC)
  n <- nrow(Qstar)
  num.regimes <- ncol(Qstar)
  num.final.Ynodes <- dim(Qstar)[3]
  
  stopifnot(num.betas == ncol(combined.summary.measures))

  finalIC <- matrix(0, nrow=n, ncol=num.betas)
  for (j in 1:num.final.Ynodes) {
    for (i in 1:num.regimes) {
      if (any(msm.weights[, i, j] > 0)) {
        m1 <- matrix(Qstar[, i, j] - m.beta[, i, j], ncol=1)   #n x 1
        for (k in 1:num.betas) {
          m2 <- combined.summary.measures[, k, i, j] # n x 1
          finalIC[, k] <- finalIC[, k] + msm.weights[, i, j] * observation.weights * (m1 * m2) 
        }
      }
    }  
  }
  if (any(abs(colSums(finalIC)) > 0.001 )) {
    msg <- capture.output(cat("final IC problem", colSums(finalIC)))
    warning(paste(msg, collapse="\n"))
  }
  IC <- IC + finalIC
  return(IC) 
}

# Normalize the influence curve matrix
NormalizeIC <- function(IC, combined.summary.measures, m.beta, msm.weights, g.ratio, observation.weights) {    
  #combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes   (num.measures=num.summary.measures + num.baseline.covariates)
  #g.ratio = g.unbounded / g.bounded : n x num.regimes x num.final.Ynodes
  n <- nrow(IC)
  num.betas <- ncol(IC)
  num.regimes <- dim(combined.summary.measures)[3]
  num.final.Ynodes <- dim(combined.summary.measures)[4]
  
  C <- array(0, dim=c(num.betas, num.betas, n))
  for (j in 1:num.final.Ynodes) {
    for (i in 1:num.regimes) {
      positive.msm.weights <- which(msm.weights[, i, j] > 0)
      tempC <- array(0, dim=c(num.betas, num.betas, n))
      for (k in positive.msm.weights) {
        if (max(m.beta[, i, j]) > (1 - 1e-6) || min(m.beta[, i, j]) < 1e-6) {
          warning("All predicted probabilities are all very close to 0 or 1. Unable to compute standard errors.")
          return(matrix(NA, nrow=num.betas, ncol=num.betas))
        }
        m.beta.temp <- m.beta[k, i, j]  
        h <- matrix(combined.summary.measures[k, , i, j], ncol=1) * msm.weights[k, i, j] * g.ratio[k, i, j]
        
        tempC[, , k] <- h %*% t(h) * m.beta.temp * (1 - m.beta.temp) * observation.weights[k] / msm.weights[k, i, j]
      }
      if (any(is.na(tempC))) stop("NA in tempC")
      C <- C + tempC
    }
  }
  C <- apply(C, c(1, 2), mean)
  if (rcond(C) < 1e-12) {
    C <- matrix(NA, nrow=num.betas, ncol=num.betas)
    warning("rcond(C) near 0, standard errors not available")
  } else {
    normalized.IC <- t(safe.solve(C, t(IC))) #IC %*% solve(C) 
    if (!any(abs(colSums(IC)) > 0.001) && !any(is.na(normalized.IC)) && any(abs(colSums(normalized.IC)) > 0.001)) {
      msg <- capture.output({
        cat("normalized IC problem", colSums(normalized.IC), "\n")
        cat("inv(C) = \n")
        print(safe.solve(C))
      })
      warning(paste(msg, collapse="\n"))
    }
  }
  return(C)
}

# Get a single regime from the regimes array
GetABar <- function(regimes, i) {
  abar <- AsMatrix(regimes[, , i]) #if there's only 1 Anode, make sure abar comes back as a matrix
  return(abar)
}

# Targeting step - update the initial fit of Q using clever covariates
UpdateQ <- function(Qstar.kplus1, logitQ, combined.summary.measures, subs, cum.g, working.msm, uncensored, intervention.match, msm.weights, gcomp, observation.weights) { 
  #logitQ, Qstar.kplus1: n x num.regimes
  #cum.g: n x num.regimes (already indexed for this node)
  #subs: n x num.regimes
  #uncensored: n x 1
  #intervention.match: n x num.regimes
  #summary.measures: num.regimes x num.summary.measures
  #baseline.covariates: names/indicies: num.baseline.covariates x 1
  #msm.weights: n x num.regimes
  #stacked.summary.measures: (n*num.regimes) x num.measures
  #combined.summary.measures: n x num.measures x num.regimes   (num.measures=num.summary.measures + num.baseline.covariates)
  #h.g.ratio: n x num.regimes x num.measures
  #observation.weights: n x 1
  n <- nrow(logitQ)
  num.regimes <- ncol(logitQ)
  off <- as.vector(logitQ)
  Y <- as.vector(Qstar.kplus1)
  
  stacked.summary.measures <- apply(combined.summary.measures, 2, rbind)
  
  weight.vec <- uncensored * observation.weights * as.vector(intervention.match) / as.vector(cum.g) * as.vector(msm.weights) #recycles uncensored and observation.weights
  f <- as.formula(paste(working.msm, "+ offset(off)"))
  data.temp <- data.frame(Y, stacked.summary.measures, off)
  newdata <- data.temp  
  if (gcomp) {
    Qstar <- plogis(logitQ)
    m <- "no Qstar fit because gcomp=TRUE (so no updating step)"
  } else {
    ctrl <- glm.control(trace=FALSE, maxit=1000)
    
    if (any(subs & weight.vec>0)) {
      SuppressGivenWarnings(m <- glm(f, data=data.temp, subset=as.vector(subs) & (weight.vec > 0), family="quasibinomial", weights=scale(weight.vec, center=FALSE), control=ctrl), GetWarningsToSuppress(TRUE)) #this should include the indicators; only include weight.vec>0 because others have NAs
      SuppressGivenWarnings(Qstar <- matrix(predict(m, newdata=newdata, type="response"), nrow=nrow(logitQ)), GetWarningsToSuppress(TRUE))  #this should NOT include the indicators  #note: could also use plogis(off + X %*% coef(m)) [but this has problems with NAs in coef(m)?]
    } else {
      Qstar <- plogis(logitQ)
      m <- "no Qstar fit because no subjects alive, uncensored, following intervention"
    }
    
  }
  indicator <- matrix(uncensored * observation.weights, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) * matrix(intervention.match, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) #I(A=rule and uncensored) * observation.weights
  h.g.ratio <- stacked.summary.measures / matrix(cum.g, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) * indicator # I() * h * observation.weights / g
  dim(h.g.ratio) <- c(n, num.regimes, ncol(h.g.ratio))
  for (i in 1:num.regimes) {
    h.g.ratio[, i, ] <- h.g.ratio[, i, ] * msm.weights[, i] #recycles msm.weights
    weight.zero.index <- msm.weights[, i] == 0
    h.g.ratio[weight.zero.index, i, ] <- 0  #cum.g is 0 so X is NA so h.g.ratio is NA when weight is 0
  }
  if (ncol(stacked.summary.measures) != dim(h.g.ratio)[3]) stop("only works if working.msm has only main terms")
  return(list(Qstar=Qstar, h.g.ratio=h.g.ratio, X=stacked.summary.measures, off=off, fit=m)) 
}

# Sometimes GLM doesn't converge and the updating step of TMLE doesn't solve the score equation (sum of TMLE influence curve not equal to zero). This function attempts to solve the score equation directly using various optimizers. [Note: this was needed more often before we made changes to the updating algorithm.]
FixScoreEquation <- function(Qstar.kplus1, h.g.ratio, uncensored, intervention.match, deterministic.list, off, X, regimes.with.positive.weight) {
  CalcScore <- function(e) {
    Qstar <- QstarFromE(e)
    ICtemp <- CalcIC(Qstar.kplus1, Qstar, h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
    return(sum(colSums(ICtemp) ^ 2)) #each column has to add to zero
  }
  
  QstarFromE <- function(e) {
    Qstar <- plogis(off + X %*% e) #X: n x (num.summary.measures + num.baseline.covariates) (which should be num.beta);  e: num.beta x 1 
    dim(Qstar) <- dim(Qstar.kplus1)
    Qstar[deterministic.list$is.deterministic, ] <- deterministic.list$Q
    return(Qstar)
  }
  
  FindMin <- function(minimizer) {
    num.tries <- 30 #increase this to solve more problems at the cost of longer run time
    init.e <- numeric(num.betas) #first try an initial estimate of epsilon=0
    for (i in 1:num.tries) {
      if (minimizer == "nlminb") {
        m <- nlminb(start=init.e, objective=CalcScore, control=list(abs.tol=max.objective, eval.max=500, iter.max=500, x.tol=1e-14, rel.tol=1e-14))
        e <- m$par
        obj.val <- m$objective
      } else if (minimizer == "optim") {
        if (num.betas == 1) {
          m <- optimize(f = CalcScore, interval=c(-1, 1) * 10^i)
          e <- m$minimum
          obj.val <- m$objective
        } else {
          m <- optim(par=init.e, fn=CalcScore, control=list(abstol=max.objective, reltol=1e-14, maxit=2000))
          e <- m$par
          obj.val <- m$value
        }
      } else if (minimizer == "nlm") {
        m <- nlm(f=CalcScore, p=init.e)
        e <- m$estimate
        obj.val <- m$minimum
      } else {
        stop("bad minimizer")
      }
      if (obj.val < max.objective) {
        m$ltmle.msg <- paste("updating step using glm failed to solve score equation; solved using", minimizer)
        return(list(e=e, solved=TRUE, m=m))
      }
      init.e <- rnorm(num.betas) #if the first try didn't work, try a random initial estimate of epsilon
    }
    return(list(e=numeric(num.betas), solved=FALSE, m="score equation not solved!")) #return Q (not updated)
  }
  max.objective <- 0.0001 ^ 2
  num.betas <- ncol(X)
  
  for (offset.lbound in c(1e-8, 0.0001, 0.001, 0.01)) {
    off <- Bound(off, qlogis(c(offset.lbound, 1-offset.lbound)))
    l <- FindMin("nlminb")
    if (! l$solved) l <- FindMin("optim")
    if (! l$solved) l <- FindMin("nlm")
    if (l$solved) break
  }
  if (! l$solved) {
    stop("all minimizers failed")
  }
  Qstar <- QstarFromE(l$e)
  return(list(Qstar=Qstar, fit=l$m))
}

# Estimate how long it will take to run ltmleMSM
EstimateTime <- function(inputs) {
  sample.size <- 50
  if (nrow(inputs$data) < sample.size) {
    message(paste("Timing estimate unavailable when n <", sample.size))
    return(NULL)
  }
  sample.index <- sample(nrow(inputs$data), size=sample.size)
  
  small.inputs <- inputs
  small.inputs$data <- inputs$data[sample.index, ]
  small.inputs$regimes <- inputs$regimes[sample.index, , , drop=F]
  small.inputs$observation.weights <- inputs$observation.weights[sample.index]
  if (is.numeric(inputs$gform)) small.inputs$gform <- inputs$gform[sample.index, , , drop=F]
  if (length(dim(inputs$msm.weights)) == 3) small.inputs$msm.weights <- inputs$msm.weights[sample.index, , , drop=F]
  start.time <- Sys.time()
  #try.result <- MainCalcs(small.inputs)
  try.result <- suppressWarnings(try(MainCalcs(small.inputs), silent=TRUE))
  if (inherits(try.result, "try-error")) {
    message("Timing estimate unavailable")
  } else {
    elapsed.time <- Sys.time() - start.time 
    est.time1 <- round(sqrt(as.double(elapsed.time, units="mins") * nrow(inputs$data) / sample.size), digits=0)
    est.time2 <- round(as.double(elapsed.time, units="mins") * nrow(inputs$data) / sample.size, digits=0)
    if (est.time2 == 0) {
      est.time.str <- "< 1 minute"
    } else if (est.time2 == 1) {
      est.time.str <- "1 minute"
    } else {
      est.time.str <- paste(est.time1, "to", est.time2, "minutes")
    }
    message("Estimate of time to completion: ", est.time.str)
  }
  return(NULL)
}

#' Get standard error, p-value, and confidence interval for one ltmle object 
#' Summarizing results from Longitudinal Targeted Maximum Likelihood Estimation
#' (ltmle)
#' 
#' These functions are methods for class \code{ltmle} or \code{summary.ltmle}
#' objects.
#' 
#' \code{summary.ltmle} returns the parameter value of the estimator, the
#' estimated variance, a 95 percent confidence interval, and a p-value.
#' 
#' \code{summary.ltmleEffectMeasures} returns the additive treatment effect for
#' each of the two objects in the \code{abar} list passed to \code{ltmle}.
#' Relative risk, and odds ratio are also returned, along with the variance,
#' confidence interval, and p-value for each.
#' 
#' \code{summary.ltmleMSM} returns a matrix of MSM parameter estimates.
#' 
#' @aliases summary.ltmle print.ltmle print.summary.ltmle summary.ltmleMSM
#' print.ltmleMSM print.summary.ltmleMSM summary.ltmleEffectMeasures
#' print.ltmleEffectMeasures print.summary.ltmleEffectMeasures
#' @param object an object of class "\code{ltmle}" or "\code{ltmleMSM}" or
#' "\code{ltmleEffectMeasures}", usually a result of a call to
#' \code{\link{ltmle}} or \code{\link{ltmleMSM}}.
#' @param x an object of class "\code{summary.ltmle}" or
#' "\code{summary.ltmleMSM}" or "\code{ltmleEffectMeasures}", usually a result
#' of a call to \code{\link{summary.ltmle}} or \code{\link{summary.ltmleMSM}}.
#' @param estimator character; one of "tmle", "iptw", "gcomp". The estimator
#' for which to get effect measures. "tmle" is valid iff the original
#' ltmle/ltmleMSM call used gcomp=FALSE. "gcomp" is valid iff the original
#' ltmle/ltmleMSM call used gcomp=TRUE
#' @param digits the number of significant digits to use when printing.
#' @param signif.stars logical. If \code{TRUE}, significance stars are printed
#' for each coefficient.
#' @param \dots further arguments passed to or from other methods.
#' @return \code{summary.ltmle} returns an object of class
#' "\code{summary.ltmle}", a list with components \item{treatment}{a list with
#' components summarizing the estimate of \code{object} \itemize{ \item
#' \code{estimate} - the parameter estimate of \eqn{E[Y_d]} \item
#' \code{std.dev} - estimated standard deviation of parameter \item
#' \code{p.value} - two-sided p-value \item \code{CI} - vector of length 2 with
#' 95 percent confidence interval } }
#' 
#' \item{call}{the matched call to \code{ltmle} for \code{object}}
#' \item{estimator}{the \code{estimator} input argument}
#' \item{variance.estimate.ratio}{ratio of the TMLE based variance estimate to
#' the influence curve based variance estimate}
#' 
#' \code{summary.ltmleEffectMeasures} returns an object of class
#' "\code{summary.ltmleEffectMeasures}", a list with same components as
#' \code{summary.ltmle} above, but also includes: \item{effect.measures}{a list
#' with components, each with the same components as \code{treatment} in
#' \code{summary.ltmle} above \itemize{ \item \code{treatment} - corresponds to
#' the first in the list \code{abar} (or \code{rule}) passed to \code{ltmle}
#' \item \code{control} - corresponds to the second in the list \code{abar} (or
#' \code{rule}) passed to \code{ltmle} \item \code{ATE} - average treatment
#' effect \item \code{RR} - relative risk \item \code{OR} - odds ratio } }
#' 
#' \code{summary.ltmleMSM} returns an object of class
#' "\code{summary.ltmleMSM}", a matrix with rows for each MSM parameter and
#' columns for the point estimate, standard error, 2.5percent confidence
#' interval, 97.5percent confidence interval, and p-value.
#' @seealso \code{\link{ltmle}}, \code{\link{summary}}
#' @examples
#' 
#' rexpit <- function(x) rbinom(n = length(x), size = 1, prob = plogis(x))
#' 
#' # Compare the expected outcomes under two counterfactual plans: Treatment plan:
#' # set A1 to 1 if W > 0, set A2 to 1 if W > 1.5, always set A3 to 1 Control plan:
#' # always set A1, A2, and A3 to 0
#' W <- rnorm(1000)
#' A1 <- rexpit(W)
#' A2 <- rexpit(W + 2 * A1)
#' A3 <- rexpit(2 * A1 - A2)
#' Y <- rexpit(W - A1 + 0.5 * A2 + 2 * A3)
#' data <- data.frame(W, A1, A2, A3, Y)
#' treatment <- cbind(W > 0, W > 1.5, 1)
#' control <- matrix(0, nrow = 1000, ncol = 3)
#' result <- ltmle(data, Anodes = c("A1", "A2", "A3"), Ynodes = "Y", abar = list(treatment, 
#'     control))
#' print(summary(result))
#' 
#' ## For examples of summary.ltmle and summary.ltmleMSM, see example(ltmle)
#' 
#' @export 
summary.ltmle <- function(object, estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...) {
  if ("control.object" %in% names(list(...))) stop("The control.object parameter has been deprecated. To obtain additive treatment effect, risk ratio, and relative risk, call ltmle with abar=list(treatment, control). See ?ltmle and ?summary.ltmleEffectMeasures.")
  if (! estimator[1] %in% c("tmle", "iptw", "gcomp")) stop("estimator should be one of: tmle, iptw, gcomp. If you are trying to use control.object, the control.object parameter has been deprecated. To obtain additive treatment effect, risk ratio, and relative risk, call ltmle with abar=list(treatment, control). See ?ltmle and ?summary.ltmleEffectMeasures.")
  if (estimator == "tmle" && object$gcomp) stop("estimator 'tmle' is not available because ltmleMSM was called with gcomp=TRUE")
  if (estimator == "gcomp" && !object$gcomp) stop("estimator 'gcomp' is not available because ltmleMSM was called with gcomp=FALSE")

  IC.variance <- var(object$IC[[estimator]])
  if (estimator=="tmle" && !is.null(object$variance.estimate)) {
    v <- max(IC.variance, object$variance.estimate) 
  } else {
    v <- IC.variance
  }
  variance.estimate.ratio=v/IC.variance
  
  if (object$binaryOutcome) {
    CIBounds <- c(0, 1)
  } else {
    CIBounds <- c(-Inf, Inf)  #could truncate at Yrange, but it's not clear that's right
  }
  treatment <- GetSummary(list(long.name=NULL, est=object$estimates[estimator], gradient=1, log.std.err=FALSE, CIBounds=CIBounds), v, n=length(object$IC[[estimator]]))
  ans <- list(treatment=treatment, call=object$call, estimator=estimator, variance.estimate.ratio=variance.estimate.ratio)
  class(ans) <- "summary.ltmle"
  return(ans)
}

# Get standard errors, p-values, confidence intervals for an ltmleEffectMeasures object: treatment EYd, control EYd, additive effect, relative risk, odds ratio
#' @rdname summary.ltmle
#' @export 
summary.ltmleEffectMeasures <- function(object, estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...) {
  info <- GetSummaryLtmleMSMInfo(object, estimator)
  beta <- info$estimate
  IC <- info$IC
  y0 <- plogis(beta[1])
  y1 <- plogis(beta[1] + beta[2])
  eff.list <- list(
    treatment=list(long.name="Treatment Estimate", est=y1, gradient=c(y1*(1-y1), y1*(1-y1)), log.std.err=FALSE, CIBounds=0:1),  
    control=list(long.name="Control Estimate", est=y0, gradient=c(y0*(1-y0), 0), log.std.err=FALSE, CIBounds=0:1), 
    ATE=list(long.name="Additive Treatment Effect", est=y1 - y0, gradient=c(y1*(1-y1) - y0*(1-y0), y1*(1-y1)), log.std.err=FALSE, CIBounds=c(-1, 1)), 
    RR=list(long.name="Relative Risk", est=y1 / y0, gradient=c(y0-y1, 1-y1), log.std.err=TRUE, CIBounds=c(0, Inf)), 
    OR=list(long.name="Odds Ratio", est=exp(beta[2]), gradient=c(0, 1), log.std.err=TRUE, CIBounds=c(0, Inf)))
  if (!object$binaryOutcome) {
    eff.list$RR <- eff.list$OR <- NULL #not valid if non-binary outcome
  }
  n <- nrow(IC)
  
  measures.IC <- lapply(eff.list, GetSummary, var(IC), n)
  if (is.null(object$variance.estimate)) {
    measures.variance.estimate <- NULL #if IC.variance.only=T
  } else {
    measures.variance.estimate <- lapply(eff.list, GetSummary, object$variance.estimate, n)  
  }
  measures.max <- measures.IC
  for (i in seq_along(measures.variance.estimate)) {
    std.dev.diff <- measures.variance.estimate[[i]]$std.dev - measures.IC[[i]]$std.dev
    if (!is.na(std.dev.diff) && (std.dev.diff > 0)) { #can be NA if all Y_d are near 0 or 1
      measures.max[[i]] <- measures.variance.estimate[[i]]
    }
  }
  
  ans <- list(call=object$call, effect.measures=measures.max, variance.estimate.ratio=info$variance.estimate.ratio, transformOutcome=object$transformOutcome, estimator=estimator)
  class(ans) <- "summary.ltmleEffectMeasures"
  return(ans) 
}

# Do some error checking and get basic info about the ltmleMSM object
GetSummaryLtmleMSMInfo <- function(object, estimator) {
  if (! estimator %in% c("tmle", "iptw", "gcomp")) stop("estimator should be one of: tmle, iptw, gcomp")
  if (estimator == "tmle") {
    if (object$gcomp) stop("estimator 'tmle' is not available because ltmleMSM was called with gcomp=TRUE")
    estimate <- object$beta
    IC <- object$IC
  } else if (estimator == "iptw") {
    estimate <- object$beta.iptw
    IC <- object$IC.iptw
  } else if (estimator == "gcomp") {
    if (!object$gcomp) stop("estimator 'gcomp' is not available because ltmleMSM was called with gcomp=FALSE")
    estimate <- object$beta
    IC <- object$IC
  }
  IC.variance <- apply(IC, 2, var)
  if (is.null(object$variance.estimate)) { 
    v <- IC.variance 
  } else {
    v <- pmax(diag(object$variance.estimate), IC.variance)
  }
  variance.estimate.ratio <- v / IC.variance
  return(list(estimate=estimate, IC=IC, variance.estimate.ratio=variance.estimate.ratio, v=v))
}

# Get summary measures for MSM parameters (standard errors, p-values, confidence intervals)
#' @rdname summary.ltmle
#' @export 
summary.ltmleMSM <- function(object, estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...) {
  info <- GetSummaryLtmleMSMInfo(object, estimator)
  estimate <- info$estimate
  v <- info$v
  n <- nrow(info$IC)
  std.dev <- sqrt(v/n)
  pval <- 2 * pnorm(-abs(estimate / std.dev))
  CI <- GetCI(estimate, std.dev)  
  cmat <- cbind(estimate, std.dev, CI, pval)
  dimnames(cmat) <- list(names(estimate), c("Estimate", "Std. Error", "CI 2.5%", "CI 97.5%", "p-value"))
  ans <- list(cmat=cmat, estimator=estimator, transformOutcome=object$transformOutcome, variance.estimate.ratio=info$variance.estimate.ratio) 
  class(ans) <- "summary.ltmleMSM"
  return(ans)
}

# Print method for summary.ltmleMSM
#' @rdname summary.ltmle
#' @export 
print.summary.ltmleMSM <- function(x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), ...) {
  cat("Estimator: ", x$estimator, "\n")
  if (x$estimator=="gcomp") {cat("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.\n")}
  printCoefmat(x$cmat, digits = digits, signif.stars = signif.stars, 
               na.print = "NA", has.Pvalue=TRUE, ...)
  if (x$transformOutcome) {
    Yrange <- attr(x$transformOutcome, "Yrange")
    cat("NOTE: The MSM is modeling the transformed outcome ( Y -", min(Yrange),
        ")/(", max(Yrange),"-", min(Yrange),")")
  }
  CheckVarianceEstimateRatio(x)
  invisible(x)
}

# Print method for summary.ltmle
#' @rdname summary.ltmle
#' @export 
print.summary.ltmle <- function(x, ...) {
  cat("Estimator: ", x$estimator, "\n")
  if (x$estimator=="gcomp") {cat("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.\n")}
  PrintCall(x$call)
  PrintSummary(x$treatment)
  CheckVarianceEstimateRatio(x)
  invisible(x)
}

# Print method for ltmleEffectMeasures
#' @rdname summary.ltmle
#' @export 
print.ltmleEffectMeasures <- function(x, ...) {
  PrintCall(x$call)
  cat("Use summary(...) to get estimates, standard errors, p-values, and confidence intervals for treatment EYd, control EYd, additive effect, relative risk, and odds ratio.\n")
  invisible(x)
}

# Print method for summary.ltmleEffectMeasures
#' @rdname summary.ltmle
#' @export 
print.summary.ltmleEffectMeasures <- function(x, ...) {
  cat("Estimator: ", x$estimator, "\n")
  if (x$estimator=="gcomp") {cat("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.\n")}
  PrintCall(x$call)
  lapply(x$effect.measures, PrintSummary)
  CheckVarianceEstimateRatio(x)
  if (x$transformOutcome) {
    Yrange <- attr(x$transformOutcome, "Yrange")
    cat("NOTE: All parameters are based on the transformed outcome ( Y -", min(Yrange),
        ")/(", max(Yrange),"-", min(Yrange),")")
  }
  invisible(x)
}

# Print a warning message if the TMLE based variance estimate is much greater than the IC based variance estimate 
CheckVarianceEstimateRatio <- function(summary.obj) {
  if (any(is.na(summary.obj$variance.estimate.ratio))) {
    warning("Unable to compute standard errors.")
    return(NULL)
  }
  if (any(summary.obj$variance.estimate.ratio > 100)) {
    warning.msg <- paste0("max(TMLE based variance estimate / IC based variance estimate) = ", floor(max(summary.obj$variance.estimate.ratio)), ".\nWhen this ratio is greater than 100, both variance estimates are less likely to be accurate.")
    warning(warning.msg)
  }
}

# Print method for ltmleMSM
#' @rdname summary.ltmle
#' @export 
print.ltmleMSM <- function(x, ...) {
  PrintCall(x$call)
  if (x$gcomp) {
    cat("GCOMP Beta Estimates: \n")
  } else {
    cat("TMLE Beta Estimates: \n")
  }
  print(x$beta)
  if (x$transformOutcome) {
    Yrange <- attr(x$transformOutcome, "Yrange")
    cat("NOTE: The MSM is modeling the transformed outcome ( Y -", min(Yrange),
        ")/(", max(Yrange),"-", min(Yrange),")")
  }
  invisible(x)
}

# Print method for ltmle
#' @rdname summary.ltmle
#' @export 
print.ltmle <- function(x, ...) {
  PrintCall(x$call)
  if (x$gcomp) {
    cat("GCOMP Estimate: ", x$estimates["gcomp"], "\n")
  } else {
    cat("TMLE Estimate: ", x$estimates["tmle"], "\n")
  }  
  invisible(x)
}

# Print a call
PrintCall <- function(cl) {
  cat("Call:\n", paste(deparse(cl), sep = "\n", collapse = "\n"), "\n\n", sep = "")
}

# Print estimate, standard error, p-value, confidence interval
PrintSummary <- function(x) {
  if (!is.null(x$long.name)) cat(x$long.name, ":\n", sep="")
  cat("   Parameter Estimate: ", signif(x$estimate, 5), "\n")
  if (x$log.std.err) {
    if (x$long.name == "Relative Risk") {
      param.abbrev <- "RR"
    } else if (x$long.name == "Odds Ratio") {
      param.abbrev <- "OR"
    } else {
      stop("unexpected x$long.name") 
    }
    cat("  Est Std Err log(", param.abbrev, "):  ", sep="")
  } else {
    cat("    Estimated Std Err:  ")
  }
  cat(signif(x$std.dev, 5), "\n")
  cat("              p-value: ", ifelse(x$pvalue <= 2*10^-16, "<2e-16",signif(x$pvalue, 5)), "\n")
  cat("    95% Conf Interval:",paste("(", signif(x$CI[1], 5), ", ", signif(x$CI[2], 5), ")", sep=""),"\n\n")
  invisible(x)
}

#Calculate estimate, standard deviation, p-value, confidence interval
GetSummary <- function(eff.list, cov.mat, n) {
  estimate <- eff.list$est
  v <- t(eff.list$gradient) %*% cov.mat %*% eff.list$gradient
  stopifnot(length(v) == 1)
  std.dev <- sqrt(v[1, 1] / n)
  
  if (eff.list$log.std.err) {
    pvalue <- 2 * pnorm(-abs(log(estimate) / std.dev))
    CI <- exp(GetCI(log(estimate), std.dev))
  } else {
    pvalue <- 2 * pnorm(-abs(estimate / std.dev))
    CI <- GetCI(estimate, std.dev)
  }
  CI <- Bound(CI, eff.list$CIBounds) 
  return(list(long.name=eff.list$long.name, estimate=estimate, std.dev=std.dev, pvalue=pvalue, CI=CI, log.std.err=eff.list$log.std.err))
}

# Calculate 95% confidence interval
GetCI <- function(estimate, std.dev) {
  x <- qnorm(0.975) * std.dev
  CI <- cbind("2.5%"=estimate - x, "97.5%"=estimate + x)
  return(CI)
}

# Parametric estimation of each g-factor
EstimateG <- function(inputs, regime.index) {
  abar <- GetABar(inputs$regimes, regime.index)
  gmat <- prob.A.is.1 <- matrix(NaN, nrow=nrow(inputs$data), ncol=length(inputs$nodes$AC))
  gmat.meanL <- prob.A.is.1.meanL <- cum.g.meanL <- cum.g.meanL.unbounded <- array(NaN, dim=c(nrow(inputs$data), length(inputs$nodes$AC), length(inputs$nodes$LY) - 1))
  uncensored <- rep(TRUE, nrow(inputs$data))
  fit <- vector("list", length(inputs$nodes$AC))
  names(fit) <- names(inputs$data)[inputs$nodes$AC]
  abar.meanL <- abar
  for (i in sseq(1, length(inputs$nodes$A))) {
    if (any(is.na(abar.meanL[, i]))) { #abar.meanL needs to be nonNA
      abar.meanL[is.na(abar.meanL[, i]), i] <- Mode(abar.meanL[, i], na.rm = TRUE)
    }
  }
  for (i in 1:length(inputs$nodes$AC)) {
    cur.node <- inputs$nodes$AC[i]
    uncensored <- IsUncensored(inputs$data, inputs$nodes$C, cur.node)
    newdata <- SetA(inputs$data, abar, inputs$nodes, cur.node)
    newdata.meanL <- SetA(inputs$data, abar.meanL, inputs$nodes, cur.node)
    deterministic.origdata <- IsDeterministic(inputs$data, cur.node, inputs$deterministic.Q.function, inputs$nodes, called.from.estimate.g=TRUE, inputs$survivalOutcome)$is.deterministic #deterministic due to death or Q.function
    deterministic.newdata <- IsDeterministic(newdata, cur.node, inputs$deterministic.Q.function, inputs$nodes, called.from.estimate.g=TRUE, inputs$survivalOutcome)$is.deterministic #deterministic due to death or Q.function - using data modified so A = abar
    if (is.numeric(inputs$gform)) {
      prob.A.is.1[, i] <- inputs$gform[, i, regime.index]  #if gform is numeric, it's a matrix of prob.A.is.1
      g.est <- list(fit="gform passed as numeric, so no estimation took place")
    } else {
      deterministic.g.list.origdata <- IsDeterministicG(inputs$data, cur.node, inputs$deterministic.g.function, inputs$nodes, using.newdata=F) #deterministic due to acnode map - using original data
      deterministic.g.list.newdata <- IsDeterministicG(newdata, cur.node, inputs$deterministic.g.function, inputs$nodes, using.newdata=T) #deterministic due to acnode map - using data modified so A = abar
      deterministic.g.origdata <- deterministic.g.list.origdata$is.deterministic
      if (inputs$stratify) {
        intervention.match <- InterventionMatch(inputs$data, abar, inputs$nodes$A, inputs$nodes$AC[i]) 
        subs <- uncensored & intervention.match & !deterministic.origdata & !deterministic.g.origdata
      } else {
        subs <- uncensored & !deterministic.origdata & !deterministic.g.origdata
      }
      
      if (all(deterministic.g.list.newdata$is.deterministic | deterministic.newdata)) {
        # all rows are set deterministically, no need to estimate
        g.est <- list(fit="all rows are set deterministically, no estimation at this node")
        #fill in prob.A.is.1 below
      } else {
        # not all rows are set deterministically
        if (any(subs)) {
          g.est <- Estimate(inputs$gform[i], data=inputs$data, subs=subs, family="quasibinomial", newdata=newdata, SL.library=inputs$SL.library.g, type="response", nodes=inputs$nodes, observation.weights=inputs$observation.weights)
          prob.A.is.1[, i] <- g.est$predicted.values
          #n x numACnodes x (numLYnodes - 1)
          #[,,k] is prob.A.is.1 with all L and Y nodes after and including LYnodes[k] set to mean of L (na.rm=T)
          prob.A.is.1.meanL[, i, ] <- PredictProbAMeanL(newdata.meanL, inputs$nodes, subs, g.est$fit, SL.XY=g.est$SL.XY)
        } else {
          msg <- paste0("ltmle failed trying to estimate ", inputs$gform[i], " because there are no observations that are\nuncensored", ifelse(inputs$stratify, ", follow abar,", ""), " and are not set deterministically due to death or deterministic.g.function or deterministic.Q.function\n")
          stop(msg)
        }
      }
      prob.A.is.1[deterministic.g.list.newdata$is.deterministic, i] <- prob.A.is.1.meanL[deterministic.g.list.newdata$is.deterministic, i, ] <- deterministic.g.list.newdata$prob1
    } 
    #prob.A.is.1 is prob(a=1), gmat is prob(a=abar)
    #cur.abar can be NA after censoring/death if treatment is dynamic
    if (cur.node %in% inputs$nodes$A) {
      cur.abar <- abar[, inputs$nodes$A == cur.node]
      cur.abar.meanL <- abar.meanL[, inputs$nodes$A == cur.node]
    } else {
      cur.abar <- cur.abar.meanL <- rep(1, nrow(inputs$data))  #if this is a cnode, abar is always 1 (uncensored)
    }
    gmat[, i] <- CalcGVec(prob.A.is.1[, i], cur.abar, deterministic.newdata)
    gmat.meanL[, i, ] <- apply(AsMatrix(prob.A.is.1.meanL[, i, ]), 2, CalcGVec, cur.abar.meanL, deterministic.newdata)
    if (any(is.na(gmat[uncensored, i]))) stop("Error - NA in g. g should only be NA after censoring. If you passed numeric gform, make sure there are no NA values except after censoring. Otherwise something has gone wrong.")
    fit[[i]] <- g.est$fit
  }
  cum.g.unbounded <- CalcCumG(gmat, c(0, 1))
  cum.g <- CalcCumG(gmat, inputs$gbounds)
  for (i in sseq(1, dim(gmat.meanL)[3])) {
    cum.g.meanL[, , i] <- CalcCumG(drop3(gmat.meanL[, , i, drop=F]), inputs$gbounds)
    cum.g.meanL.unbounded[, , i] <- CalcCumG(drop3(gmat.meanL[, , i, drop=F]), c(0,1)) 
  }
  return(list(cum.g=cum.g, cum.g.unbounded=cum.g.unbounded, cum.g.meanL=cum.g.meanL, fit=fit, prob.A.is.1=prob.A.is.1, cum.g.meanL.unbounded=cum.g.meanL.unbounded))
}

PredictProbAMeanL <- function(data, nodes, subs, fit, SL.XY) {
  #Predict the probability that A=1 if L and Y nodes are set to their mean (or median) values
  
  #probAis1.meanL is n x num.LYnodes - 1
  #probAis1.meanL[, k] is prob.A.is.1 with all L and Y nodes after and including LYnodes[k] set to mean of L
  
  #somewhat inefficient - for W A.1 L.2 A.2 L.3 A.3 Y, does P(A.1=1) setting L.3 to mean and then L.2 and L.3 to mean, but none of these can be used in P(A.1=1) because they're after A.1
  
  #A is already set to abar in data
  probAis1.meanL <- matrix(NaN, nrow(data), length(nodes$LY) - 1)
  if (ncol(probAis1.meanL) == 0) return(probAis1.meanL)
  all.LY.nodes <- sort(union(nodes$L, nodes$Y)) #not the same as nodes$LY, which removes blocks
  newdata <- data
  LYindex <- length(nodes$LY)
  for (i in length(all.LY.nodes):1) { 
    regression.node <- all.LY.nodes[i]
    L <- data[subs, regression.node]
    if (is.numeric(L) && !IsBinary(L)) {
      meanL <- mean(L, na.rm = TRUE)
    } else {
      meanL <- Mode(L, na.rm = TRUE) #for factors and binaries
    }
    newdata[, regression.node] <- meanL
    if (regression.node %in% nodes$LY[1:length(nodes$LY)-1]) {
      LYindex <- LYindex - 1
      SuppressGivenWarnings({
        if ("SuperLearner" %in% class(fit)) {
          if ("ltmle.added.constant" %in% SL.XY$rhs) {   #see ltmle:::Estimate for why this is needed
            newdata.temp <- cbind(newdata, ltmle.added.constant=1)
          } else {
            newdata.temp <- newdata
          }
          newdata.temp <- newdata.temp[, SL.XY$rhs, drop=FALSE]
          probAis1.meanL[, LYindex] <- predict(fit, newdata = newdata.temp, X=SL.XY$X, Y=SL.XY$Y, onlySL=TRUE)$pred
        } else {
          probAis1.meanL[, LYindex] <- predict(fit, newdata = newdata, type = "response")
        }
      }, "prediction from a rank-deficient fit may be misleading")
    }
  }
  if (any(is.na(probAis1.meanL[, 1]))) stop("NA in probAis1.meanL[, 1]")
  return(probAis1.meanL)
}



CalcGVec <- function(prob.A.is.1, cur.abar, deterministic.newdata) {
  g <- rep(NA, length(prob.A.is.1))
  g[!is.na(cur.abar) & cur.abar == 1] <- prob.A.is.1[!is.na(cur.abar) & cur.abar == 1]
  g[!is.na(cur.abar) & cur.abar == 0] <- 1 - prob.A.is.1[!is.na(cur.abar) & cur.abar == 0]    
  g[deterministic.newdata] <- 1  #a=abar deterministically after death or other deterministic Q
  return(g)
}

# Truncate values within supplied bounds
Bound <- function(x, bounds) {
  stopifnot(length(bounds) == 2 && !any(is.na(bounds)))
  x[x < min(bounds)] <- min(bounds)
  x[x > max(bounds)] <- max(bounds)
  return(x)
}

# Convert named nodes to indicies of nodes
NodeToIndex <- function(data, node) {
  if (! is.data.frame(data)) stop("data must be a data frame")
  if (is.numeric(node) || is.null(node)) return(node)
  if (! is.character(node)) stop("nodes must be numeric, character, or NULL")
  index <- match(node, names(data))
  if (any(is.na(index))) {
    stop(paste("\nnamed node(s) not found:", node[is.na(index)]))
  }
  return(index)
}

# Run GLM or SuperLearner
Estimate <- function(form, data, subs, family, newdata, SL.library, type, nodes, observation.weights) {
  stopifnot(type %in% c("link", "response"))
  if (form == "IDENTITY") {
    predicted.values <- data[, "Q.kplus1"]
    if (type == "link") {
      stopifnot(family %in% c("binomial", "quasibinomial"))
      predicted.values <- qlogis(Bound(predicted.values, bounds=c(0.0001, 0.9999)))
    }
    m <- "no fit because form == IDENTITY"
    return(list(predicted.values=predicted.values, fit=m))
  }
  data <- ConvertCensoringNodesToBinary(data, nodes$C) #convert factors to binaries for compatability with glm and some SL libraries
  f <- as.formula(form)
  if (any(is.na(data[subs, LhsVars(f)]))) stop("NA in Estimate")
  observation.weights <- observation.weights[subs]
  if (is.null(SL.library) || length(RhsVars(f)) == 0) { #in a formula like "Y ~ 1", call glm
    #estimate using GLM
    if (sum(subs) > 1) {
      SuppressGivenWarnings({
        m <- get.stack("glm.ltmle.memoized", mode="function", ifnotfound=glm.ltmle)(f, data=data[subs, all.vars(f), drop=F], observation.weights=observation.weights, family=family, control=glm.control(trace=FALSE, maxit=1000)) 
        predicted.values <- predict(m, newdata=newdata, type=type)
      }, GetWarningsToSuppress())
    } else {
      #glm breaks when sum(subs) == 1
      predicted.values <- rep(data[subs, LhsVars(f)], nrow(newdata))
      m <- "fit not returned because there was only 1 observation to fit"
    }
    SL.XY <- NULL
  } else {
    #estimate using SuperLearner
    if (family == "quasibinomial") family <- "binomial"
    
    rhs <- setdiff(RhsVars(f), rownames(alias(f, data=data[subs,])$Complete))  #remove aliased columns from X - these can cause problems if they contain NAs and the user is expecting the column to be dropped
    new.subs <- apply(newdata[, rhs, drop=FALSE], 1, function (x) !any(is.na(x)))  #remove NA values from newdata - these will output to NA anyway and cause errors in SuperLearner
    Y <- data[subs, LhsVars(f)]
    X <- data[subs, rhs, drop=FALSE]
    newX <- newdata[new.subs, rhs, drop=FALSE]
    if (ncol(X) == 1) {
      #SuperLearner crashes if there are screening algorithms and only one column - add a constant
      X <- cbind(X, ltmle.added.constant=1)
      newX <- cbind(newX, ltmle.added.constant=1)
      rhs <- c(rhs, "ltmle.added.constant")
    }
    SetSeedIfRegressionTesting()
    try.result <- try({
      SuppressGivenWarnings(m <- SuperLearner::SuperLearner(Y=Y, X=X, SL.library=SL.library, verbose=FALSE, family=family, newX=newX, obsWeights=observation.weights), c("non-integer #successes in a binomial glm!", "prediction from a rank-deficient fit may be misleading")) 
    })
    SL.XY <- list(X=X, Y=Y, rhs=rhs)
    GetSLStopMsg <- function(Y) ifelse(all(Y %in% c(0, 1, NA)), "", "\n Note that many SuperLeaner libraries crash when called with continuous dependent variables, as in the case of initial Q regressions with continuous Y or subsequent Q regressions even if Y is binary.")
    if (inherits(try.result, "try-error")) {
      stop(paste("\n\nError occured during call to SuperLearner:\n", form, GetSLStopMsg(Y), "\n The error reported is:\n", try.result))
    }
    if (all(is.na(m$SL.predict))) {
      stop(paste("\n\nSuperLearner returned all NAs during regression:\n", form, GetSLStopMsg(Y)))
    }
    predicted.values <- rep(NA, nrow(newdata))
    predicted.values[new.subs] <- m$SL.predict
    if (max(predicted.values, na.rm=T) > 1 || min(predicted.values, na.rm=T) < 0) {
      msg <- paste("SuperLearner returned predicted.values > 1 or < 0: [min, max] = [", min(predicted.values, na.rm=T), ",", max(predicted.values, na.rm=T), "]. Bounding to [0,1]")
      warning(msg)
      predicted.values <- Bound(predicted.values, bounds=c(0, 1))
    }
    if (type == "link") {
      stopifnot(family == "binomial")
      predicted.values <- qlogis(Bound(predicted.values, bounds=c(0.0001, 0.9999)))
    }
  }
  return(list(predicted.values=predicted.values, fit=m, SL.XY=SL.XY))
}

# This is here for memoizing
glm.ltmle <- function(f, data, family, control, observation.weights) {
  #note: observation.weights is in Estimate (the environment of f)
  return(glm(f, data=data.frame(data[, all.vars(f), drop=F], observation.weights), family=family, control=control, weights=scale(observation.weights, center=FALSE)))
}

# Calculate bounded cumulative G
CalcCumG <- function(g, gbounds) {
  cum.g <- AsMatrix(Bound(t(apply(g, 1, cumprod)), gbounds)) #AsMatrix to fix problems where apply returns a vector
  return(cum.g)
}

# Determine which patients are following specified treatment regime (abar)
#return vector of [numObservations x 1] I(A==abar) from Anodes[1] to the Anode just before cur.node
# note: if calling from outside ltmle:::, cur.node needs to be the node index, not a string!
InterventionMatch <- function(data, abar, Anodes, cur.node) {
  intervention.match <- XMatch(data, abar, Anodes, cur.node, all, default=TRUE)
  return(intervention.match)
}

# Determine which patients are uncensored
#return vector of [numDataRows x 1] I(C=uncensored) from Cnodes[1] to the Cnode just before cur.node
# note: if calling from outside ltmle:::, cur.node needs to be the node index, not a string!
IsUncensored <- function(data, Cnodes, cur.node) {
  if (! all(sapply(data[, Cnodes], is.factor))) stop("something has gone wrong in ltmle:::IsUncensored - all Cnodes should have been converted to factors")
  uncensored <- XMatch(data, Xbar="uncensored", Cnodes, cur.node, all, default=TRUE)
  return(uncensored)
}

# Determine which patients have died or have Q set deterministically by user function before cur.node
# return list:
#    is.deterministic: vector of [numObservations x 1] - true if patient is already dead before cur.node or set by deterministic.Q.function
#    Q.value: vector of [which(is.deterministic) x 1] - value of Q
IsDeterministic <- function(data, cur.node, deterministic.Q.function, nodes, called.from.estimate.g, survivalOutcome) {
  #set Q.value to 1 if previous y node is 1
  if (survivalOutcome) {
    is.deterministic <- XMatch(data, Xbar=1, nodes$Y, cur.node, any, default=FALSE) #deterministic if any previous y node is 1
  } else {
    is.deterministic <- rep(FALSE, nrow(data))
  }
  
  #get Q values from deterministic.Q.function
  default <- list(is.deterministic=is.deterministic, Q.value=1)
  if (is.null(deterministic.Q.function)) return(default)
  #put this in a try-catch?
  det.list <- deterministic.Q.function(data=data, current.node=cur.node, nodes=nodes, called.from.estimate.g=called.from.estimate.g)
  if (is.null(det.list)) return(default)
  if (called.from.estimate.g) {
    #it's ok if Q.value isn't returned if called.from.estimate.g
    if (!is.list(det.list) || !("is.deterministic" %in% names(det.list)) || !(length(det.list) %in% 1:2)) stop("deterministic.Q.function should return a list with names: is.deterministic, Q.value")
  } else {
    if (!is.list(det.list) || !setequal(names(det.list), c("is.deterministic", "Q.value")) || length(det.list) != 2) stop("deterministic.Q.function should return a list with names: is.deterministic, Q.value")
  }
  
  if (! length(det.list$Q.value) %in% c(1, length(which(det.list$is.deterministic)))) stop("the length of the 'Q.value' element of deterministic.Q.function's return argument should be either 1 or length(which(det.list$is.deterministic))")
  
  #check that these observations where Q.value is 1 due to death (previous y is 1) aren't set to anything conflicting by deterministic.Q.function
  Q.value.from.function <- rep(NA, nrow(data))
  Q.value.from.function[det.list$is.deterministic] <- det.list$Q.value
  set.by.function.and.death <- is.deterministic & det.list$is.deterministic
  if (any(Q.value.from.function[set.by.function.and.death] != 1)) {
    stop(paste("inconsistent deterministic Q at node:", names(data)[cur.node]))
  }
  finalY <- data[, max(nodes$Y)]
  inconsistent.rows <- (det.list$Q.value %in% c(0,1)) & (det.list$Q.value != finalY[det.list$is.deterministic]) & !is.na(finalY[det.list$is.deterministic])
  if (any(inconsistent.rows)) stop(paste("At node:",names(data)[cur.node], "deterministic.Q.function is inconsistent with data - Q.value is either 0 or 1 but this does not match the final Y node value\nCheck data rows:", paste(head(rownames(data)[det.list$is.deterministic][inconsistent.rows]), collapse=" ")))
  
  #return combined values
  Q.value <- rep(NA, nrow(data))
  Q.value[is.deterministic] <- 1
  Q.value[det.list$is.deterministic] <- det.list$Q.value
  is.deterministic <- is.deterministic | det.list$is.deterministic
  Q.value <- Q.value[is.deterministic]
  if (any(is.na(c(is.deterministic, Q.value)))) stop("NA in is.deterministic or Q.value")
  return(list(is.deterministic=is.deterministic, Q.value=Q.value))
}

# Determine which patients have an Anode value which is deterministic 
# For example, deterministic.g.function may be used to specify that once a patient starts treatment, they stay on treatment and this should be taken into consideration during estimation of G
IsDeterministicG <- function(data, cur.node, deterministic.g.function, nodes, using.newdata) {
  default <- list(is.deterministic=rep(FALSE, nrow(data)), prob1=NULL)
  if (is.null(deterministic.g.function)) return(default)
  #put this in a try-catch?
  det.list <- deterministic.g.function(data=data, current.node=cur.node, nodes=nodes)
  if (is.null(det.list)) return(default)
  if (!is.list(det.list) || !setequal(names(det.list), c("is.deterministic", "prob1")) || length(det.list) != 2) stop("deterministic.g.function should return a list with names: is.deterministic, prob1")
  if (! length(det.list$prob1) %in% c(1, length(which(det.list$is.deterministic)))) stop("the length of the 'prob1' element of deterministic.g.function's return argument should be either 1 or length(which(det.list$is.deterministic))")
  
  inconsistent.rows <- (det.list$prob1 %in% c(0,1)) & (det.list$prob1 != data[det.list$is.deterministic, cur.node]) & !is.na(data[det.list$is.deterministic, cur.node])
  if (any(inconsistent.rows)) {
    err.msg <- paste("At node:",names(data)[cur.node], "deterministic.g.function is inconsistent with data - prob1 is either 0 or 1 but this does not match the node value.\nCheck data rows:", paste(head(rownames(data)[det.list$is.deterministic][inconsistent.rows]), collapse=" "))
    if (using.newdata) {
      err.msg <- paste(err.msg, "\n This error occured while calling deterministic.g.function on data where Anodes are set to abar.")
      cat("deterministic.g.function is inconsistent with data.\nAfter setting Anodes to abar, the data looks like this:\n")
      print(head(data[det.list$is.deterministic[inconsistent.rows], ]))
    }
    stop(err.msg)
  }
  return(det.list)
}

#Utility function called by IsUncensored, IsDeterministic - compares history in d within Xnodes prior to cur.node to Xbar
#
#any.all should be either the function 'any' or the function 'all'
#default: value to return if value is NA or there are no nodes before cur.node (TRUE for InterventionMatch and IsUncensored because NA indicates person matches intervention/is uncensored until they died; FALSE for IsDeterministic because NA indicates person was alive until they were censored)
XMatch <- function(data, Xbar, Xnodes, cur.node, any.all, default) {
  if (!any(Xnodes < cur.node)) return(rep(default, nrow(data)))
  last.Xnode.index <- which.max(Xnodes[Xnodes < cur.node])
  
  Xnodes.subset <- Xnodes[1:last.Xnode.index]
  if (identical(Xbar, 1) || identical(Xbar, "uncensored")) {
    Xbar.subset <- Xbar 
  } else {
    Xbar.subset <- Xbar[, 1:last.Xnode.index]
  } 
  d.subset <- data[, Xnodes.subset, drop=FALSE]
  matches <- apply(d.subset == Xbar.subset, 1, any.all)
  matches[is.na(matches)] <- default
  return(matches)
}

# Calculate the TMLE influence curve for one node
CalcIC <- function(Qstar.kplus1, Qstar, h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight) {
  n <- nrow(Qstar)
  num.regimes <- ncol(Qstar)
  num.betas <- dim(h.g.ratio)[3] #h.g.ratio: n x num.regimes x num.betas
  
  IC <- matrix(0, nrow=n, ncol=num.betas)
  for (i in regimes.with.positive.weight) {
    index <- uncensored & intervention.match[, i]
    if (any(h.g.ratio[index, i, ] != 0)) {
      regimeIC <- matrix(0, nrow=n, ncol=num.betas)
      regimeIC[index, ] <- (Qstar.kplus1[index, i] - Qstar[index, i]) * h.g.ratio[index, i, ]
      IC <- IC + regimeIC
    }
  }
  return(IC)
}

#Set the Anodes of d to abar and Cnodes to uncensored (up to and including cur.node - cur.node itself is included for consistency checking in DeterministicG)
SetA <- function(data, abar, nodes, cur.node) {
  Anode.index <- nodes$A <= cur.node
  data[, nodes$A[Anode.index]] <- abar[, Anode.index]
  
  Cnode.index <- nodes$C <= cur.node
  data[, nodes$C[Cnode.index]] <- factor(rep("uncensored", nrow(data))) #recycled
  return(data)
}

# Return the left hand side variable of formula f as a character
LhsVars <- function(f) {
  f <- as.formula(f)
  return(as.character(f[[2]]))
}

# Return the right hand side variables of formula f as a character vector
RhsVars <- function(f) {
  f <- as.formula(f)
  return(all.vars(f[[3]]))
}

# Error checking for inputs
CheckInputs <- function(data, nodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, deterministic.Q.function, observation.weights, gcomp) {
  stopifnot(length(dim(regimes)) == 3)
  num.regimes <- dim(regimes)[3]
  if (!all(is.null(GetLibrary(SL.library, "Q")), is.null(GetLibrary(SL.library, "g")))) {
    if (!requireNamespace("SuperLearner")) stop("SuperLearner package is required if SL.library is not NULL")
  } 
  #each set of nodes should be sorted - otherwise causes confusion with gform, Qform, abar
  if (is.unsorted(nodes$A, strictly=TRUE)) stop("Anodes must be in increasing order")
  if (is.unsorted(nodes$C, strictly=TRUE)) stop("Cnodes must be in increasing order")
  if (is.unsorted(nodes$L, strictly=TRUE)) stop("Lnodes must be in increasing order")
  if (is.unsorted(nodes$Y, strictly=TRUE)) stop("Ynodes must be in increasing order")
  if (is.unsorted(final.Ynodes, strictly=TRUE)) stop("final.Ynodes must be in increasing order")
  
  if (length(nodes$L) > 0) {
    if (min(nodes$L) < min(nodes$AC)) stop("Lnodes are not allowed before A/C nodes. If you want to include baseline nodes, include them in data but not in Lnodes")
    if (max(nodes$L) > max(nodes$Y)) stop("Lnodes are not allowed after the final Y node")
  }
  if (min(nodes$Y) < min(nodes$AC)) stop("Ynodes are not currently allowed before A/C nodes.")
  
  all.nodes <- c(nodes$A, nodes$C, nodes$L, nodes$Y)
  if (length(all.nodes) > length(unique(all.nodes))) stop("A node cannot be listed in more than one of Anodes, Cnodes, Lnodes, Ynodes")
  if (is.null(nodes$Y)) stop("Ynodes cannot be null")
  if (is.null(nodes$AC)) stop("Anodes and Cnodes cannot both be null")
  
  if (min(all.nodes) < ncol(data)) {
    if (!all((min(all.nodes):ncol(data)) %in% all.nodes)) {
      stop("All nodes after the first of A-, C-, L-, or Ynodes must be in A-, C-, L-, or Ynodes")
    }
  }
  for (reserved.name in c("observation.weights", "Q.kplus1", "Qstar")) {
    #these might cause conflicts
    if (reserved.name %in% names(data)) stop(paste(reserved.name, "is reserved and may not be used as a column name of data"))
  }
  
  #If gform is NULL, it will be set by GetDefaultForm; no need to check here
  if (!is.null(gform)) {
    if (is.character(gform)) {
      if (length(gform) != length(nodes$AC)) stop("length(gform) != length(c(Anodes, Cnodes))")
      for (i in 1:length(gform)) {
        if (LhsVars(gform[i]) != names(data)[nodes$AC[i]]) {
          stop("The LHS of gform[", i, "] should be the name of the ", i, "th A or C node")
        }
        parents <- if(nodes$AC[i] > 1) {
          names(data)[1:(nodes$AC[i]-1)]
        } else {
          NULL
        }
        if (!all(RhsVars(gform[i]) %in% parents)) {
          stop("Some nodes in gform[", i, "] are not parents of ", LhsVars(gform[i]))
        }
      }
    } else {
      if (! is.numeric(gform)) stop("gform should be a character vector or numeric")
      if (nrow(gform) != nrow(data)) stop("if gform is numeric, it should have the same number of rows as data")
      if (ncol(gform) != length(nodes$AC)) stop("if gform is numeric, it should have the same number of columns as length(c(Anodes, Cnodes))")
      if (length(dim(gform)) != 3 || dim(gform)[3] != num.regimes) stop("if gform is numeric, dim[3] should be num.regimes")
      if (!is.null(deterministic.g.function)) stop("if gform is numeric, deterministic.g.function must be NULL")
      if (max(gform, na.rm=T) > 1 || min(gform, na.rm=T) < 0) stop("if gform is numeric, all values should be probabilities")
    }
  }
  
  #If Qform is NULL, it will be set by GetDefaultForm; no need to check here
  if (!is.null(Qform)) {
    if (! is.character(Qform)) stop("Qform should be a character vector")
    if (length(Qform) != length(nodes$LY)) {
      stop("length of Qform is not equal to number of L/Y nodes")
    }
    for (i in 1:length(Qform)) {
      if (length(names(Qform[i])) == 0) stop("Each element of Qform must be named. The name must match the name of the corresponding L/Y node in data.")
      if (names(Qform[i]) != names(data)[nodes$LY[i]]) stop("The name of each element of Q must match the name of the corresponding L/Y node in data.")
      if (Qform[i] != "IDENTITY") {
        #This is only meant to be used in a call by ltmle:::EstimateVariance
        if (LhsVars(Qform[i]) != "Q.kplus1") stop("LHS of each Qform should be Q.kplus1")
        parents <- names(data)[1:(nodes$LY[i]-1)]
        if (!all(RhsVars(Qform[i]) %in% parents)) {
          stop("Some nodes in Qform[", i, "] are not parents of ", names(Qform[i]))
        }    
      }
    }
  }
  
  if (length(gbounds) != 2) stop("gbounds should have length 2")
  if (! (is.null(deterministic.g.function) || is.function(deterministic.g.function))) {
    stop("deterministic.g.function should be a function or NULL")
  }
  
  if (! all(unlist(data[, nodes$A]) %in% c(0, 1, NA))) stop("in data, all Anodes should be binary")
  #note: Cnodes are checked in ConvertCensoringNodes
  
  all.Y <- unlist(data[, nodes$Y])
  
  binaryOutcome <- all(all.Y %in% c(0, 1, NA))
  
  if (binaryOutcome) {
    if (is.null(survivalOutcome)) {
      if (length(nodes$Y) == 1) {
        survivalOutcome <- FALSE #doesn't matter 
      } else {
        stop("All Ynodes are 0, 1, or NA; the outcome is treated as binary. The 'survivalOutcome' argument must be specified if there are multiple Ynodes.")
      }
    }    
    if (!is.null(Yrange) && !is.equal(Yrange, c(0L, 1L))) {
      stop("All Ynodes are 0, 1, or NA, but Yrange is something other than NULL or c(0, 1)")
    }
  } else {
    if (is.null(survivalOutcome)) {
      survivalOutcome <- FALSE
    }    
    if (survivalOutcome) {
      stop("When survivalOutcome is TRUE, all Ynodes should be 0, 1, or NA")
    } 
  }
  
  for (i in nodes$Y) {
    uncensored <- IsUncensored(data, nodes$C, cur.node=i)
    deterministic <- IsDeterministic(data, cur.node=i, deterministic.Q.function=NULL, nodes, called.from.estimate.g=FALSE, survivalOutcome)$is.deterministic #pass deterministic.Q.function=NULL so we're only picking up deaths (if surivalOutcome=FALSE, deterministic will all be FALSE)
    if (any(is.na(data[deterministic, i])) || ! all(data[deterministic, i] == 1)) stop("For survival outcomes, once a Ynode jumps to 1 (e.g. death), all subsequent Ynode values should be 1.")    
    if (any(is.na(data[uncensored, i]))) stop("Ynodes may not be NA except after censoring")
  }
  
  if (! is.equal(dim(regimes)[1:2], c(nrow(data), length(nodes$A)))) stop("Problem with abar or regimes:\n   In ltmleMSM, regimes should have dimensions n x num.Anodes x num.regimes\n   In ltmle, abar should be a matrix with dimensions n x num.Anodes or a vector with length num.Anodes")
  stopifnot(num.regimes == nrow(summary.measures))
  if (!all(regimes %in% c(0, 1, NA))) stop("all regimes should be binary")
  for (i in seq_along(nodes$A)) {
    cur.node <- nodes$A[i]
    uncensored <- IsUncensored(data, nodes$C, cur.node)
    deterministic <- IsDeterministic(data, cur.node, deterministic.Q.function, nodes, called.from.estimate.g=TRUE, survivalOutcome)$is.deterministic
    if (any(is.na(regimes[uncensored & !deterministic, i, ]))) {
      stop("NA in regimes/abar not allowed (except after censoring/death)")
    }
  }
  
  if ((length(dim(summary.measures)) != 3) || ! is.equal(dim(summary.measures)[c(1, 3)], c(num.regimes, length(final.Ynodes)))) stop("summary.measures should be an array with dimensions num.regimes x num.summary.measures x num.final.Ynodes")
  if (class(working.msm) != "character") stop("class(working.msm) must be 'character'")
  if (LhsVars(working.msm) != "Y") stop("the left hand side variable of working.msm should always be 'Y' [this may change in future releases]")
  if (!is.vector(observation.weights) || length(observation.weights) != nrow(data) || any(is.na(observation.weights)) || any(observation.weights < 0) || max(observation.weights) == 0) stop("observation.weights must be NULL or a vector of length nrow(data) with no NAs, no negative values, and at least one positive value")
 return(list(survivalOutcome=survivalOutcome, binaryOutcome=binaryOutcome))
}

TransformOutcomes <- function(data, nodes, Yrange) {
  all.Y <- unlist(data[, nodes$Y])
  transformOutcome <- FALSE
  if (!is.null(Yrange)) {
    #if Yrange was specified
    rng <- range(all.Y, na.rm=TRUE)
    if (min(rng) < min(Yrange) || max(rng) > max(Yrange)) {
      #Truncate if Y vals are outside Yrange
      message("Some Ynodes are not in [Yrange[1], Yrange[2]], Y values are truncated")
      data[,nodes$Y][data[,nodes$Y] < min(Yrange)]<- min(Yrange)
      data[,nodes$Y][data[,nodes$Y] > max(Yrange)] <- max(Yrange)       
    } 
    #Then transform
    transformOutcome <- TRUE
  } else {
    #if Yrange was not specified, get it
    Yrange <- range(all.Y, na.rm=TRUE)
    if (min(Yrange) < 0 || max(Yrange) > 1) {
      #And see if we need to transform
      transformOutcome <- TRUE
      message("Some Ynodes are not in [0, 1], and Yrange was NULL, so all Y nodes are being\ntransformed to (Y-min.of.all.Ys)/range.of.all.Ys") 
    }
  }
  if (transformOutcome) {
    attr(transformOutcome, 'Yrange') <- Yrange 
    data[,nodes$Y] <- (data[, nodes$Y]-min(Yrange))/diff(Yrange) 
  }
  return(list(data=data, transformOutcome=transformOutcome))
}

# Set all nodes (except Y) to NA after death or censoring; Set Y nodes to 1 after death
CleanData <- function(data, nodes, deterministic.Q.function, survivalOutcome, showMessage=TRUE) {
  #make sure binaries have already been converted before calling this function
  is.nan.df <- function (x) {
    y <- if (length(x)) {
      do.call("cbind", lapply(x, "is.nan"))
    } else {
      matrix(FALSE, length(row.names(x)), 0)
    }
  }
  is.na.strict <- function (x) is.na(x) & !is.nan.df(x)  #only for data.frames
  changed <- FALSE
  ua <- rep(TRUE, nrow(data))  #uncensored and alive
  if (ncol(data) == 1) return(data)
  deterministic.Q.function.depends.on.called.from.estimate.g <- length(grep("called.from.estimate.g", as.character(body(deterministic.Q.function)))) > 0
  for (i in 1:(ncol(data)-1)) {
    if (any(is.na(data[ua, 1:i]))) stop("NA values are not permitted in data except after censoring or a survival event")
    is.deterministic <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=TRUE, survivalOutcome=survivalOutcome)$is.deterministic #check determinisitic including node i 
    
    if (deterministic.Q.function.depends.on.called.from.estimate.g) {
      is.deterministic.Q <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=FALSE, survivalOutcome=survivalOutcome)$is.deterministic 
      if (any(is.deterministic[ua] & !is.deterministic.Q[ua])) stop("Any row set deterministic by deterministic.Q.function(..., called.from.estimate.g=TRUE) must imply that the row is also set deterministic by deterministic.Q.function(..., called.from.estimate.g=FALSE)") #det.Q.fun(T) should imply det.Q.fun(F)
    }
    
    ua[ua] <- !is.deterministic[ua]
    if (any(is.na(ua))) stop("internal ltmle error - ua should not be NA in CleanData")
    if (! all(is.na.strict(data[is.deterministic, setdiff((i+1):ncol(data), nodes$Y), drop=FALSE]))) {
      data[is.deterministic, setdiff((i+1):ncol(data), nodes$Y)] <- NA #if deterministic, set all nodes except Y to NA
      changed <- TRUE
    }
    
    if (i %in% nodes$C) {
      censored <- data[, i] == "censored" & ua
      if (! all(is.na.strict(data[censored, (i+1):ncol(data), drop=FALSE]))) {
        data[censored, (i+1):ncol(data)] <- NA  #if censored, set all nodes (including Y) to NA
        changed <- TRUE
      }
      ua[ua] <- !censored[ua] 
      if (any(is.na(ua))) stop("internal ltmle error - ua should not be NA in CleanData")
    } 
  }
  if (changed && showMessage) {
    message("Note: for internal purposes, all nodes after a censoring event are set to NA and \n all nodes (except Ynodes) are set to NA after Y=1 if survivalFunction is TRUE (or if specified by deterministic.Q.function).\n Your data did not conform and has been adjusted. This may be relevant if you are \n writing your own deterministic function(s) or debugging ltmle.")
  }
  return(data)
}

# Get the default Q or g formula - each formula consists of all parent nodes except censoring and event nodes [also except A nodes if stratifying]
GetDefaultForm <- function(data, nodes, is.Qform, stratify, survivalOutcome, showMessage) {
  if (is.Qform) {
    lhs <- rep("Q.kplus1", length(nodes$LY))
    node.set <- nodes$LY
  } else {
    lhs <- names(data)[nodes$AC]
    node.set <- nodes$AC
  }
  if (stratify) {
    stratify.nodes <- c(nodes$C, nodes$A)
  } else {
    stratify.nodes <- c(nodes$C)
  }
  if (survivalOutcome) {
    stratify.nodes <- c(stratify.nodes, nodes$Y)
  }
  form <- NULL
  for (i in seq_along(node.set)) {
    cur.node <- node.set[i]
    if (cur.node == 1) {
      form[i] <- paste(lhs[i], "~ 1")  #no parent nodes
    } else {
      parent.node.names <- names(data)[setdiff(1:(cur.node - 1), stratify.nodes)]
      if (length(parent.node.names) == 0) {
        form[i] <- paste(lhs[i], "~ 1")
      } else {
        form[i] <- paste(lhs[i], "~", paste(parent.node.names, collapse=" + "))     
      }
    }
    names(form)[i] <- names(data)[cur.node]
  }
  
  if (showMessage) {
    #Prints formulas with automatic wrapping thanks to print.formula
    message(ifelse(is.Qform, "Qform", "gform"),
            " not specified, using defaults:")
    lapply(seq_along(form), function(i, names) {
      message("formula for ", names[i], ":")
      #Using print on a formula because it nicely wraps
      message(capture.output(print(as.formula(form[i]), showEnv=FALSE)))
    }, names=names(form))
    message("")
  }
  return(form)
}

# Organize nodes
CreateNodes <- function(data, Anodes, Cnodes, Lnodes, Ynodes) {  
  Anodes <- NodeToIndex(data, Anodes)
  Cnodes <- NodeToIndex(data, Cnodes)
  Lnodes <- NodeToIndex(data, Lnodes)
  Ynodes <- NodeToIndex(data, Ynodes)
  
  nodes <- list(A=Anodes, C=Cnodes, L=Lnodes, Y=Ynodes, AC=sort(c(Anodes, Cnodes)))
  nodes$LY <- CreateLYNodes(data, nodes, check.Qform=FALSE)
  return(nodes)
}

# Get the LY nodes but don't include "blocks" of L/Y nodes uninterrupted by A/C nodes
CreateLYNodes <- function(data, nodes, check.Qform, Qform) {
  LYnodes <- sort(c(nodes$L, nodes$Y))
  #if there are no A/C nodes between two or more LY nodes, only the first LY node in the block is considered an LY node
  nodes.to.remove <- NULL
  if (length(LYnodes) > 1) {
    for (i in 1:(length(LYnodes) - 1)) {
      cur.node <- LYnodes[i]
      next.node <- LYnodes[i + 1]
      if (! any(cur.node:next.node %in% nodes$AC)) {
        nodes.to.remove <- c(nodes.to.remove, next.node)
      }
    }
  }
  new.LYnodes <- setdiff(LYnodes, nodes.to.remove)
  if (check.Qform) {
    removed.Qform.index <- NULL
    for (i in nodes.to.remove) {
      index <- which(names(Qform) == names(data)[i])
      if (length(index) > 0) {
        removed.Qform.index <- c(removed.Qform.index, index)
      }
    }
    if (! is.null(removed.Qform.index)) {
      message("L/Y nodes (after removing blocks)  : ", names(data)[new.LYnodes], "\n")
      message("Qform names                        : ", names(Qform), "\n")
      message(paste("The following nodes are not being considered as L/Y nodes because they are part of a block of L/Y nodes. They are being dropped from Qform:\n"), paste(names(Qform)[removed.Qform.index], "\n", collapse=" "))
      Qform <- Qform[-removed.Qform.index]
    }
    return(list(LYnodes=new.LYnodes, Qform=Qform))
  }
  return(new.LYnodes)
}

# SL.library can be a character vector of library or a list with two separate vectors, one for Q and one for g
GetLibrary <- function(SL.library, estimate.type) {
  if (is.null(names(SL.library))) return(SL.library)
  if (! identical(sort(names(SL.library)), sort(c("Q", "g")))) stop("If SL.library has names, it must have two names: Q and g")
  if (! estimate.type %in% c("Q", "g")) stop("bad estimate.type")
  return(SL.library[[estimate.type]])
}

GetMsmWeights <- function(inputs) {
  n <- nrow(inputs$data)
  num.regimes <- dim(inputs$regimes)[3]
  stopifnot(num.regimes >= 1)
  num.final.Ynodes <- length(inputs$final.Ynodes)
  if (is.equal(inputs$msm.weights, "empirical")) {
    #default is probability of following abar given alive, uncensored; conditioning on past treatment/no censoring, but not L, W; duplicates get weight 0
    msm.weights <- matrix(nrow=num.regimes, ncol=num.final.Ynodes)
    
    for (j in 1:num.final.Ynodes) {
      final.Ynode <- inputs$final.Ynodes[j]
      inputs.subset <- SubsetInputs(inputs, final.Ynode)
      uncensored <- IsUncensored(inputs.subset$data, inputs.subset$nodes$C, cur.node=final.Ynode)
      if (dim(inputs.subset$regimes)[2] > 0) {
        is.duplicate <- duplicated(inputs.subset$regimes, MARGIN=3)
      } else {
        is.duplicate <- c(FALSE, rep(TRUE, num.regimes - 1))  #in case there are C nodes but no A nodes before a Ynode
      }
      for (i in 1:num.regimes) {
        if (is.duplicate[i]) {
          msm.weights[i, j] <- 0
        } else {
          intervention.match <- InterventionMatch(inputs.subset$data, abar=GetABar(inputs.subset$regimes, i), inputs.subset$nodes$A, cur.node=final.Ynode)
          msm.weights[i, j] <- sum(uncensored & intervention.match) / nrow(inputs.subset$data)
        } 
      }
    }
  } else if (is.null(inputs$msm.weights)) {
    msm.weights <- array(1, dim=c(n, num.regimes, num.final.Ynodes))
  } else {
    msm.weights <- inputs$msm.weights
  }
  if (is.equal(dim(msm.weights), c(num.regimes, num.final.Ynodes))) {
    msm.weights <- array(rep(msm.weights, each=n), dim=c(n, num.regimes, num.final.Ynodes))
  } else if (! is.equal(dim(msm.weights), c(n, num.regimes, num.final.Ynodes))) {
    stop("dim(msm.weights) should be c(n, num.regimes, num.final.Ynodes) or c(num.regimes, num.final.Ynodes)")
  }
  if (any(is.na(msm.weights)) || any(msm.weights < 0)) stop("all msm.weights must be >= 0 and not NA")
  return(msm.weights)
}

# Converts a general formula to a main terms formula and combine summary measures with baseline covariates
# Ex: If working.msm is "Y ~ X1*X2", convert to "Y ~ -1 + S1 + S1 + S3 + S4" where 
# S1 is 1 (intercept), S2 is X1, S3 is X2, S4 is X1:X2
ConvertToMainTerms <- function(data, msm, summary.measures, nodes) {
  baseline.column.names <- names(data)[seq(1, min(c(nodes$A, nodes$L, nodes$C, nodes$Y)) - 1)]
  summary.column.names <- colnames(summary.measures)
  rhs.vars <- RhsVars(msm)
  if (length(intersect(baseline.column.names, summary.column.names)) > 0) stop("Baseline covariate columns of data and columns of summary.measures may not have the same name")
  if (!all(rhs.vars %in% c(baseline.column.names, summary.column.names))) stop("All right hand side variables in working.msm must be either column names of summary.measures or column names of baseline covariates")
  baseline.column.names <- intersect(baseline.column.names, rhs.vars)
  baseline.data <- data[, baseline.column.names, drop=FALSE]
  num.regimes <- dim(summary.measures)[1]
  num.summary.measures <- dim(summary.measures)[2]
  num.final.Ynodes <- dim(summary.measures)[3]
  n <- nrow(data)
  for (j in 1:num.final.Ynodes) {
    for (i in 1:num.regimes) {
      combined.summary.measures <- model.matrix(as.formula(msm), data.frame(Y=1, baseline.data, matrix(summary.measures[i, , j], nrow=n, ncol=num.summary.measures, byrow=TRUE, dimnames=list(NULL, colnames(summary.measures)))))
      if (i == 1 && j == 1) {
        #initialize here now that we know how many columns there are
        main.terms.summary.measures <- array(dim=c(n, ncol(combined.summary.measures), num.regimes, num.final.Ynodes))
        beta.names <- colnames(combined.summary.measures) #this is the same for all i and j
      }
      main.terms.summary.measures[, , i, j] <- combined.summary.measures
    }
  }
  colnames(main.terms.summary.measures) <- paste("S", 1:ncol(main.terms.summary.measures), sep="") #temp names
  main.terms.msm <- paste("Y ~ -1 +", paste(colnames(main.terms.summary.measures), collapse=" + ")) #formula using temp names 
  return(list(msm=main.terms.msm, summary.measures=main.terms.summary.measures, beta.names=beta.names, baseline.column.names=baseline.column.names))
}

# Convert censoring nodes stored as binaries into factors (factors are recommended but binaries are currently accepted)
ConvertCensoringNodes <- function(data, Cnodes, has.deterministic.functions=FALSE) {
  error.msg <- "in data, all Cnodes should be factors with two levels, 'censored' and 'uncensored'\n See ?BinaryToCensoring \n (binary is also accepted, where 0=censored, 1=uncensored, but this is not recommended)"
  for (i in Cnodes) {
    col <- data[, i]
    if (is.numeric(col)) {
      if (! all(col %in% c(0, 1, NA))) stop(error.msg)
      data[, i] <- BinaryToCensoring(is.uncensored=col)
      if (has.deterministic.functions) warning("Censoring nodes have been converted from binaries to factors - see ?BinaryToCensoring.\n Note that if you are writing your own deterministic.g.function or deterministic.Q.function that censoring nodes are converted to factors\n before these functions are called.")
    } else if (is.factor(col)) {
      if (! all(levels(col) %in% c("censored", "uncensored"))) {
        stop("all levels of data[, Cnodes] should be in censored, uncensored (NA should not be a level)")
      }
      #no action required
    } else {
      stop(error.msg)
    } 
  }
  return(data)
}

#Before passing data to SuperLearner, convert factors to binary
ConvertCensoringNodesToBinary <- function(data, Cnodes) {
  CensoringToBinary <- function(x) {
    if (! all(levels(x) %in% c("censored", "uncensored"))) {
      stop("all levels of data[, Cnodes] should be in censored, uncensored (NA should not be a level)")
    }
    b <- rep(NA_integer_, length(x))
    b[x == "censored"] <- 0L
    b[x == "uncensored"] <- 1L
    return(b)
  }
  
  error.msg <- "in data, all Cnodes should be factors with two levels, 'censored' and 'uncensored' \n (binary is also accepted, where 0=censored, 1=uncensored, but is not recommended)"
  for (i in Cnodes) {
    col <- data[, i]
    if (is.numeric(col)) {
      if (! all(col %in% c(0, 1, NA))) stop(error.msg)
    } else if (is.factor(col)) {
      data[, i] <- CensoringToBinary(col)
    } else {
      stop(error.msg)
    } 
  }
  return(data)
}

# We don't want to show all of the warnings 
GetWarningsToSuppress <- function(update.step=FALSE) {
  warnings.to.suppress <- c("glm.fit: fitted probabilities numerically 0 or 1 occurred", "prediction from a rank-deficient fit may be misleading") 
  if (update.step) {
    # It's ok if the updating step doesn't converge, we'll fix it in FixScoreEquation 
    warnings.to.suppress <- c(warnings.to.suppress, "glm.fit: algorithm did not converge")
  }
  return(warnings.to.suppress)
}

SetSeedIfRegressionTesting <- function() {
  #if comparing outputs between different versions of ltmle, we need to sync random numbers 
  #before calling SuperLearner or FixScoreEquation since these use random numbers
  seed <- Sys.getenv("LTMLE.REGRESSION.TESTING.SEED") 
  stopifnot(length(seed) == 1)
  if (seed != "") {
    seed <- as.numeric(seed)
    stopifnot(is.finite(seed))
    set.seed(seed)
    cat("set seed 0.9-7\n") 
  }
  invisible(NULL)
}

Default.SL.Library <- list("SL.glm",
                           "SL.stepAIC",
                           "SL.bayesglm", 
                           c("SL.glm", "screen.corP"), 
                           c("SL.step", "screen.corP"), 
                           c("SL.step.forward", "screen.corP"), 
                           c("SL.stepAIC", "screen.corP"), 
                           c("SL.step.interaction", "screen.corP"), 
                           c("SL.bayesglm", "screen.corP")
)  

