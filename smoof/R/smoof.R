#' smoof: Single and Multi-Objective Optimization test functions.
#'
#' The \pkg{smoof} R package provides generators for huge set of single- and
#' multi-objective test functions, which are frequently used in the literature
#' to benchmark optimization algorithms. Moreover the package provides methods
#' to create arbitrary objective functions in an object-orientated manner, extract
#' their parameters sets and visualize them graphically.
#'
#' @section Some more details:
#' Given a set of criteria \eqn{\mathcal{F} = \{f_1, \ldots, f_m\}} with each
#' \eqn{f_i : S \subseteq \mathbf{R}^d \to \mathbf{R} , i = 1, \ldots, m} being an
#' objective-function, the goal in \emph{Global Optimization (GO)} is to find the best
#' solution \eqn{\mathbf{x}^* \in S}. The set \eqn{S} is termed the \emph{set of
#' feasible soluations}. In the case of only a single objective function \eqn{f},
#' - which we want to restrict ourself in this brief description - the goal is to
#' minimize the objective, i. e., \deqn{\min_{\mathbf{x}} f(\mathbf{x}).}
#' Sometimes we may be interested in maximizing the objective function value, but
#' since \eqn{min(f(\mathbf{x})) = -\min(-f(\mathbf{x}))}, we do not have to tackle
#' this separately.
#' To compare the robustness of optimization algorithms and to investigate their behaviour
#' in different contexts, a common approach in the literature is to use \emph{artificial
#' benchmarking functions}, which are mostly deterministic, easy to evaluate and given
#' by a closed mathematical formula.
#' A recent survey by Jamil and Yang lists 175 single-objective benchmarking functions
#' in total for global optimization [1]. The \pkg{smoof} package offers implementations
#' of a subset of these functions beside some other functions as well as
#' generators for large benchmarking sets like the noiseless BBOB2009 function set [2]
#' or functions based on the multiple peaks model 2 [3].
#'
#' @references
#' [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
#' functions for global optimization problems, Int. Journal of Mathematical
#' Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150-194 (2013).
#' [2] Hansen, N., Finck, S., Ros, R. and Auger, A. Real-Parameter Black-Box
#' Optimization Benchmarking 2009: Noiseless Functions Definitions. Technical report
#' RR-6829. INRIA, 2009.
#' [3] Simon Wessing, The Multiple Peaks Model 2, Algorithm Engineering Report
#' TR15-2-001, TU Dortmund University, 2015.
#'
#' @docType package
#' @name smoof-package
NULL
