##' Sparsity by Worst-Case Quadratic Penalties
##'
##' This package is designed to fit accurately several popular
##' penalized linear regression models using the algorithm proposed in
##' Grandvalet, Chiquet and Ambroise (submitted) by solving quadratic
##' problems with increasing size.
##'
##' @section Features:
##'
##' At the moment, two \code{R} fitting functions are available:
##' \enumerate{%
##'
##' \item the \code{\link{elastic.net}} function, which solves a family of
##' linear regression problems penalized by a mixture of
##' \eqn{\ell_1}{l1} and \eqn{\ell_2}{l2} norms. It notably includes
##' the LASSO (Tibshirani, 1996), the adaptive-LASSO (Zou, 2006), the
##' Elastic-net (Zou and Hastie, 2006) or the Structured Elastic-net
##' (Slawski et al., 2010). See examples as well as the available
##' \code{demo(quad_enet)}.
##'
##' \item the \code{\link{bounded.reg}} function, which fits a linear model
##' penalized by a mixture of \eqn{\ell_\infty}{infinity} and
##' \eqn{\ell_2}{l2} norms. It owns the same versatility as the
##' \code{elastic.net} function regarding the \eqn{\ell_2}{l2} norm,
##' yet the \eqn{\ell_1}{l1}-norm is replaced by the infinity
##' norm. Check \code{demo(quad_breg)} and examples.}
##'
##' The problem commonly solved for these two functions writes
##' \if{latex}{\deqn{% \hat{\beta}_{\lambda_1,\lambda_2} = \arg
##' \min_{\beta} \frac{1}{2} (y - X \beta)^T (y - X \beta) + \lambda_1
##' \|D \beta \|_{q} + \frac{\lambda_2}{2} \beta^T S \beta, }}
##' \if{html}{\out{ <center> &beta;<sup>hat</sup>
##' <sub>&lambda;<sub>1</sub>,&lambda;<sub>2</sub></sub> =
##' argmin<sub>&beta;</sub> 1/2 RSS(&beta) + &lambda;<sub>1</sub>
##' &#124; D &beta; &#124;<sub>q</sub> + &lambda;/2 <sub>2</sub>
##' &beta;<sup>T</sup> S &beta;, </center> }}
##' \if{text}{\deqn{beta.hat(lambda1, lambda2) = argmin_beta 1/2
##' RSS(beta) + lambda1 |D beta|_q + lambda2 beta' S beta,}} where
##' \eqn{q=1}{q=1} for \code{elastic.net} and
##' \eqn{q=\infty}{q=infinity} for \code{bounded.reg}.  The diagonal
##' matrix \eqn{D}{D} allows different weights for the first part of
##' the penalty. The structuring matrix \eqn{S}{S} can be used to
##' introduce some prior information regarding the predictors. It is
##' provided via a positive semidefinite matrix.
##'
##' The S4 objects produced by the fitting procedures own the
##' classical methods for linear model in \code{R}, as well as methods
##' for plotting, (double) cross-validation and for the stability
##' selection procedure of Meinshausen and Buhlmann (2010).
##'
##' All the examples of this documentation have been included to the
##' package source, in the 'examples' directory. Some (too few!)
##' routine testing scripts using the \pkg{testhat} package are also
##' present in the 'tests' directory, where we check basic
##' functionalities of the code, especially the reproducibility of the
##' Lasso/Elastic-net solution path with the \pkg{lars},
##' \pkg{elasticnet} and \pkg{glmnet} packages.  We also check the
##' handling of runtime errors or unstabilities.
##'
##' @section Algorithm:
##'
##' The general strategy of the algorithm relies on maintaining an
##' active set of variables, starting from a vector of zeros. The
##' underlying optimization problem is solved only on the activated
##' variables, thus handling with small smooth problems with
##' increasing size. Hence, by considering a decreasing grid of values
##' for the penalty \eqn{\lambda_1}{lambda1} and fixing
##' \eqn{\lambda_2}{lambda2}, we may explore the whole path of
##' solutions at a reasonable numerical cost, providing that
##' \eqn{\lambda_1}{lambda1} does not end up too small.
##'
##' For the \eqn{\ell_1}{l1}-based methods (available in the
##' \code{elastic.net} function), the size of the underlying problems
##' solved is related to the number of nonzero coefficients in the
##' vector of parameters. With the \eqn{\ell_\infty}{infinity}-norm,
##' (available in the \code{boundary.reg} function), we do not produce
##' sparse estimator. Nevertheless, the size of the systems solved
##' along the path deals with the number of unbounded variables for
##' the current penalty level, which is quite smaller than the number
##' of predictors for a reasonable \eqn{\lambda_1}{lambda1}. The same
##' kind of proposal was made in Zhao, Rocha and Yu (2009).
##'
##' Underlying optimization is performed by direct resolution of
##' quadratic sub problems, which is the main purpose of this
##' package. This strategy is thoroughly exposed in Grandvalet,
##' Chiquet and Ambroise (submitted). Still, we also implemented the
##' popular and versatile proximal (FISTA) approaches for routine
##' checks and numerical comparisons. A coordinate descent approach is
##' also included, yet only for the \code{elastic.net} fitting
##' procedure.
##'
##' The default setting uses the quadratic approach that gives its
##' name to the package. It has been optimized to be the method of
##' choice for small and medium scale problems, and produce very
##' accurate solutions. However, the first order methods (coordinate
##' descent and FISTA) can be interesting in situations where the
##' problem is close to singular, in which case the Cholesky
##' decomposition used in the quadratic solver can be computationally
##' unstable. Though it is extremely unlikely for
##' \code{\link{elastic.net}} -- and if so, we encourage the user to
##' send us back any report of such an event --, this happens at times
##' with \code{\link{bounded.reg}}. Regarding this issue, we let the
##' possibility for the user to run the optimization of the
##' \code{\link{bounded.reg}} criterion in a (hopefully) 'bulletproof'
##' mode: using mainly the fast and accurate quadratic approach, it
##' switches to the slower but more robust proximal resolution when
##' unstability is detected.
##'
##' @section Technical remarks:
##'
##' Most of the numerical work is done in C++, relying on the
##' \pkg{RcppArmadillo} package. We also provide a (double)
##' cross-validation procedure and functions for stabilty selection,
##' both using the multi-core capability of the computer, through the
##' \pkg{parallel} package. This feature is not available for Windows
##' user, though. Finally, note that the plot methods enjoy some
##' (still very few) of the capabilities of the \pkg{ggplot2} package.
##'
##' We hope to enrich \pkg{quadrupen} with other popular fitting
##' procedures and develop other statistical tools, particularly
##' towards bootstrapping and model selection purpose. Sparse matrix
##' encoding is partially supported at the moment, and will hopefully
##' be thoroughly available in the future, thanks to upcoming updates
##' of the great \pkg{RcppArmadillo} package.
##'
##' @name quadrupen-package
##' @aliases quadrupen
##' @docType package
##' @author Julien Chiquet \email{julien.chiquet@@genopole.cnrs.com}
##'
##' @references
##' Yves Grandvalet, Julien Chiquet and Christophe Ambroise,
##' \href{http://arxiv.org/abs/1210.2077}{Sparsity by Worst-case Quadratic
##' Penalties}, arXiv preprint, 2012.
##'
##' \itemize{%
##'
##' \item Nicolas Meinshausen and Peter Buhlmann. Stability Selection,
##' JRSS(B), 2010.
##' \item  Martin Slawski, Wolfgang zu Castell, and Gerhard
##' Tutz. Feature selection guided by structural information, AOAS,
##' 2010.
##' \item Peng Zhao, Guillerme Rocha and Bin Yu. The composite
##' absolute penalties family for grouped and hierarchical variable
##' selection, The Annals of Statistics, 2009.
##' \item  Hui Zou. The Adaptive Lasso and Its Oracle Properties,
##' JASA, 2006.
##' \item Hui Zou and Trevor Hastie. Regularization and variable
##' selection via the elastic net, JRSS(B), 2006.
##' \item  Robert Tibshirani. Regression Shrinkage and Selection
##' via the Lasso, JRSS(B), 1996.
##' }
##' @import Matrix parallel Rcpp methods ggplot2 reshape2 scales grid
##' @useDynLib quadrupen
NULL
