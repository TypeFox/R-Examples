##' Distributed gradient boosting based on the \pkg{mboost} package.
##'
##' The parboost package implements distributed gradient boosting based on
##' the \pkg{mboost} package. When should you use parboost instead of mboost?
##' There are two use cases:
##' 1. The data takes too long to fit as a whole
##' 2. You want to bag and postprocess your boosting models to get a
##'     more robust ensemble
##' parboost is designed to scale up component-wise functional gradient
##' boosting in a distributed memory environment by splitting the
##' observations into disjoint subsets. Alternatively, parboost can
##' generate and use bootstrap samples of the original data. Each cluster
##' node then fits a boosting model to its subset of the data. These
##' boosting models are combined in an ensemble, either with equal
##' weights, or by fitting a (penalized) regression model on the
##' predictions of the individual models on the complete data. All other
##' functionality of mboost is left untouched for the moment.
##'
##' @name parboost
##' @docType package
##' @import parallel iterators
##' @importFrom doParallel registerDoParallel
##' @importFrom glmnet cv.glmnet
##' @importFrom mboost mboost blackboost selected mstop boost_control
##' Gaussian Binomial Poisson
##' @importFrom plyr laply
##' @importFrom party ctree_control
##' @importFrom caret createResample createFolds
##' @references Peter Buehlmann and Bin Yu (2003),
##' Boosting with the L2 loss: regression and classification.
##' \emph{Journal of the American Statistical Association}, \bold{98},
##' 324--339.
##'
##' Peter Buehlmann and Torsten Hothorn (2007),
##' Boosting algorithms: regularization, prediction and model fitting.
##' \emph{Statistical Science}, \bold{22}(4), 477--505.
##'
##' Torsten Hothorn, Peter Buehlmann, Thomas Kneib, Mattthias Schmid and
##' Benjamin Hofner (2010), Model-based Boosting 2.0. \emph{Journal of
##' Machine Learning Research}, \bold{11}, 2109--2113.
##'
##' Yoav Freund and Robert E. Schapire (1996),
##' Experiments with a new boosting algorithm.
##' In \emph{Machine Learning: Proc. Thirteenth International Conference},
##' 148--156.
##'
##' Jerome H. Friedman (2001),
##' Greedy function approximation: A gradient boosting machine.
##' \emph{The Annals of Statistics}, \bold{29}, 1189--1232.
##'
##' Benjamin Hofner, Andreas Mayr, Nikolay Robinzonov and Matthias Schmid
##' (2012). Model-based Boosting in R: A Hands-on Tutorial Using the R
##' Package mboost. \emph{Department of Statistics, Technical Report No. 120}.\cr
##' \url{http://epub.ub.uni-muenchen.de/12754/}
##'
##' T. Hothorn, P. Buehlmann, T. Kneib, M. Schmid, and B. Hofner
##' (2013). mboost: Model-Based Boosting, R package version 2.2-3,
##' \url{http://CRAN.R-project.org/package=mboost}.
NULL

##' Benchmark Problem Friedman 2
##'
##' Dataset taken from \pkg{mlbench}. The inputs are 4 independent variables uniformly
##' distributed over the ranges
##' \deqn{0 \le x1 \le 100}
##' \deqn{40 \pi \le x2 \le 560 \pi}
##' \deqn{0 \le x3 \le 1}
##' \deqn{1 \le x4 \le 11}
##' The outputs are created according to the formula
##' \deqn{y = (x1^2 + (x2 x3 - (1/(x2 x4)))^2)^{0.5} + e}
##' where \emph{e} is \eqn{N(0,sd)}.
##'
##' @docType data
##' @keywords datasets
##' @format A data frame with 100 rows and 5 variables
##' @name friedman2
##' @source \url{http://cran.r-project.org/web/packages/mlbench/index.html}
##' @references Breiman, Leo (1996) Bagging predictors. Machine Learning 24, pages
##' 123-140.
##'
##' Friedman, Jerome H. (1991) Multivariate adaptive regression
##' splines. The Annals of Statistics 19 (1), pages 1-67.
##'
##' Friedrich Leisch & Evgenia Dimitriadou (2010). mlbench: Machine
##'  Learning Benchmark Problems. R package version 2.1-1.
NULL
