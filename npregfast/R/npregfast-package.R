#' \code{npregfast}: Nonparametric Estimation of 
#' Regression Models with Factor-by-Curve Interactions. 
#'
#'
#' This package provides a method for obtain nonparametric estimates of regression
#' models using local polynomial kernel smoothers or splines. Particular 
#' features of the package
#' are facilities for fast smoothness estimation, and the calculation of their 
#' first and second derivative. Users can define the smoothers parameters. 
#' Confidence intervals calculation is provided by bootstrap methods. 
#' Binning techniques were applied to speed up computation in the estimation 
#' and testing processes.
#'
#' @name npregfast
#' @docType package
#' @details \tabular{ll}{ Package: \tab npregfast\cr Type: \tab Package\cr
#' Version: \tab 1.3.0\cr Date: \tab 2016-04-29\cr License: \tab MIT + file LICENSE\cr}
#'
#'\code{npregfast} is designed along lines similar to those of other \code{R} 
#'regression packages. The main function of the library is \code{frfast} 
#'which, by default, fits a nonparametric regression model based on local 
#'polynomial kernel smoothers. Note that through the argument \code{formula}
#' users can decide to fit a model by taking or not taking the interaction 
#' into account and by the argument \code{formula} it is posible to select 
#' the type of smoother: kernel or splines. Numerical and graphical summaries 
#' of the fitted object can be 
#' obtained by using the generic functions, \code{print.frfast}, 
#' \code{summary.frfast} and \code{plot.frfast}. Another of these generic 
#' functions is \code{predict.frfast}, which takes a fitted model of the 
#' \code{frfast} class and, given a new data set of values of the covariate, 
#' produces predictions.
#' As mentioned above, this package can be used to fit models taking into 
#' account factor-by-curve interactions. In this framework, it will be 
#' necessary to ascertain if the factor produces an effect on the response 
#' and thus, there is a interaction or, in contrast, the estimated regression 
#' curves are equal. To this end, the package provides the \code{globaltest} 
#' function which answers this question through a bootstrap-based test. 
#' If the factor results significant, then \code{plotdiff()} enables the user 
#' to obtain a graphical representation that shows the differences between 
#' the estimated curves (estimate, first or second derivative) for any set of 
#' two levels of the factor. Additionally, with \code{critical()} it is possible
#'  to obtain the value of the covariate that maximises the estimate and 
#'  first derivative of the function and the value of the covariate that equals 
#'  the second derivative to zero, for each of these levels. Again, to test if 
#'  these estimated points are equal for all levels, the package provides the 
#'  \code{localtest} function. Note that, to compare these points between 
#'  any set of two levels, a confidence interval for the difference can be 
#'  obtained by applying \code{criticaldiff()}.
#'
#'
#'
#'
#' For a listing of all routines in the NPRegfast package type:
#' \code{library(help="npregfast")}. 
#' 
#' View a \href{http://sestelo.shinyapps.io/npregfast}{demo Shiny app}
#' or see the full \href{https://github.com/sestelo/npregfast}{README} on GitHub.
#' 
#' 
#' @author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#' 
#' @references 
#' Efron, B. (1979). Bootstrap methods: another look at the jackknife. 
#' Annals of Statistics, 7, 1--26.
#' 
#' Efron, E. and Tibshirani, R. J. (1993). An introduction to the Bootstrap. 
#' Chapman and Hall, London.
#' 
#' Huxley, J. S. (1924). Constant differential growth-ratios and their 
#' significance. Nature, 114:895--896.
#' 
#' Sestelo, M. (2013). Development and computational implementation of 
#' estimation and inference methods in flexible regression models. 
#' Applications in Biology, Engineering and Environment. PhD Thesis, Department
#' of Statistics and O.R. University of Vigo.
#' 
#' Sestelo, M. and Roca-Pardinas, J. (2011). A new approach to estimation of 
#' length-weight relationship of \eqn{Pollicipes}  \eqn{pollicipes} 
#' (Gmelin, 1789) on the Atlantic coast of Galicia (Northwest Spain): some 
#' aspects of its biology and management. Journal of Shellfish Research, 
#' 30(3), 939--948.
#' 
#' Wand, M. P. and Jones, M. C. (1995). Kernel Smoothing. Chapman & Hall, London.
#'
#'
#'
#'



NULL