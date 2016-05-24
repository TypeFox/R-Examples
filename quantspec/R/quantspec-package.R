#' @include quantspec-package.R
NULL

################################################################################
#' Quantile-Based Spectral Analysis of Time Series
#'
#' Methods to determine, smooth and plot quantile
#' periodograms for univariate and (since v1.2-0) multivariate time series.
#'
#' @details
#'  \tabular{ll}{
#'    \cr Package: \tab quantspec
#'    \cr Type:    \tab Package
#'    \cr Version: \tab 1.2-1
#'    \cr Date:    \tab 2016-03-27
#'    \cr License: \tab GPL (>= 2)
#'  }
#'
#' @section Contents:
#' The \pkg{quantspec} package contains a hierachy of S4 classes with
#' corresponding methods and functions serving as constructors. The following
#' class diagrams provide an overview on the structure of the package. In the
#' first and second class diagram the classes implementing the estimators are
#' shown. In the first diagram the classes related to periodogram-based
#' estimation are displayed:  
#'
#' \if{html}{\figure{main-mv.png}{options: width=960}}
#' \if{latex}{\figure{main-mv.pdf}{options: width=12cm}}
#' 
#' In the second diagram the classes related to lag window-based
#' estimation are displayed:
#' 
#' \if{html}{\figure{main2-mv.png}{options: width=768}}
#' \if{latex}{\figure{main2-mv.pdf}{options: width=8cm}}
#'
#' In the third class diagram the classes implementing model quantities are
#' displayed. A relation to the ``empirical classes'' is given via the fact that
#' the quantile spectral densities are computed by simulation of quantile
#' periodograms and a common abstract superclass \code{QSpecQuantity} which
#' is used to provide a common interface to quantile spectral quantities.
#'
#' \if{html}{\figure{csd-mv.png}{options: width=768}}
#' \if{latex}{\figure{csd-mv.pdf}{options: width=12cm}}
#'
#' Besides the object-oriented design a few
#' auxiliary functions exists. They serve as parameters or are mostly for
#' internal use. A more detailed description of the framework can be found in
#' the paper on the package (Kley, 2016).
#'
#' @section Organization of the source code / files in the \code{/R} folder:
#' All of the source code related to the specification of a certain class is
#' contained in a file named \code{Class-[Name_of_the_class].R}. This includes,
#' in the following order,
#' \enumerate{
#'   \item all roxygen \code{@@include} to insure the correctly generated
#'         collate for the DESCRIPTION file.
#'   \item \code{\\setClass} preceded by a meaningful roxygen documentation.
#'   \item specification of an \code{initialize} method, where appropriate.
#'   \item all accessor and mutator method (i. e., getter and setter); first
#'         the ones returning attributes of the object, then the ones returning
#'         associated objects.
#'   \item constructors; use generics if there is more than one of them.
#'   \item \code{show} and \code{plot} methods.
#' }
#'
#' @section Coding Conventions:
#' To improve readability of the software and documentation this package was
#' written in the spirit of the ``Coding conventions of the Java Programming
#' Language'' (Oracle, 2015). In particular, the naming conventions for classes
#' and methods have been adopted, where ``Class names should be nouns, in mixed
#' case with the first letter of each internal word capitalized.'' and
#' ``Methods should be verbs, in mixed case with the first letter lowercase,
#' with the first letter of each internal word capitalized.''
#'
#' @section Naming Conventions for the Documentation:
#' To reflect the structure of the contents of the package in the documentation
#' file, the following system for naming of the sections is adopted:
#' \itemize{
#'   \item Documentation of an S4 class is named as the name of the class
#'         followed by ``-class''. [cf. \code{\link{QuantilePG-class}}]
#'  \item Documentation of a constructor for an S4-class is named as
#'         the name of the class followed by ``-constructor''.
#'         [cf. \code{\link{QuantilePG-constructor}}]
#'  \item Documentation of a method dispaching to an object of a certain
#'         S4 class is named by the name of the method, followed by ``-'',
#'         followed by the name of the Class.
#'         [cf. \code{\link{getValues-QuantilePG}}]
#' }
#'
#' @name quantspec-package
#' @aliases quantspec
#' @docType package
#' @author Tobias Kley
#'
#' @import graphics
#' @import methods
#' @import stats4
#' @importFrom Rcpp evalCpp
#' 
#' @useDynLib quantspec

#'
#' @references
#' Kley, T. (2014a). Quantile-Based Spectral Analysis: Asymptotic Theory and
#' Computation. Ph.D. Dissertation, Ruhr University Bochum.
#' \url{http://www-brs.ub.ruhr-uni-bochum.de/netahtml/HSS/Diss/KleyTobias/}.
#'
#' Kley, T. (2016). Quantile-Based Spectral Analysis in an Object-Oriented
#' Framework and a Reference Implementation in R: The quantspec Package.
#' Journal of Statistical Software, \bold{70}(3), 1--27.
#'
#' Dette, H., Hallin, M., Kley, T. & Volgushev, S. (2015).
#' Of Copulas, Quantiles, Ranks and Spectra: an \eqn{L_1}{L1}-approach to
#' spectral analysis. \emph{Bernoulli}, \bold{21}(2), 781--831.
#' [cf. \url{http://arxiv.org/abs/1111.7205}]
#'
#' Kley, T., Volgushev, S., Dette, H. & Hallin, M. (2016).
#' Quantile Spectral Processes: Asymptotic Analysis and Inference.
#' \emph{Bernoulli}, \bold{22}(3), 1770--1807.
#' [cf. \url{http://arxiv.org/abs/1401.8104}]
#' 
#' Barunik, J. & Kley, T. (2015).
#' Quantile Cross-Spectral Measures of Dependence between Economic Variables.
#' [cf. \url{http://arxiv.org/abs/1510.06946}]
#'
#' Oracle (2015). Coding conventions of the Java Programming Language.
#' \url{http://www.oracle.com/technetwork/java/codeconvtoc-136057.html}.
#' Accessed 2015-03-25.
#'
NULL

# Taken from quantreg-package and adapted.
".onAttach" <- function(lib, pkg) {
  if(interactive() || getOption("verbose"))
    packageStartupMessage("Package quantspec loaded.\n     To cite, see citation(\"quantspec\").\n     For demos, see demo(package = \"quantspec\").")
}
