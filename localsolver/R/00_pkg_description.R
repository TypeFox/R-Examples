#' R API to LocalSolver
#'
#' The package converts R data onto input and data for LocalSolver, executes optimization and exposes 
#' optimization results as R data.
#' 
#' \emph{LocalSolver} (\url{http://www.localsolver.com/}) is an optimization engine developed 
#' by \emph{Innovation24} (\url{http://www.innovation24.fr/}). It is designed to solve large-scale 
#' mixed-variable non-convex optimization problems.
#' 
#' The localsolver package is result of cooperation of \emph{WLOG Solutions} (\url{http://www.wlogsolutions.com/en/}) 
#' in collaboration Decision Support and Analysis Division at \emph{Warsaw School of Economics} 
#' (\url{http://www.sgh.waw.pl/en/}).
#' 
#' The localsolver package allows for solving optimization problems in R by using the LocalSolver program. 
#' It combines efficiency and easiness of model formulation, which characterize the LocalSolver program, 
#' with the elasticity of data in R. The model formulation in LSP language is provided by a user as a string
#' or text file, and the input data - as a list. The solution is also a list, with chosen decision variables,
#' constraints and objective function values. Thus the problem can be solved many times (for example, by 
#' iterating with different parameter settings) and the result is provided in a form which makes further 
#' processing or visualization in R fast and easy.
#' 
#' It is highly recommended that the LocalSolver software is installed before the R package.\cr
#' The test version (limited to 100 decisions and 1000 variables) can be downloaded from the website:
#'   \url{http://www.localsolver.com/download.html}\cr
#' The access to the trial license for full version requires contacting the \emph{Innovation24} company. For further
#' information please visit the website: \url{http://www.localsolver.com/}
#'   
#' Bugs and feature requests can be reported by e-mail: \emph{Walerian Sokolowski <walerian.sokolowski@@wlogsolutions.com>}. We are also 
#' looking forward for any future package development suggestions.
#'  
#' @docType package
#' @name localsolver
NULL


#' Data set for google_machine_reassignment demo.
#' 
#' @docType data
#' @keywords datasets
#' @name google_example_1_data
NULL
