#' SolverOptions
#'
#' Documentation for the SolverOptions parameter
#'
#' SolverOptions (optional) is a list containing: 'DDEsolver', 'tol', 'hbsize', 'method' and 'atol'. DDEsolver equal to 'deSolve' or 'PBS' sets the R package used to solve the DDEs. The 'tol' option sets the relative tolerances and 'hbsize' sets the size of the history buffer. The remaining two items, 'method' and 'atol' set the numerical integration scheme and the absolute tolerance if 'DDEsolver'='deSolve' (PBS does not have these options). If solverOptions is not specified at all, or if only some of the options are specified, the default values: list('DDEsolver'='deSolve', 'tol'=1e-7, 'hbsize'=1e3, 'method'=lsoda, 'atol'=1e-7, 'dt'=0.1) will be used.
#' @aliases SolverOptions
#' @name SolverOptions
NULL
