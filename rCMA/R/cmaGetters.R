#
# Getters and Setters
#
######################################################################################
#' rCMA Getters and Setters.
#'
#' Get or set various elements of CMA-ES Java object \code{cma}. \cr
#' \cr
#' \code{cmaSetDimension} sets the problem dimension (only prior to \code{\link{cmaInit}}) \cr
#' \code{cmaGetDimension} returns the problem dimension \cr
#' \code{cmaSetPopulationSize} sets the population size (only prior to \code{\link{cmaInit}}) \cr
#' \code{cmaGetPopulationSize} returns the population size  \cr
#' \code{cmaSetInitialX} set the mean vector for the initial population (only prior to \code{\link{cmaInit}}) \cr
#' \code{cmaGetInitialX} returns the mean vector for the initial population  \cr
#' \code{cmaSetCountEval} sets the counter for fitness function evaluations (only prior to \code{\link{cmaInit}}) \cr
#' \code{cmaGetCountEval} returns the counter for fitness function evaluations  \cr
#'
#' @param cma CMA-ES Java object, created with \code{\link{cmaNew}}
#' @param i   a parameter of type integer
#' @param p   a parameter of type long
#' @param initialX  either a double or a double vector of length  \code{\link{cmaGetDimension}}
#'
#' @return none for the setters, the requested element(s) for the getters
#' @seealso   \code{\link{cmaSetStopFitness}}, \code{\link{cmaNew}}, \code{\link{cmaInit}}
#' @export
cmaSetDimension <- function(cma,i) {
  rJava::.jcall(cma,,"setDimension",as.integer(i));
}
#' @rdname cmaSetDimension
#' @export
cmaGetDimension <- function(cma) {
  rJava::.jcall(cma,"I","getDimension");
}
#' @rdname cmaSetDimension
#' @export
cmaSetPopulationSize <- function(cma,i) {
  rJava::.jcall(cma,,"setPopulationSize",as.integer(i));
}
#' @rdname cmaSetDimension
#' @export
cmaGetPopulationSize <- function(cma) {
  rJava::.jcall(cma,"I","getPopulationSize");
}
#' @rdname cmaSetDimension
#' @export
cmaSetInitialX <- function(cma,initialX) {
  rJava::.jcall(cma,,"setInitialX",initialX);
}
#' @rdname cmaSetDimension
#' @export
cmaGetInitialX <- function(cma) {
  rJava::.jcall(cma,"[D","getInitialX");
}
#' @rdname cmaSetDimension
#' @export
cmaSetCountEval <- function(cma,p) {
  rJava::.jcall(cma,"J","setCountEval",rJava::.jlong(p));
}
#' @rdname cmaSetDimension
#' @export
cmaGetCountEval <- function(cma) {
  rJava::.jcall(cma,"J","getCountEval");
}

######################################################################################
#' rCMA Stop Conditions.
#'
#' Set various stop conditions of CMA-ES Java object \code{cma} (only prior to \code{\link{cmaInit}}). \cr
#' \cr
#' \code{cmaSetStopFitness} sets the stop condition: fitness function below \code{d} (default: DOUBLE.MinValue)  \cr
#' \code{cmaSetStopMaxFunEvals} sets the stop condition: max number of fitness function evaluations  \cr
#' \code{cmaSetStopTolFun} sets the stop condition: delta of fitness function below \code{d} (default: 1e-12)  \cr
#'
#' @param cma CMA-ES Java object, created with \code{\link{cmaNew}}
#' @param p   a parameter of type long
#' @param d   a parameter of type double
#'
#' @note
#'    If your fitness can become negative, you need to set \code{cmaSetStopFitness} to a value different
#'    from the default to prevent premature stopping.
#'
#'    The properties file (read by \code{\link{cmaNew}}) can be used to set further stop conditions. 
#'    If they are not set, the following defaults are active:
#'    \tabular{ccc}{
#'      name            \tab  default setting \tab  meaning    \cr
#'      stopTolFunHist  \tab  1e-13           \tab  similar to stopTolFun, see CMA-ES Javadoc for details     \cr
#'      stopTolX        \tab  0.0             \tab  stop if seacrh steps become smaller than stopTolX     \cr
#'      stopTolXfactor  \tab  0.0             \tab  stop if search steps become smaller than stopTolXFactor * initial step size    \cr
#'      stopMaxIter     \tab  +Inf            \tab  stop if number of iterations (generations) are greater     
#'    }
#'  
#' @seealso   \code{\link{cmaSetDimension}}, \code{\link{cmaNew}}, \code{\link{cmaInit}}
#' @export
cmaSetStopFitness <- function(cma,d) {
  rJava::.jcall(cma,,"setOptionsStopFitness",as.double(d));
}
#' @rdname cmaSetStopFitness
#' @export
cmaSetStopMaxFunEvals <- function(cma,p) {
  rJava::.jcall(cma,,"setOptionsStopMaxFunEvals",rJava::.jlong(p));
}
#' @rdname cmaSetStopFitness
#' @export
cmaSetStopTolFun <- function(cma,d) {
  rJava::.jcall(cma,,"setOptionsStopTolFun",as.double(d));
}

