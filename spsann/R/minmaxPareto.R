#' Pareto minimum and maximum values
#' 
#' Compute the minimum and maximum attainable values of the objective functions that compose a multi-objective 
#' combinatorial optimization problem.
#' 
#' @inheritParams optimACDC
#' @inheritParams spJitter
#' 
#' @param osc A list of objects of class \code{OptimizedSampleConfiguration} (OSC). Each OSC of the list must
#' be named after the objective function with which it has been optimized. For example,
#' \code{osc = list(CORR = osc_corr, DIST = osc_dist)}.
#' 
#' @details 
#' \subsection{Multi-objective combinatorial optimization problems}{
#' A method of solving a multi-objective combinatorial optimization problem (MOCOP) is to aggregate the 
#' objective functions into a single \emph{utility function}. In \pkg{spsann}, the aggregation 
#' is performed using the \emph{weighted sum method}, which incorporates in the weights the preferences of 
#' the user regarding the relative importance of each objective function.
#' 
#' The weighted sum method is affected by the relative magnitude of the different function values. The
#' objective functions implemented in \pkg{spsann} have different units and orders of magnitude.
#' The consequence is that the objective function with the largest values may have a numerical dominance 
#' during the optimization. In other words, the weights may not express the true preferences of the user,
#' resulting that the meaning of the utility function becomes unclear because the optimization will favour the
#' objective function which is numerically dominant.
#' 
#' A reasonable solution to avoid the numerical dominance of any objective function is to scale the objective
#' functions so that they are constrained to the same approximate range of values. Several 
#' function-transformation methods can be used for this end and \pkg{spsann} has four of them available.
#' 
#' The \emph{upper-lower-bound approach} requires the user to inform the maximum (nadir point) and minimum 
#' (utopia point) absolute function values. The resulting function values will always range between 0 and 1.
#' 
#' The \emph{upper-bound approach} requires the user to inform only the nadir point, while the utopia
#' point is set to zero. The upper-bound approach for transformation aims at equalizing only the upper bounds
#' of the objective functions. The resulting function values will always be smaller than or equal to 1.
#' 
#' In most cases, the absolute maximum and minimum values of an objective function cannot be calculated 
#' exactly. If the user is uncomfortable with guessing the nadir and utopia points, there an option
#' for using \emph{numerical simulations}. It consists of computing the function value for many random system
#' configurations. The mean function value obtained over multiple simulations is used to set the nadir point, 
#' while the the utopia point is set to zero. This approach is similar to the upper-bound approach, but the 
#' function values will have the same orders of magnitude only at the starting point of the optimization.
#' Function values larger than one are likely to occur during the optimization. We recommend the user to avoid
#' this approach whenever possible because the effect of the starting configuration on the optimization as a 
#' whole usually is insignificant or arbitrary.
#' 
#' The \emph{upper-lower-bound approach} with the minimum and maximum \emph{attainable} values of the objective
#' functions that compose the MOCOP, also known as the \emph{Pareto minimum and maximum values}, is the most 
#' elegant solution to scale the objective functions. However, it is the most time consuming. It works as
#' follows:
#' 
#' \enumerate{
#' \item Optimize a sample configuration with respect to each objective function that composes the MOCOP;
#' \item Compute the function value of every objective function that composes the MOCOP for every optimized
#' sample configuration;
#' \item Record the minimum and maximum absolute function values attained for each objective function that
#' composes the MOCOP -- these are the Pareto minimum and maximum.
#' }
#' 
#' For example, consider \bold{ACDC}, a MOCOP composed of two objective functions: \bold{CORR} and 
#' \bold{DIST}. The minimum absolute attainable value of \bold{CORR} is obtained when the sample configuration 
#' is optimized with respect to only \bold{CORR}, i.e. when the evaluator and generator objective functions 
#' are the same (see the intersection between the second line and second column in the table below). This is 
#' the Pareto minimum of \bold{CORR}. It follows that the maximum absolute attainable value of \bold{CORR} is 
#' obtained when the sample configuration is optimized with regard to only \bold{DIST}, i.e. when the 
#' evaluator function is difference from the generator function (see the intersection between the first row
#' and the second column in the table below). This is the Pareto maximum of \bold{CORR}. The same logic 
#' applies for finding the Pareto minimum and maximum of \bold{DIST}.
#' 
#' \tabular{rll}{
#' \emph{Evaluator} \tab \emph{Generator} \tab      \cr
#'                  \tab DIST             \tab CORR \cr
#' DIST             \tab 0.5              \tab 8.6  \cr
#' CORR             \tab 6.4              \tab 0.3  \cr
#' }
#' 
#' }
#' 
#' @return 
#' A data frame with the Pareto minimum and maximum values.
#' 
#' @references
#' Arora, J. \emph{Introduction to optimum design}. Waltham: Academic Press, p. 896, 2011.
#'
#' Marler, R. T.; Arora, J. S. Survey of multi-objective optimization methods for engineering. 
#' \emph{Structural and Multidisciplinary Optimization}, v. 26, p. 369-395, 2004.
#' 
#' Marler, R. T.; Arora, J. S. Function-transformation methods for multi-objective optimization. 
#' \emph{Engineering Optimization}, v. 37, p. 551-570, 2005.
#'
#' Marler, R. T.; Arora, J. S. The weighted sum method for multi-objective optimization: new insights.
#' \emph{Structural and Multidisciplinary Optimization}, v. 41, p. 853-862, 2009.
#' 
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[spsann]{optimACDC}}, \code{\link[spsann]{SPAN}}
#' @export
#' @examples 
#' \dontrun{
#' # This example takes more than 5 seconds
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' covars <- meuse.grid[, c(1, 2)]
#' 
#' # CORR
#' schedule <- scheduleSPSANN(initial.acceptance = 0.1, chains = 1, 
#'                            x.max = 1540, y.max = 2060, x.min = 0,
#'                            y.min = 0, cellsize = 40)
#' set.seed(2001)
#' osc_corr <- optimCORR(points = 10, candi = candi, covars = covars, 
#'                       schedule = schedule)
#' 
#' # DIST
#' set.seed(2001)
#' osc_dist <- optimDIST(points = 10, candi = candi, covars = covars,
#'                       schedule = schedule)
#' 
#' # PPL
#' set.seed(2001)
#' osc_ppl <- optimPPL(points = 10, candi = candi, schedule = schedule)
#' 
#' # MSSD
#' set.seed(2001)
#' osc_mssd <- optimMSSD(points = 10, candi = candi, schedule = schedule)
#' 
#' # Pareto
#' pareto <- minmaxPareto(osc = list(DIST = osc_dist, CORR = osc_corr,
#'                                   PPL = osc_ppl, MSSD = osc_mssd),
#'                        candi = candi, covars = covars)
#' pareto
#' }
# FUNCTION - MAIN ##############################################################
minmaxPareto <-
  function (osc, candi, covars) {
    
    obj <- c("CORR", "DIST", "PPL", "MSSD")
    if (!all(names(osc) %in% obj == TRUE)) {
      idx <- which(names(osc) %in% obj == FALSE)
      message(paste("'", names(osc)[idx], "'", " not recognized as a valid name\n", sep = ""))
    } else {
      idx <- match(names(osc), obj)
      osc <- osc[idx]
    }
    
    # Convert numeric covariates into factor covariates
    if (pedometrics::anyFactor(covars) && !pedometrics::allFactor(covars)) {
      id <- which(!sapply(covars, is.factor))
      message(paste("converting ", length(id), " numeric covariates into factor covariates...", sep = ""))
      covars[, id] <- pedometrics::stratify(x = covars[, id], n = nrow(osc[[1]][["points"]]))
    }
    
    # Get objective function parameters
    strata.type <- osc$CORR$objective$strata.type
    use.coords <- osc$CORR$objective$use.coords
    
    # Compute objective function values
    obj_corr <- sapply(1:length(osc), function (i) objCORR(
      osc[[i]]$points, covars = covars, candi = candi, strata.type = strata.type, use.coords = use.coords))
    obj_dist <- sapply(1:length(osc), function (i) objDIST(
      osc[[i]]$points, covars = covars, candi = candi, strata.type = strata.type, use.coords = use.coords))
    
    if (all(c("PPL", "MSSD") %in% names(osc))) {
      
      # Get objective function parameters
      lags <- osc$PPL$objective$lags
      criterion <- osc$PPL$objective$criterion
      pairs <- osc$PPL$objective$pairs
      # x.max <- osc$PPL@spsann$jitter$x[2]
      # x.min <- osc$PPL@spsann$jitter$x[1]
      # y.max <- osc$PPL@spsann$jitter$y[2]
      # y.min <- osc$PPL@spsann$jitter$y[1]
      
      # Compute objective function values
      obj_ppl <- sapply(1:length(osc), function (i) 
        objPPL(osc[[i]]$points, candi = candi, lags = lags, criterion = criterion, pairs = pairs))
               # ,x.max = x.max, y.max = y.max, x.min = x.min, y.min = y.min
      obj_mssd <- sapply(1:length(osc), function (i) objMSSD(osc[[i]]$points, candi = candi))
      
      # Prepare output
      res <- data.frame(
        CORR = obj_corr, DIST = obj_dist, PPL = obj_ppl, MSSD = obj_mssd, 
        row.names = c("CORR", "DIST", "PPL", "MSSD"))
    } else {
      res <- data.frame(CORR = obj_corr, DIST = obj_dist, row.names = c("CORR", "DIST"))
    }
    return (res)
  }
