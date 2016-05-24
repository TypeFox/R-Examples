################################################################################

##########################################
##                                      ##
##    Wrapper for cniperPointPlot.R     ##
##                                      ##
##                                      ##
##########################################

##############################################################
#' Wrapper function for cniperPointPlot - Computation and Plot
#'  of Cniper Contamination and Cniper Points
#'
#' The wrapper takes most of arguments to the cniperPointPlot
#' function by default and gives a user possibility to run the
#' function with low number of arguments
#'
#' @param fam object of class L2ParamFamily
#'
#' @param ... additional parameters (in particular to be passed on to \code{plot})
#'
#' @param lower the lower end point of the contamination interval
#'
#' @param upper the upper end point of the contamination interval
#'
#' @param rescale the flag for rescaling the axes for better view of the plot
#'
#' @param with.legend the flag for showing the legend of the plot
#'
#' @param withCall the flag for the call output
#'
#' @return invisible(NULL)
#
#' @section Details: Calls \code{cniperPointPlot} with suitably chosen defaults; if \code{withCall == TRUE}, the call to \code{cniperPointPlot} is returned.
#'
#' @export
#' @rdname CniperPointPlotWrapper
#'
#' @examples
#' fam <- GammaFamily()
#' CniperPointPlot(fam=fam, main = "Gamma", lower = 0, upper = 5, withCall = FALSE)
##############################################################

##@fam - parameter family
## lower - left point of the x-axis
## upper - right point of the x-axis
## alpha.trsp - optional transparency of the plot
## with.legend - optional legend indicator
## withCall - optional indicator of the function call
#
CniperPointPlot <- function(fam,...
                                  ,lower = getdistrOption("DistrResolution")
                                  ,upper=1-getdistrOption("DistrResolution")
                                  ,with.legend = TRUE
                                  ,rescale = FALSE
                                  ,withCall = TRUE){
  ###
  ### 1. grab the dots (and manipulate it within the wrapper function)
  ###
  ###
  ### do something to fix the good default arguments
  ###
  mc <- as.list(match.call(expand.dots = FALSE))[-1]
  dots <- mc$"..."
  if(is.null(mc$lower)) lower <- getdistrOption("DistrResolution")
  if(is.null(mc$upper)) upper <- 1-getdistrOption("DistrResolution")
  if(is.null(mc$with.legend)) mc$with.legend <- TRUE
  if(is.null(mc$withCall)) mc$withCall <- TRUE
  if(missing(fam)) stop("Argument 'fam' must be given as argument to 'CniperPointPlot'")
  ###
  ### 2. build up the argument list for the (powerful/fullfledged)
  ### graphics/diagnostics function;
  ##

  ## Scaling of the axes
  print(fam)
  scaleList <- rescaleFunction(fam, FALSE, rescale)
  print(scaleList)

  argsList <- c(list(L2Fam = substitute(fam)
                   ,data = substitute(NULL)
                   ,neighbor = substitute(ContNeighborhood(radius = 0.5))
                   ,risk = substitute(asMSE())
                   ,lower = substitute(lower)
                   ,upper = substitute(upper)
                   ,n = substitute(101)
                   ,withMaxRisk = substitute(TRUE)
                   ,scaleN = substitute(9)
                   ,cex.pts = substitute(1)
                   ,col.pts = substitute(par("col"))
                   ,pch.pts = substitute(19)
                   ,jitter.fac = substitute(1)
                   ,with.lab = substitute(FALSE)
                   ,lab.pts = substitute(NULL)
                   ,lab.font = substitute(NULL)
                   ,alpha.trsp = substitute(NA)
                   ,which.lbs = substitute(NULL)
                   ,which.Order  = substitute(NULL)
                   ,return.Order = substitute(FALSE)
                   ,adj = 0.5
                   ,cex.main = substitute(1.5)
                   ,cex.lab = substitute(1.5)
                   ,main = ""#"Outlyingness Plot"
                   ,xlab=substitute("Dirac point")
                   ,ylab=substitute("Asymptotic Risk difference (classic - robust)")
                   ,bty = substitute("o")
                   ), scaleList)
  print(argsList)
  ##parameter for plotting
  if(mc$with.legend)
  {
    argsList$col.main <- "black"
    argsList$col.lab <- "black"
  }
  else
  {
    argsList$col.main <- "white"
    argsList$col.lab <- "white"
  }

  args <- .merge.lists(argsList, dots)
  print(args)
  ###
  ### 3. build up the call but grab it and write it into an object
  ###
  cl <- substitute(do.call(cniperPointPlot,args0), list(args0=args))
  ### manipulate it so that the wrapper do.call is ommitted
  cl0 <- as.list(cl)[-1]
  mycall <- c(cl0[1],unlist(cl0[-1]))
  mycall <- as.call(mycall)
  ###
  ### 4. evaluate the call (i.e., produce the graphic)
  ###
  eval(mycall)
  ###
  ### 5. return the call (if withCall==TRUE)
  ###
  if(mc$withCall) print(mycall)

}
#CniperPointPlot(fam=fam, main = "Gamma", lower = 0, upper = 5, withCall = FALSE)
