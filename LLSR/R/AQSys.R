####################################################################################################################
#' @import rootSolve
#' @import graphics
#' @import stats
#' @import svDialogs
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
####################################################################################################################
require(rootSolve)
require(svDialogs)
# options(warn = 1)
####################################################################################################################
#'Merchuk's Equation to fit Binodal Experimental Data
#' @rdname AQSys
#' @name AQSys
#' @description .
#' @export
#' @seealso \itemize{
#' \item \code{\link{AQSys.default}}
#' \item \code{\link{AQSys.plot}}
#' \item \code{\link{AQSys.tielines}}
#' \item \code{\link{AQSys.crpt}}
#' \item \code{\link{AQSysOthmer}}
#' \item \code{\link{AQSysBancroft}}
#' }
####################################################################################################################
AQSys <- function(XYdt,...)
  UseMethod("AQSys")
####################################################################################################################
#' @rdname AQSys
#' @title Merchuk's nonlinear Equation
#' @description Perform a nonlinear regression fit using several mathemmatical descriptors in order to determine the
#' equation's parameters.
#' @details The function returns functions parameters after fitting experimental data to the equations listed in AQSysList().
#' @param mathDesc - Character String specifying the nonlinear empirical equation to fit data.
#' The default method uses Merchuk's equation. Other mathematical descriptors can be listed using AQSysList().
#' @param XYdt - Binodal Experimental data that will be used in the nonlinear fit
#' @param ... Additional optional arguments. None are used at present.
#' @method AQSys default
#' @export
#' @return A list containing three data.frame variables with all data parsed from the worksheet and parameters calculated
#' through the available mathematical descriptions.
#' @examples
#' #Populating variable XYdt with binodal data
#' XYdt <- peg4kslt[,1:2]
#' #Fitting XYdt using Merchuk's function
#' AQSys(XYdt)
AQSys.default <- function(XYdt, mathDesc = "merchuk",...) {
  # each switch option calls a correspondent equation to fit XYdt
  # equations are functions declared in AQSysFormulas.R
  switch(
    mathDesc,
    merchuk = {
      ans <- mrchk(XYdt)
    },
    murugesan = {
      ans <- mrgsn(XYdt)
    },
    tello = {
      ans <- tello(XYdt)
    },
    AQSys.err("0")
  )
  # return fitting parameters ans statistical data
  return(ans)
}
####################################################################################################################
# MERCHUK PLOT TEST FUNCTION
#' @rdname AQSys.plot
#' @title Dataset and Fitted Function plot
#' @description The function returns a plot after fitting a dataset to a given equation.
#' @details This version uses the plot function and return a regular bidimensional plot.
#' @method AQSys plot
#' @export AQSys.plot
#' @export
#' @param mathDesc - Character String specifying the nonlinear empirical equation to fit data. The default method uses
#' Merchuk's equation. Other possibilities can be seen in AQSysList().
#' @param ... Additional optional arguments. None are used at present.
#' @param XYdt Binodal Experimental data that will be used in the nonlinear fit
#' @param xlbl Plot's Horizontal axis label.
#' @param ylbl Plot's Vertical axis label.
#' @param main Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param col Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param type Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cex Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexlab Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexaxis Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexmain Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexsub Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param xmax Maximum value for the Horizontal axis' value
#' @param ymax Maximum value for the Vertical axis' value
#' @param HR Magnify Plot's text to be compatible with High Resolution size [type:Boulean]
#' @param NP Number of points used to build the fitted curve. Default is 100. [type:Integer]
#' @param clwd Plot's axis line width
#' @param save Optimize Plot's elements to be compatible with High Resolution size [type:Boulean]
#' @return A plot containing the experimental data, the correspondent curve for the binodal in study and the curve's raw XY data.
#' @examples
#' #Populating variable XYdt with binodal data
#' XYdt <- peg4kslt[,1:2]
#' #Plot XYdt using Merchuk's function
#' #
#' AQSys.plot(XYdt)
#' #
AQSys.plot <-
  function  (XYdt, xlbl = "", ylbl = "", main = NULL, col = "blue", type = "p",
             cex = 1, cexlab = 1, cexaxis = 1, cexmain = 1, cexsub = 1,
             xmax = 0.4, ymax = 0.5, HR = FALSE, NP = 100, mathDesc = "merchuk",
             clwd = NULL, save = FALSE, ...)
  {
    #
    if (save == TRUE) {
      # Open dialog to get filename string
      filename <-
        paste(dlgInput(message = "Enter the figure filename:")$res, ".png", sep = "")
      # Check if filename is invalid and quite if so
      if (filename == ".png") {
        stop("Filename is NULL or INVALID.", call. = TRUE)
      }
      # Get user choice for a directory to save the plot
      savePath <- dlgDir()$res
      # Check if path is invalid and quite if so
      if (savePath == "") {
        stop("Path is NULL or INVALID.", call. = TRUE)
      }else{
        savePath <- paste(savePath, filename, sep = .Platform$file.sep)
      }
      # set plot resolution based on user choice
      if (HR == TRUE) {
        png(savePath, width = 5235, height = 3240, units = "px")
      }else{
        png(savePath, width = 1745, height = 1080, units = "px")
      }
    }
    #
    # set graph parameters to export plots in High Quality
    if (is.null(clwd)) {
      clwd <- AQSysHR(HR)
    }else{
      AQSysHR(HR)
    }
    # select which model will be used to generate the plot
    switch(
      mathDesc,
      merchuk = {
        # fit data using chosen equation and get coefficients
        CoefSET <- summary(mrchk(XYdt))$coefficients[, 1]
        # set the equation that will be used to plot the phase diagram
        Fn <- AQSys.mathDesc("merchuk")
      },
      murugesan = {
        # fit data using chosen equation and get coefficients
        CoefSET <- summary(mrgsn(XYdt))$coefficients[, 1]
        # set the equation that will be used to plot the phase diagram
        Fn <- AQSys.mathDesc("murugesan")
      },
      tello = {
        # fit data using chosen equation and get coefficients
        CoefSET <- summary(tello(XYdt))$coefficients[, 1]
        # set the equation that will be used to plot the phase diagram
        Fn <- AQSys.mathDesc("tello")
      },
      # if user selects an option not available, it triggers an error (check AQSys.err.R for details)
      AQSys.err("0")
    )
    #plot phase diagram using experimental data and with previously selected parameters
    plot(
      XYdt, xlab = xlbl, ylab = ylbl, main = main, col = col, type = type,
      cex = cex, cex.lab = cexlab, cex.axis = cexaxis, cex.main = cexmain,
      cex.sub = cexsub, xlim = c(0,xmax), ylim = c(0,ymax)
    )
    # mass fraction range of bottom-rich component (min is 0, max is 1)
    x <- sort(runif(NP,0.001,xmax))
    #
    #SWITCH WORKS ONLY FOR THREE PARAMETER'S EQUATIONS.
    #IF NECESSARY A HIGHER NUMBER, INSERT CONDITIONAL BELOW.
    #Maybe change it to have as input the whole coefficient set?
    #a<-summary(mrchk(peg4kslt[,1:2]))$coefficients[,1]
    #
    # add curve generated using regression parameters
    rawdt <- curve(Fn(CoefSET,x),
                   add = TRUE, n = NP)
    names(rawdt) <- c("XC","YC")
    rawdt <- as.data.frame(rawdt)
    # make available data from fitted curve to user. Function returns it silently
    # but user can get data using simple assign '<-'
    invisible(rawdt)
    # Set ticks Thickness
    axis(side = 1, lwd = clwd)
    axis(side = 2, lwd = clwd)
    #
    if (save == TRUE) {
      invisible(dev.off())
    }
    #
  }
####################################################################################################################
