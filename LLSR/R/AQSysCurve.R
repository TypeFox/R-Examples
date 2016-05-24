####################################################################################################################
#' @import svDialogs
####################################################################################################################
require(svDialogs)
####################################################################################################################
#' @rdname AQSysCurve
#' @title This functions plot a curve based in the chosen model and its parameters.
#' @description The function returns a plot after using the parameters and model given by the user.
#' @details The function owns predefined set of equations that can be seen below and must be used, with adequated parameters,
#' to return a plot which represent the chosen model.
#' @export AQSysCurve
#' @param mathDesc Equation to be used: merchuk, murugesan [type:string]
#' @param param Model's parameters [type::data.frame]
#' @param xlbl Plot's Horizontal axis label.
#' @param ylbl Plot's Vertical axis label.
#' @param col Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cex Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexlab Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexaxis Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexmain Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexsub Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param xmax Maximum value for the Horizontal axis' value (bottom-rich component)  [type:double]
#' @param mpl Multiples curves overlayed in a single plot. Default is FALSE. [type::LOGIC]
#' @param clwd Plot's axis line width
#' @param save Optimize Plot's elements to be compatible with High Resolution size [type:Boulean]
#' @param HR Magnify Plot's text to be compatible with High Resolution size [type:Boulean]
#' @param NP Number of points used to build the fitted curve. Default is 100. [type:Integer]
#' @inheritParams graphics::plot.default
#' @return A plot using the input model within the chosen interval and the curve's raw XY data.
#' If no interval is selected, xmax = 0.4.
#' @examples
#' \dontrun{
#' AQSysCurve("murugesan", as.data.frame(c(0.90, -3.48, 2.92)), mpl = TRUE, col = "red")
#' }
####################################################################################################################
AQSysCurve <-
  function  (mathDesc, param, xlbl = "", ylbl = "", main = NULL, col = "black", type = "p",
             cex = 1, cexlab = 1, cexaxis = 1, cexmain = 1, cexsub = 1, xmax = 0.4, mpl = FALSE,
             HR = FALSE, NP = 100, clwd = NULL, save = FALSE)
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
    # unlist and convert parameters to double
    param <- as.double(unlist(param))
    # mass fraction range of bottom-rich component (min is 0, max is 1)
    x <- sort(runif(NP, 0.001, xmax))
    # select which model will be used to generate the plot
    switch(
      mathDesc,
      # Fn receives a function correspondent to users choice.
      # check AQSys.mathDesc.R for more details.
      merchuk = {
        Fn <- AQSys.mathDesc("merchuk")
      },
      murugesan = {
        Fn <- AQSys.mathDesc("murugesan")
      },
      tello = {
        Fn <- AQSys.mathDesc("tello")
      },
      # if user selects an option not available, it triggers an error
      # (check AQSys.err.R for details)
      AQSys.err("0")
    )
    # generate curve using selected Fn function and store data in rawdt
    rawdt <- curve(
      Fn(param, x), col = col, add = mpl,
      xlim = c(0, xmax), xlab = xlbl, ylab = ylbl, lwd = clwd
    )
    # set headers for output variable to be returned
    names(rawdt) <- c("XC", "YC")
    # convert output to dataframe variable
    rawdt <- as.data.frame(rawdt)
    # make data available to user
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
