####################################################################################################################
#' @import svDialogs
#' @import ggtern
####################################################################################################################
#' @rdname AQSysTern
#' @title This functions plot a ternary plot based in the chosen model and its parameters.
#' @description The function returns a ternary plot after using the parameters and model given by the user.
#' @details The function have predefined set of equations that can be seen below and must be used, with adequated parameters,
#' to return a plot which represent the chosen model.
#' @export AQSysTern
#' @param XYdt Binodal Experimental data that will be used in the nonlinear fit. [type::data.frame]
#' @param color Set data point's color. All color's names R knows about can be found in \code{\link[grDevices]{colors}}. [type:string]
#' @param shape Set of aesthetic mappings created by \code{\link[ggplot2]{aes}} or \code{\link[ggplot2]{aes_}}.
#' @param title Plot's Title. Default is NULL, for no title. [type:string]
#' @param style Plot's Style.
#' @param xlbl Plot's Bottom-enriched Component axis label. [type:string]
#' @param ylbl Plot's Upper-enriched Component axis label. [type:string]
#' @param zlbl Plot's Water Fraction axis label. [type:string]
#' @param x_arrow_lbl Plot's Bottom-enriched Component arrow label. [type:string]
#' @param y_arrow_lbl Plot's Upper-enriched Component arrow label. [type:string]
#' @param z_arrow_lbl Plot's Water Fraction arrow label. [type:string]
#' @param HR Magnify Plot's text to be compatible with High Resolution size [type:Boulean]
#' @param silent Perform functions taks without returning variables or requesting Input. If TRUE, wdir and filename are mandatory. [type:Boulean]
#' @param save Magnify Plot's text to be compatible with High Resolution size [type:Boulean]
#' @param single TRUE if a single series will be plot and FALSE if otherwise. If FALSE, series variable must be provided. [type:Boulean]
#' @param wdir Set working directory in which plot's file will be saved. Set as "" to save in the current working directory. [type:string]
#' @param filename Set a name for the plot's file. [type:string]
#' @param series A data.frame containin series names in the first columns and its respective colors in the second column. [type:data.frame]
#' @param ltitle plot's legend title. [type:string]
#' @param wlabel Label on ternary arrows. [type:string]
#' @param tcolor T axis color. [type:string]
#' @param lcolor L axis color. [type:string]
#' @param rcolor R axis color. [type:string]
#' @return A ternary plot using the input model within the chosen interval and the curve's raw XY data.
#' @examples
#' \dontrun{
#' AQSysTern(peg4kslt[1:2])
#' }
####################################################################################################################
AQSysTern <-
  function  (XYdt, color = "black", shape = 1, title = NULL, style = "bw", wlabel = "(%, m/m)",
             xlbl = "X", ylbl = "Y", zlbl = "Z", tcolor = "red2", lcolor = "darkgoldenrod2", rcolor = "blue4",
             x_arrow_lbl = xlbl, y_arrow_lbl = ylbl, z_arrow_lbl = zlbl,
             HR = FALSE, silent = FALSE, save = silent, single = TRUE,
             wdir = NULL, filename = NULL, series = NULL, ltitle = NULL)
  {
    # Solution for "no visible binding for global variable VAR" and "no visible global function definition for"
    x <- y <- z <- serie <- NULL
    labs <-
      geom_point <-
      element_line <-
      element_line <-
      scale_shape_manual <-
      scale_colour_manual <- element_text <- NULL
    #
    #
    #
    #
    #
    if ((ncol(XYdt) != 2) && (single  == TRUE)) {
      #
      stop("Dataset parameter must be a data.frame, or list, with exactly two columns.")
      #
    } else if ((silent == TRUE) &&
               (is.null(wdir) || is.null(filename))) {
      #
      stop("Input parameters filename and wdir must be set if silent == TRUE.")
      #
    } else if (single == FALSE) {
      #
      if (ncol(XYdt) <= 2) {
        stop("At leas two (x,y) pairs are necessary to plot multiple curves (single == FALSE).")
      }
      #
      n_sys_pairs <- ncol(XYdt[3:length(XYdt)]) / 2
      #
      if ((n_sys_pairs * 2) %% 2 != 0) {
        #
        stop(
          "Number of data columns must be a multiple of two. A (x,y) pair is necessary for each curve to be plotted."
        )
        #
      } else {
        #
        data.temp <- NULL
        data.melted <-
          data.frame(
            x = double(), y = double(), z = double(), serie = character()
          )
        #
        if (is.null(series)) {
          srn <- NULL
          crn <- NULL
          shn <- NULL
          #
          for (i in 1:(n_sys_pairs + 1)) {
            # Serie Name
            srn <- c(srn, paste("Serie", i))
            # Colour Number
            crn <- c(crn, i)
            # Shape Number
            shn <- c(shn, i * 3)
          }
          #
          series <- data.frame(srn, crn, shn)
        }
        #
        for (i in seq(1, length(XYdt), 2)) {
          #
          data.temp <-
            data.frame(XYdt[i:(i + 1)], 1 - rowSums(XYdt[i:(i + 1)]), series[((i + 1) / 2), 1])
          names(data.temp) <- c("x", "y", "z", "serie")
          data.melted <- rbind(data.melted, data.temp)
          #
        }
        tern_data <- data.melted
        #
        series_colors <- t(series[2])
        names(series_colors) <- t(series[1])
        series_shapes <- t(series[3])
        names(series_shapes) <- names(series_colors)
      }
    } else if (single == TRUE) {
      #
      tern_data <- data.frame(XYdt[1:2], 1 - rowSums(XYdt[1:2]))
      names(tern_data) <- c("x", "y", "z")
      #
    }
    #
    #
    #
    #
    #
    if (save == TRUE) {
      if (HR == TRUE) {
        image_format <- ".svg"
      }else{
        image_format <- ".png"
      }
      #
      if (is.null(filename)) {
        # Get user choice for a filename to save the plot
        filename <-
          dlgInput(message = "Enter the figure filename:")$res
      }
      # complete filename with the appropriated extension
      filename <- paste(filename, image_format, sep = "")
      # Check if filename is invalid and quite if so
      if (filename == image_format) {
        stop("Filename is NULL or INVALID.", call. = TRUE)
      }
      #
      #
      if (is.null(wdir)) {
        # Get user choice for a directory to save the plot
        wdir <- dlgDir()$res
      }
      # Check if path is invalid and quite if so
      if ((wdir == "") && (silent == FALSE)) {
        #
        stop("Path is NULL or INVALID.", call. = TRUE)
        #
      }else if ((wdir == "") && (silent == TRUE)) {
        #
        wdir <- getwd()
        wdir <- paste(wdir, filename, sep = .Platform$file.sep)
        #
      }else{
        #
        wdir <- paste(wdir, filename, sep = .Platform$file.sep)
        #
      }
    }
    #
    #
    #
    #
    #
    components_limits <- round(max(tern_data[2]) / .92, 2)
    #
    tern_image <-
      ggtern(tern_data, aes(z, x, y)) +
      theme_light() +
      tern_limits(T = components_limits,L = 1,R = components_limits) +
      Tlab(xlbl) +
      Llab(ylbl) +
      Rlab(zlbl) +
      Tarrowlab(x_arrow_lbl) +
      Larrowlab(y_arrow_lbl) +
      Rarrowlab(z_arrow_lbl) +
      theme_showarrows() +
      Wlab(wlabel) +
      labs(title = title) +
      theme(legend.position = "bottom")
    #
    if (single == FALSE) {
      tern_image <-
        tern_image + geom_point(size = 2, aes(shape = serie, colour = serie)) +
        # scale_shape_discrete(name=ltitle) +
        # scale_colour_discrete(name=ltitle) +
        scale_colour_manual(name = ltitle, values = series_colors) +
        scale_shape_manual(name = ltitle, values = series_shapes)
      # scale_colour_manual(values = series_colors) +
      #theme(legend.title = element_blank())
    } else {
      tern_image <- tern_image + geom_point(shape = shape, colour = color)
    }
    #
    #
    #
    #
    #
    if (style == "bw") {
      tern_image <- tern_image +
        theme(
          text = element_text(size = 18),
          tern.axis.line = element_line(
            size = 2, linetype = 1, color = "black"
          ),
          tern.axis.arrow = element_line(color = "black"),
          tern.axis.ticks = element_line(color = "black"),
          tern.axis.title = element_text(color = "black"),
          tern.axis.arrow.text = element_text(color = "black"),
          tern.axis.text = element_text(color = "black"),
          tern.panel.grid = element_line(
            size = 1, linetype = 2, colour = "black"
          ),
          tern.panel.grid.major = element_line(color = "black"),
          tern.panel.grid.minor = element_line(color = "black")
        )
    } else if (style == "custom") {
      tern_image <- tern_image +
        theme(
          tern.axis.line = element_line(size = 2, linetype = 1),
          tern.axis.line.T = element_line(color = tcolor),
          tern.axis.line.L = element_line(color = lcolor),
          tern.axis.line.R = element_line(color = rcolor),
          
          tern.axis.title.T = element_text(color = tcolor),
          tern.axis.title.L = element_text(color = lcolor),
          tern.axis.title.R = element_text(color = rcolor),
          
          tern.axis.arrow.T = element_line(color = tcolor),
          tern.axis.arrow.L = element_line(color = lcolor),
          tern.axis.arrow.R = element_line(color = rcolor),
          
          tern.axis.arrow.text.T = element_text(color = tcolor),
          tern.axis.arrow.text.L = element_text(color = lcolor),
          tern.axis.arrow.text.R = element_text(color = rcolor),
          
          tern.axis.text.T = element_text(color = tcolor),
          tern.axis.text.L = element_text(color = lcolor),
          tern.axis.text.R = element_text(color = rcolor),
          
          tern.panel.grid.minor.T = element_line(
            color = tcolor, size = .1, linetype = 2
          ),
          tern.panel.grid.minor.L = element_line(
            color = lcolor, size = .1, linetype = 2
          ),
          tern.panel.grid.minor.R = element_line(
            color = rcolor, size = .1, linetype = 2
          ),
          
          tern.panel.grid.major.T = element_line(
            color = tcolor, size = .1, linetype = 2
          ),
          tern.panel.grid.major.L = element_line(
            color = lcolor, size = .1, linetype = 2
          ),
          tern.panel.grid.major.R = element_line(
            color = rcolor, size = .1, linetype = 2
          )
        )
    } else if (style == "rgbw") {
      tern_image <- tern_image + theme_rgbw() +
        theme(
          tern.axis.line = element_line(size = 2, linetype = 1),
          tern.panel.grid = element_line(size = 1, linetype = 2)
        )
    } else{
      tern_image <- tern_image + theme_bw() +
        theme(
          tern.axis.line = element_line(size = 2, linetype = 1),
          tern.panel.grid = element_line(size = 1, linetype = 2)
        )
    }
    #
    #
    #
    #
    #
    if (silent == FALSE) {
      print(tern_image)
    }
    if (save == TRUE) {
      ggsave(
        filename = wdir, plot = tern_image, width = 21.14 / 2, height = 14.39 / 2
      )
    }
    #
  }
#
