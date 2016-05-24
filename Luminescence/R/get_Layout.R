#' Collection of layout definitions
#'
#' This helper function returns a list with layout definitions for homogeneous
#' plotting.
#'
#' The easiest way to create a user-specific layout definition is perhaps to
#' create either an empty or a default layout object and fill/modify the
#' definitions (\code{user.layout <- get_Layout(data = "empty")}).
#'
#' @param layout \code{\link{character}} or \code{\link{list}} object
#' (required): name of the layout definition to be returned. If name is
#' provided the respective definition is returned. One of the following
#' supported layout definitions is possible: \code{"default"},
#' \code{"journal.1"}, \code{"small"}, \code{"empty"}. User-specific layout
#' definitions must be provided as a list object of predefined structure, see
#' details.
#' @return A list object with layout definitions for plot functions.
#' @section Function version: 0.1
#' @author Michael Dietze, GFZ Potsdam (Germany)
#' @examples
#'
#' ## read example data set
#' data(ExampleData.DeValues, envir = environment())
#'
#' ## show structure of the default layout definition
#' layout.default <- get_Layout(layout = "default")
#' str(layout.default)
#'
#' ## show colour definitions for Abanico plot, only
#' layout.default$abanico$colour
#'
#' ## set Abanico plot title colour to orange
#' layout.default$abanico$colour$main <- "orange"
#'
#' ## create Abanico plot with modofied layout definition
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  layout = layout.default)
#'
#' ## create Abanico plot with predefined layout "journal"
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  layout = "journal")
#'
#' @export
get_Layout <- function(
  layout
) {

  ## pre-defined layout selections
  if(is.character(layout) == TRUE & length(layout) == 1) {

    if(layout == "empty") {

      # empty layout definition -------------------------------------------------
      layout = list(
        abanico = list(
          font.type = list(
            main    = character(1),
            xlab1   = character(1),
            xlab2   = character(1),
            ylab    = character(1),
            zlab    = character(1),
            xtck1   = character(1),
            xtck2   = character(1),
            xtck3   = character(1),
            ytck    = character(1),
            ztck    = character(1),
            mtext   = character(1),
            summary = character(1), # optionally vector
            stats   = character(1), # optionally vector
            legend  = character(1) # optionally vector
          ),
          font.size = list(
            main    = numeric(1),
            xlab1   = numeric(1),
            xlab2   = numeric(1),
            xlab3   = numeric(1),
            ylab    = numeric(1),
            zlab    = numeric(1),
            xtck1   = numeric(1),
            xtck2   = numeric(1),
            xtck3   = numeric(1),
            ytck    = numeric(1),
            ztck    = numeric(1),
            mtext   = numeric(1),
            summary = numeric(1), # optionally vector
            stats   = numeric(1), # optionally vector
            legend  = numeric(1)  # optionally vector
          ),
          font.deco = list(
            main    = character(1),
            xlab1   = character(1),
            xlab2   = character(1),
            xlab3   = character(1),
            ylab    = character(1),
            zlab    = character(1),
            xtck1   = character(1),
            xtck2   = character(1),
            xtck3   = character(1),
            ytck    = character(1),
            ztck    = character(1),
            mtext   = character(1),
            summary = character(1), # optionally vector
            stats   = character(1), # optionally vector
            legend  = character(1) # optionally vector
          ),
          colour = list(
            main    = numeric(1), # plot title colour
            xlab1   = numeric(1), # left x-axis label colour
            xlab2   = numeric(1), # right x-axis label colour
            xlab3   = numeric(1), # right x-axis label colour
            ylab    = numeric(1), # y-axis label colour
            zlab    = numeric(1), # z-axis label colour
            xtck1   = numeric(1), # left x-axis tick colour
            xtck2   = numeric(1), # right x-axis tick colour
            xtck3   = numeric(1), # right x-axis tick colour
            ytck    = numeric(1), # y-axis tick colour
            ztck    = numeric(1), # z-axis tick colour
            mtext   = numeric(1), # subheader text colour
            summary = numeric(1), # statistic summary colour
            stats   = numeric(1), # value statistics colour
            legend  = numeric(1), # legend colour
            centrality = numeric(1), # Centrality line colour
            value.dot  = numeric(1), # De value dot colour
            value.bar  = numeric(1), # De value error bar colour
            value.rug  = numeric(1), # De value rug colour
            poly.line  = numeric(1), # polygon line colour
            poly.fill  = numeric(1), # polygon fill colour
            bar.line   = numeric(1), # polygon line colour
            bar.fill   = numeric(1), # polygon fill colour
            kde.line   = numeric(1),
            kde.fill   = numeric(1),
            grid.major = numeric(1),
            grid.minor = numeric(1),
            border     = numeric(1),
            background = numeric(1)),
          dimension = list(
            figure.width    = numeric(1), # figure width in mm
            figure.height   = numeric(1), # figure height in mm
            margin          = numeric(4), # margin sizes in mm
            main.line       = numeric(1), # line height in %
            xlab1.line      = numeric(1), # line height in %
            xlab2.line      = numeric(1), # line height in %
            xlab3.line      = numeric(1), # line height in %
            ylab.line       = numeric(1), # line height in %
            zlab.line       = numeric(1), # line height in %
            xtck1.line      = numeric(1), # line height in %
            xtck2.line      = numeric(1), # line height in %
            xtck3.line      = numeric(1), # line height in %
            ytck.line       = numeric(1), # line height in %
            ztck.line       = numeric(1), # line height in %
            xtcl1           = numeric(1), # tick length in %
            xtcl2           = numeric(1), # tick length in %
            xtcl3           = numeric(1), # tick length in %
            ytcl            = numeric(1), # tick length in %
            ztcl            = numeric(1), # tick length in %
            rugl            = numeric(1), # rug length in %
            mtext           = numeric(1), # line height in %
            summary.line    = numeric(1) # line height in %
          )),
        kde = list(
          font.type = list(
            ## run font_import() when using user-defined fonts for the first time
            main   = character(1),
            xlab   = character(1),
            ylab1  = character(1),
            ylab2  = character(1),
            xtck   = character(1),
            ytck1  = character(1),
            ytck2  = character(1),
            stats  = character(1), # optionally vector
            legend = character(1) # optionally vector
          ),
          font.size = list(
            main   = numeric(1),
            xlab   = numeric(1),
            ylab1  = numeric(1),
            ylab2  = numeric(1),
            xtck   = numeric(1),
            ytck1  = numeric(1),
            ytck2  = numeric(1),
            stats  = numeric(1), # optionally vector
            legend = numeric(1) # optionally vector
          ),
          font.deco = list(
            main   = character(1),
            xlab   = character(1),
            ylab1  = character(1),
            ylab2  = character(1),
            xtck   = character(1),
            ytck1  = character(1),
            ytck2  = character(1),
            stats  = character(1), # optionally vector
            legend = character(1) # optionally vector
          ),
          colour = list(
            main   = numeric(1), # plot title colour
            xlab   = numeric(1), # x-axis label colour
            ylab1  = numeric(1), # primary y-axis label colour
            ylab2  = numeric(1), # secondary y-axis label colour
            xtck   = numeric(1), # x-axis tick colour
            ytck1  = numeric(1), # primary y-axis tick colour
            ytck2  = numeric(1), # secondary y-axis tick colour
            box    = numeric(1), # plot frame box line colour
            mtext  = numeric(1), # subheader text colour
            stats  = numeric(1), # statistic summary colour
            kde.line        = numeric(1), # KDE line colour
            kde.fill        = numeric(1), # KDE fill colour
            value.dot       = numeric(1), # De value dot colour
            value.bar       = numeric(1), # De value error bar colour
            value.rug       = numeric(1), # De value rug colour
            mean            = numeric(1), # mean line colour
            median          = numeric(1), # median line colour
            mean.weighted   = numeric(1), # weighted mean line colour
            median.weighted = numeric(1), # weighted median line colour
            kdemax          = numeric(1), # KDE max line colour
            poly.line       = numeric(1), # polygon line colour
            poly.fill       = numeric(1), # polygon fill colour
            background      = numeric(1)),
          dimension = list(
            figure.width    = numeric(1), # figure width in mm
            figure.height   = numeric(1), # figure height in mm
            margin          = numeric(4), # margin sizes in mm
            main.line       = numeric(1), # line height in %
            xlab.line       = numeric(1), # line height in %
            ylab1.line      = numeric(1), # line height in %
            ylab2.line      = numeric(1), # line height in %
            xtck.line       = numeric(1), # line height in %
            ytck1.line      = numeric(1), # line height in %
            ytck2.line      = numeric(1), # line height in %
            xtcl            = numeric(1), # tick length in %
            ytcl1           = numeric(1), # tick length in %
            ytcl2           = numeric(1), # tick length in %
            stats.line      = numeric(1) # line height in %
          )
        )
      )
    } else if(layout == "default") {

      # default layout definition -----------------------------------------------
      layout = list(
        abanico = list(
          font.type = list(
            main    = "",
            xlab1   = "",
            xlab2   = "",
            ylab    = "",
            zlab    = "",
            xtck1   = "",
            xtck2   = "",
            xtck3   = "",
            ytck    = "",
            ztck    = "",
            mtext   = "",
            summary = "", # optionally vector
            stats   = "", # optionally vector
            legend  = "" # optionally vector
          ),
          font.size = list(
            main    = 12,
            xlab1   = 12,
            xlab2   = 12,
            xlab3   = 12,
            ylab    = 12,
            zlab    = 12,
            xtck1   = 12,
            xtck2   = 12,
            xtck3   = 12,
            ytck    = 12,
            ztck    = 12,
            mtext   = 10,
            summary = 10, # optionally vector
            stats   = 10, # optionally vector
            legend  = 10 # optionally vector
          ),
          font.deco = list(
            main    = "bold",
            xlab1   = "normal",
            xlab2   = "normal",
            xlab3   = "normal",
            ylab    = "normal",
            zlab    = "normal",
            xtck1   = "normal",
            xtck2   = "normal",
            xtck3   = "normal",
            ytck    = "normal",
            ztck    = "normal",
            mtext   = "normal",
            summary = "normal", # optionally vector
            stats   = "normal", # optionally vector
            legend  = "normal" # optionally vector
          ),
          colour = list(
            main    = 1, # plot title colour
            xlab1   = 1, # left x-axis label colour
            xlab2   = 1, # right x-axis label colour
            xlab3   = 1, # right x-axis label colour
            ylab    = 1, # y-axis label colour
            zlab    = 1, # z-axis label colour
            xtck1   = 1, # left x-axis tick colour
            xtck2   = 1, # right x-axis tick colour
            xtck3   = 1, # right x-axis tick colour
            ytck    = 1, # y-axis tick colour
            ztck    = 1, # z-axis tick colour
            mtext   = 1, # subheader text colour
            summary = 1, # statistic summary colour
            stats   = 1, # value statistics colour
            legend  = 1, # legend colour
            centrality = 1, # Centrality line colour
            value.dot  = 1, # De value dot colour
            value.bar  = 1, # De value error bar colour
            value.rug = 1, # De value rug colour
            poly.line  = NA, # polygon line colour
            poly.fill  = adjustcolor("grey75", alpha.f = 0.6), # polygon fill colour
            bar.line   = NA, # polygon line colour
            bar.fill   = "grey60", # bar fill colour
            kde.line   = 1,
            kde.fill   = NA,
            grid.major = "grey80",
            grid.minor = "none",
            border     = 1,
            background = NA),
          dimension = list(
            figure.width    = "auto", # figure width in mm
            figure.height   = "auto", # figure height in mm
            margin = c(10, 10, 10, 10), # margin sizes in mm
            main.line       = 100, # line height in %
            xlab1.line      = 90, # line height in %
            xlab2.line      = 90, # line height in %
            xlab3.line      = 90, # line height in %
            ylab.line       = 100, # line height in %
            zlab.line       = 70, # line height in %
            xtck1.line      = 100, # line height in %
            xtck2.line      = 100, # line height in %
            xtck3.line      = 100, # line height in %
            ytck.line       = 100, # line height in %
            ztck.line       = 100, # line height in %
            xtcl1           = 100, # tick length in %
            xtcl2           = 100, # tick length in %
            xtcl3           = 100, # tick length in %
            ytcl            = 100, # tick length in %
            ztcl            = 100, # tick length in %
            rugl            = 100, # rug length in %
            mtext           = 100, # line height in %
            summary.line    = 100 # line height in %
          )),
        kde = list(
          font.type = list(
            ## run font_import() when using user-defined fonts for the first time
            main   = "",
            xlab   = "",
            ylab1  = "",
            ylab2  = "",
            xtck   = "",
            ytck1  = "",
            ytck2  = "",
            stats  = "", # optionally vector
            legend = "" # optionally vector
          ),
          font.size = list(
            main   = 14,
            xlab   = 12,
            ylab1  = 12,
            ylab2  = 12,
            xtck   = 12,
            ytck1  = 12,
            ytck2  = 12,
            stats  = 12, # optionally vector
            legend = 12 # optionally vector
          ),
          font.deco = list(
            main   = "bold",
            xlab   = "normal",
            ylab1  = "normal",
            ylab2  = "normal",
            xtck   = "normal",
            ytck1  = "normal",
            ytck2  = "normal",
            stats  = "normal", # optionally vector
            legend = "normal" # optionally vector
          ),
          colour = list(
            main   = 1, # plot title colour
            xlab   = 1, # x-axis label colour
            ylab1  = 1, # primary y-axis label colour
            ylab2  = 1, # secondary y-axis label colour
            xtck   = 1, # x-axis tick colour
            ytck1  = 1, # primary y-axis tick colour
            ytck2  = 1, # secondary y-axis tick colour
            box    = 1, # plot frame box line colour
            mtext  = 1, # subheader text colour
            stats  = "#2062B3", # statistic summary colour
            kde.line        = "#2062B3", # KDE line colour
            kde.fill        = "none", # KDE fill colour
            value.dot       = 1, # De value dot colour
            value.bar       = 1, # De value error bar colour
            value.rug       = 1, # De value rug colour
            mean   = 1, # mean line colour
            median = 1, # median line colour
            mean.weighted   = 1, # weighted mean line colour
            median.weighted = 1, # weighted median line colour
            kdemax = 1, # KDE max line colour
            poly.line       = NA, # polygon line colour
            poly.fill       = "grey80", # polygon fill colour
            background      = NULL),
          dimension = list(
            figure.width    = 50, # figure width in mm
            figure.height   = 50, # figure height in mm
            margin = c(10, 10, 10, 10), # margin sizes in mm
            main.line       = 100, # line height in %
            xlab.line       = 100, # line height in %
            ylab1.line      = 100, # line height in %
            ylab2.line      = 100, # line height in %
            xtck.line       = 100, # line height in %
            ytck1.line      = 100, # line height in %
            ytck2.line      = 100, # line height in %
            xtcl            = 100, # tick length in %
            ytcl1           = 100, # tick length in %
            ytcl2           = 100, # tick length in %
            stats.line      = 100 # line height in %
          )
        )
      )
    } else if(layout == "journal") {

      # journal layout definition -----------------------------------------------
      layout = list(
        abanico = list(
          font.type = list(
            main    = "",
            xlab1   = "",
            xlab2   = "",
            ylab    = "",
            zlab    = "",
            xtck1   = "",
            xtck2   = "",
            xtck3   = "",
            ytck    = "",
            ztck    = "",
            mtext   = "",
            summary = "", # optionally vector
            stats   = "", # optionally vector
            legend  = "" # optionally vector
          ),
          font.size = list(
            main    = 8,
            xlab1   = 7,
            xlab2   = 7,
            xlab3   = 7,
            ylab    = 7,
            zlab    = 7,
            xtck1   = 7,
            xtck2   = 7,
            xtck3   = 7,
            ytck    = 7,
            ztck    = 7,
            mtext   = 6,
            summary = 6, # optionally vector
            stats   = 6, # optionally vector
            legend  = 6 # optionally vector
          ),
          font.deco = list(
            main    = "bold",
            xlab1   = "normal",
            xlab2   = "normal",
            xlab3   = "normal",
            ylab    = "normal",
            zlab    = "normal",
            xtck1   = "normal",
            xtck2   = "normal",
            xtck3   = "normal",
            ytck    = "normal",
            ztck    = "normal",
            mtext   = "normal",
            summary = "normal", # optionally vector
            stats   = "normal", # optionally vector
            legend  = "normal" # optionally vector
          ),
          colour = list(
            main    = 1, # plot title colour
            xlab1   = 1, # left x-axis label colour
            xlab2   = 1, # right x-axis label colour
            xlab3   = 1, # right x-axis label colour
            ylab    = 1, # y-axis label colour
            zlab    = 1, # z-axis label colour
            xtck1   = 1, # left x-axis tick colour
            xtck2   = 1, # right x-axis tick colour
            xtck3   = 1, # right x-axis tick colour
            ytck    = 1, # y-axis tick colour
            ztck    = 1, # z-axis tick colour
            mtext   = 1, # subheader text colour
            summary = 1, # statistic summary colour
            stats   = 1, # value statistics colour
            legend  = 1, # legend colour
            centrality = 1, # Centrality line colour
            value.dot  = 1, # De value dot colour
            value.bar  = 1, # De value error bar colour
            value.rug  = 1, # De value rug colour
            poly.line  = NA, # polygon line colour
            poly.fill  = adjustcolor("grey75", alpha.f = 0.6), # polygon fill colour
            bar.line   = NA, # polygon line colour
            bar.fill   = "grey60", # bar fill colour
            kde.line   = 1,
            kde.fill   = NA,
            grid.major = "grey80",
            grid.minor = "none",
            border     = 1,
            background = NA),
          dimension = list(
            figure.width    = 100, # figure width in mm
            figure.height   = 100, # figure height in mm
            margin = c(10, 10, 10, 10), # margin sizes in mm
            main.line       = 70, # line height in %
            xlab1.line      = 30, # line height in %
            xlab2.line      = 65, # line height in %
            xlab3.line      = 30, # line height in %
            ylab.line       = 30, # line height in %
            zlab.line       = 40, # line height in %
            xtck1.line      = 50, # line height in %
            xtck2.line      = 50, # line height in %
            xtck3.line      = 50, # line height in %
            ytck.line       = 70, # line height in %
            ztck.line       = 70, # line height in %
            xtcl1           = 50, # tick length in %
            xtcl2           = 50, # tick length in %
            xtcl3           = 50, # tick length in %
            ytcl            = 50, # tick length in %
            ztcl            = 70, # tick length in %
            rugl            = 70, # rug length in %
            mtext           = 100, # line height in %
            summary.line    = 70, # line height in %
            pch             = 50  # point size in %
          )),
        kde = list(
          font.type = list(
            ## run font_import() when using user-defined fonts for the first time
            main   = "",
            xlab   = "",
            ylab1  = "",
            ylab2  = "",
            xtck   = "",
            ytck1  = "",
            ytck2  = "",
            stats  = "", # optionally vector
            legend = "" # optionally vector
          ),
          font.size = list(
            main   = 14,
            xlab   = 12,
            ylab1  = 12,
            ylab2  = 12,
            xtck   = 12,
            ytck1  = 12,
            ytck2  = 12,
            stats  = 12, # optionally vector
            legend = 12 # optionally vector
          ),
          font.deco = list(
            main   = "bold",
            xlab   = "normal",
            ylab1  = "normal",
            ylab2  = "normal",
            xtck   = "normal",
            ytck1  = "normal",
            ytck2  = "normal",
            stats  = "normal", # optionally vector
            legend = "normal" # optionally vector
          ),
          colour = list(
            main   = 1, # plot title colour
            xlab   = 1, # x-axis label colour
            ylab1  = 1, # primary y-axis label colour
            ylab2  = 1, # secondary y-axis label colour
            xtck   = 1, # x-axis tick colour
            ytck1  = 1, # primary y-axis tick colour
            ytck2  = 1, # secondary y-axis tick colour
            box    = 1, # plot frame box line colour
            mtext  = 1, # subheader text colour
            stats  = "#2062B3", # statistic summary colour
            kde.line        = "#2062B3", # KDE line colour
            kde.fill        = "none", # KDE fill colour
            value.dot       = 1, # De value dot colour
            value.bar       = 1, # De value error bar colour
            value.rug       = 1, # De value rug colour
            mean   = 1, # mean line colour
            median = 1, # median line colour
            mean.weighted   = 1, # weighted mean line colour
            median.weighted = 1, # weighted median line colour
            kdemax = 1, # KDE max line colour
            poly.line       = NA, # polygon line colour
            poly.fill       = "grey80", # polygon fill colour
            background      = NULL),
          dimension = list(
            figure.width    = 50, # figure width in mm
            figure.height   = 50, # figure height in mm
            margin = c(10, 10, 10, 10), # margin sizes in mm
            main.line       = 100, # line height in %
            xlab.line       = 100, # line height in %
            ylab1.line      = 100, # line height in %
            ylab2.line      = 100, # line height in %
            xtck.line       = 100, # line height in %
            ytck1.line      = 100, # line height in %
            ytck2.line      = 100, # line height in %
            xtcl            = 100, # tick length in %
            ytcl1           = 100, # tick length in %
            ytcl2           = 100, # tick length in %
            stats.line      = 100 # line height in %
          )
        )
      )
    } else {
      stop("Layout definition not supported!")
    }
  } else if(is.list(layout) == TRUE) {

    ## user-specific layout definition assignment
    layout <- layout
  }

  ## return layout parameters
  return(layout)
}
