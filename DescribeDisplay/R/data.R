#' Load describe display
#' Retrieve output of from describe display plugin
#'
#' Also performs some conversion of data structures to more
#' conveient form so that other functions do not have to repeatedly
#' recompute.  Some of these conversions could probably be moved into
#' the Describe Display plugin, but it may be easier to just do them
#' on the R side..
#'
#' @param path file path
#' @return object of class dd
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @seealso \code{\link{dd_example}} for an easier way of loading example
#'   files
#' @keywords  manip
#' @export
dd_load <- function(path) {
  dd <- source(path)$value
  class(dd) <- c(dd_plot_class(dd$type), "dd")
  dd$colormap$foreground <- sapply(dd$colormap$foregroundColors,
    function(x) do.call(rgb, as.list(x))
  )
  dd$colormap$foregroundColors <- NULL
  cols <- dd$ncols %||% 1
  dd$dim <- c(dd$nplots / cols, cols)
  dd$plots <- lapply(1:dd$nplots, function(n) dd_clean_plot(dd, n))

  dd
}

#' Load example describe display file
#' Load example describe display file included with package.
#'
#' These are mainly used for testing.
#'
#' @param name name of example
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords internal
#' @export
#' @examples
#' a <- dd_example("xyplot")
dd_example <- function(name) {
  file <- paste(name, ".R", sep = "")
  path <- system.file("examples", file, package = "DescribeDisplay")

  if (!file.exists(path)) {
    stop("Cannot find example ", name, call. = FALSE)
  }

  dd_load(path)
}

#' Clean plot data structure
#' Cleans up plot data structure into consistent, easy to use data structure
#'
#' @param dd dd object
#' @param n plot number
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords internal
#' @importFrom scales expand_range
dd_clean_plot <- function(dd, n = 1) {
  names(dd$plots[[n]]) <- gsub("plot", "", names(dd$plots[[n]]))
  plot <- c(
    list(
      points = dd_points(dd, n),
      edges = dd_edges(dd, n)
    ),
    dd$plots[[n]][c("type", "projection", "params")]
  )

  if (!is.null(plot$projection)) {
    plot$baseline <- if (plot$projection == "1D plot") {
       0
    } else {
      (min(plot$points$y) - 0.05 * abs(min(plot$points$y)))
    }
  }

  if (identical(dd$plots[[n]]$scale, c(0.7, 0.7))) {
    plot$xscale <- expand_range(range(plot$points$x), 0.1)
    plot$yscale <- expand_range(range(plot$points$y), 0.1)
  } else if (
    (!is.null(dd$plots[[n]]$tformLims)) &
    (sum(dd$plots[[n]]$tformLims[1:2]) == 0)
  ) {
    plot$xscale <- range(dd$plots[[n]]$planarLims[1:2])
    plot$yscale <- range(dd$plots[[n]]$planarLims[3:4])

    if (diff(plot$yscale) == 0 ) {
      plot$yscale <- expand_range(range(plot$points$y), 0.1)
      print("THIRD")
    }

  } else {
    plot$xscale <- dd$plots[[n]]$tformLims[1:2]
    plot$yscale <- dd$plots[[n]]$tformLims[3:4]
  }

  if (!is.null(dd$plots[[n]]$stickylabels)) {
    labels <- do.call(
      rbind,
      lapply(dd$plots[[n]]$stickylabels, as.data.frame)
    )
    labels <- cbind(
      plot$points[labels$index + 1, c("x", "y")],
      label = labels$label
    )
    rl <- (labels$x - plot$xscale[1]) / diff(plot$xscale) < 0.5
    # tb <- (labels$y - plot$yscale[1]) / diff(plot$yscale) < 0.5
    labels$left <- ifelse(rl, 0, 1)
    labels$top <-  1 #ifelse(tb, 0, 1)

    labels$x <- labels$x + (-1 + 2 * rl) * 0.01 * diff(plot$xscale)
    #labels$y <- labels$y + (-1 + 2 * tb) * 0.01 * diff(plot$yscale)
    plot$labels <- labels
  }

  class(plot) <- c(plot$type, dd_plot_class(plot$projection), "ddplot")

  plot
}

#' Describe display points data
#' Retrieves the describe display points data for the given plot number.
#'
#' @param dd list of values from describe display
#' @param n plot number, defaults to first plot
#' @return data frame suitable for plotting
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords internal
dd_points <- function(dd, n = 1) {
  df <- as.data.frame(dd$plots[[n]]$points)
  df$hidden <- df$hidden != 0
  cmap <- dd$colormap$foreground
  hiddencolour <- do.call(rgb, as.list(dd$colormap$hiddenColor))
  # Remap point aesthetics to R appropriate values
  df$col <- factor(
    ifelse(df$hidden, hiddencolour, cmap[df$color + 1]),
    levels = c(rev(cmap), hiddencolour)
  )

  if (is.null(df$glyphtype)) {
    df$pch <- df$cex <- rep(1, nrow(df))
  } else {
    df$pch <- c(18, 3, 4, 1, 0, 16, 15)[df$glyphtype + 1]
    df$cex <- (df$glyphsize + 1) / 6
  }

  rownames(df) <- df$index %||% seq_len(nrow(df))

  dfNames <- intersect(names(df), c("x", "y", "col", "pch", "cex", "hidden"))
  df[order(!df$hidden), dfNames]
}

#' Describe display edge data
#' Retrieves the describe display edge data for the given plot number.
#'
#' @param dd list of values from describe display
#' @param n plot number, defaults to first plot
#' @return data frame suitable for plotting
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords internal
dd_edges <- function(dd, n = 1) {
  if (is.null(dd$plots[[n]]$edges)) return()
  df <- do.call(rbind, lapply(dd$plots[[n]]$edges, as.data.frame))

  # Remap edge aesthetics to appropriate values
  df$col <- dd$colormap$foreground[df$color + 1]
  df$lwd <- (df$lwd + 1) / 2
  df$lty <- rep(1, 6)[df$ltype + 1]

  # Return only visible edges
  df <- df[!df$hidden, c("src", "dest", "col", "lwd", "lty")]

  points <- dd_points(dd, n)
  src <- points[as.character(df$src), c("x", "y")]
  names(src) <- c("src.x", "src.y")
  dest <- points[as.character(df$dest), c("x", "y")]
  names(dest) <- c("dest.x", "dest.y")

  cbind(src, dest, df)
}

#' Describe display plot class
#' Compute valid R class name for given plot type
#'
#' @param projection type of projection that should be used
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords internal
dd_plot_class <- function(projection) {
  gsub("\\s+", "", tolower(projection))
}

#' Describe display plot defaults
#' Gather overall plot defaults for specified plot
#'
#' @param dd list of values from describe display
#' @param n plot number, defaults to first plot
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords internal
dd_defaults <- function(dd, n = 1) {
  list(
    main = dd$title,
    xlab = dd$plots[[n]]$plotparams$xlab %||% "",
    ylab = dd$plots[[n]]$plotparams$ylab %||% "",
    axes = FALSE
  )
}

#' Describe display tour axis
#' Return representation of axes for specified plot
#'
#' @param plot list of information of a plot
#' @author Hadley Wickham h.wickham [at] gmail.com
#' @keywords internal
dd_tour_axes <- function(plot) {
  if (is.null(plot$params[["F"]])) return()


  if (plot$projection == "1D Tour") {
    proj <- matrix(plot$params[["F"]], ncol = 1)
    colnames(proj) <- "x"
  } else {
    proj <- matrix(plot$params[["F"]], ncol = 2, byrow = FALSE)
    colnames(proj) <- c("x", "y")
  }

  lbls <- plot$params$labels

  ranges <- do.call(rbind,  plot$params$ranges)
  df <- data.frame(proj, label = lbls, range = ranges)

  if (plot$projection == "2D Tour") {
    df$r <- with(df, sqrt(x^2 + y^2)) # nolint
    df$theta <- atan2(df$y, df$x)
  } else {
    df <- df[nrow(df):1, ]
  }

  df
}

#' @export
print.dd <- function(x, ...) str(x)
