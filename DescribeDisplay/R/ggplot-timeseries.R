#' Compact time series data
#' A time series plot is made from scatterplots with a common x axis.
#' This function pulls the correct information out of the data.
#'
#' @param data data to display
#' @author Barret Schloerke \email{schloerke@@gmail.com}
#' @keywords internal
#' @importFrom reshape2 dcast
compact_timeseries <- function(data){
  dfx <- data.frame(
    data$plots[[1]]$points[, c("col", "pch", "cex")],
    value = data$plots[[1]]$points$x,
    variable = data$plots[[1]]$params$xlabel,
    id = 1:nrow(data$plots[[1]]$points)
  )

  df <- do.call(rbind, lapply(data$plots, function(p) {
      aes <- p$points[, c("col", "pch", "cex")]
      data.frame(
        aes,
      value = p$points$y,
        variable = p$params$ylabel,
        id = 1:nrow(p$points)
      )
    }))

  df <- dcast(df, id + ... ~ variable)
  dfx <- dcast(dfx, id + ... ~ variable)

  cPCI <- c("cex", "pch", "col", "id")
  namesDf <- names(df)
  namesInDf <- namesDf %in% cPCI
  df <- cbind(
    df[, namesInDf],
     dfx[, setdiff(names(dfx), cPCI)],
    df[, setdiff(namesDf, cPCI)]
  )
  colnames(df)[sum(namesInDf) + 1] <- data$plots[[1]]$params$xlabel

  return(df)

}


#' Create nice plots for a time series
#' Create nice looking plots complete with axes using ggplot.  Produces graphics with a uniform x axis.
#'
#' @param data to display
#' @param edges Boolean operator to tell whether to try to force the edges or not.  Will not work to remove the edges.
#' @param ... (currently) unused arguments
#' @author Barret Schloerke \email{schloerke@@gmail.com}
#' @keywords hplot
#' @export
#' @examples
#' library(ggplot2)
#' print(ggplot(dd_example("timeseries")))
#' print(ggplot(dd_example("timeseries"), edges = TRUE))
ggplot.timeseries <- function(data, edges = FALSE, ...){
  #cat("\nggplot.timeseries\n")

  x <- y <- NULL

  df <- compact_timeseries(data)

  data.par <- df[, colnames(df) %in% c("cex", "pch", "col", "id")]

  df <- df[, setdiff(colnames(df), colnames(data.par))]

  ## time series, one column no 1d plots
  grid <- expand.grid(x = 1, y = 1:ncol(df))

    grid <- subset(grid, x != y)

  all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
    xcol <- grid[i, "x"]
    ycol <- grid[i, "y"]

    data.frame(
      xvar = names(df)[ycol],
      yvar = names(df)[xcol],
      x = df[, xcol], y = df[, ycol], df, data.par
    )
  }))



  all$xvar <- factor(all$xvar, levels = names(df))
  all$yvar <- factor(all$yvar, levels = names(df))

    aesString <- aes_string(x = "x", y = "y", group = "xvar")
  class(aesString) <- "uneval"

  p <- ggplot(all, aesString) + facet_grid(xvar ~ yvar, scales = "free") +
    scale_colour_identity() +
      scale_size_identity() +
      scale_shape_identity() +
      scale_linetype_identity() +
      theme(title = element_text(data$title)) +
      scale_x_continuous(all[1, "yvar"]) +
      scale_y_continuous("") +
      geom_point(
        data = all,
        aes_string(size = "cex * 4", colour = "col", shape = "pch")
      )

  if (data$showDirectedEdges | data$showUndirectedEdges | edges == TRUE)
    p <- p +
      geom_path(
        data = all,
        aes_string(x = "x", y = "y", size = "cex", colour = "col")
      )

  p
}
