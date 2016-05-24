#'  Structure Heatmap with Missing Value Demonstration
#' @details Structure heatmap is like a normal heatmap, but is particulary useful
#' when they are missing values in the data matrix. Default color were carefully
#' chosen so normally it is a good choice for your data. However, you are still 
#' encouraged to play around with it.
#' @param data a data frame or matrix, possibly with missing values denoted by NA
#' @param high.col color for high values, can be a number or a color name, default
#' is "steelblue.
#' @param low.col color for high values, can be a number or a color name, default
#' is "white".
#' @param missing.col color for missing values, can be a number or a color nam,
#' default is "gold"
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis .
#' 
#' @export
#' @examples
#' data(hqmr.data)
#' # use a subset of the hqmr.data
#' # notice the gold chunks which represent missing values
#' subdata <- hqmr.data[1000:1200, 1:30]
#' HeatStruct(subdata)
#' # change colors for high.col, low.col and missing.col
#' HeatStruct(subdata, low.col = "blue", high.col = "red", missing.col = "black")
HeatStruct <- function (data, high.col = "steelblue", low.col = "white", 
                        missing.col = "gold", xlab = "", ylab = "" ) {
  rescale <- function (x) {
    rng <- range(x, na.rm = TRUE)
    (x - rng[1])/(rng[2] - rng[1])
  }
  cu <- function (a, b) {
    if (length(a) == 0) 
      return(b)
    if (length(b) == 0) 
      return(a)
    cbind(a, b[setdiff(names(b), names(a))])
  }
  ggp <- function (data, vars = names(data), ...) {
    scaled <- as.data.frame(lapply(data[, vars], rescale))
    data <- cu(scaled, data)
    data$ROWID <- 1:nrow(data)
    molten <- melt(data, m = vars)
    ggplot(molten, aes_string(x = "variable", y = "value", group = "ROWID"), 
           ...)
  }
  if(!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  p <- ggp(data) +
    aes_string(y = "ROWID", fill = "value", x = "variable") +
    geom_tile() +
    scale_y_continuous(expand = c(0, 1)) +
    scale_fill_continuous(low = low.col, high = high.col, 
                          na.value = missing.col, guide = "colorbar") +
    labs(x = xlab, y = ylab)
  p
}
