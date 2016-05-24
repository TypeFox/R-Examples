#' @title Plot STRUCTURE Results
#' @description Plot Q-matrix from a call to \code{\link{structure}} or 
#'   \code{\link{clumpp}}.
#' 
#' @param q.mat matrix or data.frame of assignment probabilities.
#' @param pop.col column number identifying original population designations.
#' @param prob.col column number of first assignment probabilities to first 
#'  group. It is assumed that the remainder of columns 
#'  (\code{prob.col:ncol(q.mat)}) contain all assignment probabilities.
#' @param sort.probs logical. Sort individuals by probabilities within 
#'   populations? If \code{FALSE} individuals will be plotted as in \code{q.mat}.
#' @param label.pops logical. Label the populations on the plot?
#' @param col colors to use for each group.
#' @param horiz logical. Plot bars horizontally.
#' @param legend.position the position of the legend (\code{"top", "left", 
#'   "right", "bottom"}, or two-element numeric vector).
#' 
#' @return invisibly, the ggplot object
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{structure}}, \code{\link{clumpp}}
#'
#' @importFrom ggplot2 ggplot aes_string geom_area ylab theme geom_vline 
#'   scale_x_continuous xlab coord_flip scale_fill_manual element_blank
#' @importFrom reshape2 melt
#' @importFrom RColorBrewer brewer.pal
#' @export
#' 
structurePlot <- function(q.mat, pop.col = 3, prob.col = 4, sort.probs = TRUE,
                          label.pops = TRUE, col = NULL, horiz = TRUE,
                          legend.position = c("top", "left", "right", "bottom", "none")) {
  
  legend.position <- match.arg(legend.position)
  
  # convert q.mat to sorted data.frame
  prob.cols <- prob.col:ncol(q.mat)
  qm <- data.frame(q.mat)
  qm[, pop.col] <- factor(
    qm[, pop.col], 
    levels = sort(unique(qm[, pop.col]), decreasing = horiz)
  )
  sort.cols <- c(pop.col, if(sort.probs) rev(prob.cols) else NULL)
  i <- do.call(order, qm[, sort.cols, drop = FALSE])
  qm <- qm[i, ]
  qm$x <- 1:nrow(qm)
  
  # Get population frequencies, centers and dividing points
  pop.freq <- table(qm[, pop.col])
  levels(qm[, pop.col]) <- paste(
    levels(qm[, pop.col]), "\n(n = ", pop.freq, ")", sep = ""
  )
  pop.cntr <- tapply(qm$x, qm[, pop.col], mean)
  pop.div <- tapply(qm$x, qm[, pop.col], min)[-1] - 0.5
  
  # Create data.frame for plotting
  df.cols <- colnames(qm)[c(pop.col, prob.cols)]
  df.cols <- c("x", df.cols)
  df <- melt(qm[, df.cols], id.vars = c(1, 2),
             variable.name = "Group", value.name = "probability")
  colnames(df)[1:2] <- c("x", "population")
  df <- df[order(-as.numeric(df$Group), df$probability), ]
  
  # Plot stacked bar graphs
  g <- ggplot(df, aes_string("x", "probability")) + 
    geom_area(aes_string(fill = "Group"), stat = "identity") +
    ylab("Pr(Group Membership)") +
    theme(
      axis.ticks.x = element_blank(),
      legend.position = legend.position,
      legend.title = element_blank()
    )
  if(label.pops) {
    g <- g + geom_vline(xintercept = pop.div) +
      scale_x_continuous(name = "", breaks = pop.cntr, labels = names(pop.cntr))
  } else {
    g <- g + xlab("") + theme(axis.text.x = element_blank())
  }
  if(horiz) g <- g + coord_flip()
  if(!is.null(col)) g <- g + scale_fill_manual(values = col)
  
  print(g)
  invisible(g)
}