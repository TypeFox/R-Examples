#' Plot sample time series from a portfolio simulation
#'
#' @param X  Object to plot. Should be a list of outputs from
#' \code{\link{meta_sim}}.
#' @param ylim Y axis limits.
#' @param x_axis Should an x axis be added?
#' @param y_axis  Should a y axis be added?
#' @param rate If \code{TRUE} then the first difference (rate of change) will be
#' plotted. If \code{FALSE} then the raw data will be plotted.
#' @param lwd Line width of the lines.
#' @param y_axis_ticks Location of the y-axis tick marks, if you want to specify
#' them.
#' @param start_new_plots On which elements of the list \code{X} should new
#' panels be started? A numeric vector.
#' @param labels Labels for the panels.
#' @param burn Burn in period to discard.
#' @param add_lm Add a regression trend line?
#' @param cols Colours for the lines. A vector of character.
#' @param ... Anything else to pass to \code{\link[graphics]{plot.default}}
#' @export
#' @return
#' A plot, possibly with multiple panels.
#' @examples
#' w_plans <- list()
#' w_plans[[1]] <- c(5, 1000, 5, 1000, 5, 5, 1000, 5, 1000, 5)
#' w_plans[[2]] <- c(5, 5, 5, 1000, 1000, 1000, 1000, 5, 5, 5)
#' w_plans[[3]] <- c(rep(1000, 4), rep(5, 6))
#' w_plans[[4]] <- rev(w_plans[[3]])
#' w <- list()
#' for(i in 1:4) { # loop over plans
#'  w[[i]] <- list()
#'  for(j in 1:2) { # loop over trials
#'    w[[i]][[j]] <- matrix(w_plans[[i]], nrow = 1)
#'  }
#' }
#'
#' cons_arma_ts <- list()
#' arma_env_params <- list(mean_value = 16, ar = 0.1, sigma_env = 2, ma = 0)
#' for(i in 1:4) {
#'   use_cache <- ifelse(i == 1, FALSE, TRUE)
#'   cons_arma_ts[[i]] <- meta_sim(b = w[[i]][[1]], n_pop = 10, env_params =
#'     arma_env_params, env_type = "arma", assess_freq = 5,
#'     use_cache = use_cache)
#' }
#' cols <- RColorBrewer::brewer.pal(5, "Dark2")
#' par(mfrow = c(2, 1))
#' plot_sp_A_ts(cons_arma_ts, ylim = c(0000, 12400),
#'   start_new_plots = c(1, 3),
#'   labels = c("Balanced response diversity",
#'     "ignore", "Unbalanced response diversity", "ignore"), cols = cols)

plot_sp_A_ts <- function(X, ylim, x_axis = TRUE, y_axis = TRUE, rate = FALSE,
  lwd = 1.7, y_axis_ticks = NULL, start_new_plots = 1, labels = NULL, burn = 30,
  add_lm = FALSE, cols, ...) {
  #A_range <- ldply(X, function(x) range(rowSums(x$A[-burn, ])))
  burn <- 1:burn
  for(i in 1:length(X)){
    if(i %in% start_new_plots) {
      plot(1,1,ylim = ylim, xlim = c(1, 70), type = "n",
        xlab = "", ylab = "", xaxt = "n", axes = FALSE, yaxs = "i", ...)
      if(!is.null(labels)) {
        x1 <- par("usr")[1]
        x2 <- par("usr")[2]
        y1 <- par("usr")[3]
        y2 <- par("usr")[4]
        text(x1 + (x2-x1)*0.005, y2 - (y2-y1)*0.17, labels = labels[i], pos = 4)
      }
      if(y_axis) {
        if(is.null(y_axis_ticks)) {
          axis(2, col=  "grey50", at = pretty(axTicks(2), n = 2))
        } else {
          axis(2, col=  "grey50", at = y_axis_ticks)
        }
      }
    }

    x <- X[[i]]$A[-burn, ]
    port.x <- rowSums(x)
    if(!rate) {
    lines(1:length(port.x), port.x, col = cols[i], lwd = lwd, lty = 1)
    } else {
    lines(2:length(port.x), diff(log(port.x)), col = cols[i], lwd = lwd, lty = 1)
    if(add_lm) {
      m <- lm(diff(log(port.x)) ~ c(1:length(diff(log(port.x)))))
      abline(m, col = cols[i], lty = 2, lwd = lwd)
    }
    }
    box(col = "grey50")
  }
  if(x_axis)
    axis(1, col=  "grey50")
}
