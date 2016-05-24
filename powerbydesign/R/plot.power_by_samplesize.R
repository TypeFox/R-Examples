#' Plot Power Results
#'
#' Plots the power-by-samplesize data.frames returned by the power functions.
#'
#' @param x list, results of a power function
#' @param crit_power numeric, critical power value one is looking for (adds a line
#'                            and the critical sample size to the plot)
#' @param plot_dir character, name of an existing directory where the power plots
#'                            should be saved to (leaving it  NULL will plot to
#'                            the default device)
#' @param ... additional parameters, not used at the moment
#' @seealso \code{\link{boot.power.anova}}
#'
#' @export
plot.power_by_samplesize <- function(
  x,
  crit_power = NULL,
  plot_dir = NULL,
  ...
) {

  if (! class(x) == "power_by_samplesize") stop("x must be an object of class power_by_samplesize")

  if (! is.null(plot_dir) && ! dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }

  for (effect_name in names(x)) {
    if (! is.null(plot_dir)) {grDevices::png(paste0(plot_dir,"/",stringr::str_replace_all(effect_name,":","."),".png"),400,400)}

    dat <- x[[effect_name]]

    n_total <- power <- NULL # for R CMD check

    p <- ggplot2::ggplot(dat) +
      ggplot2::theme_gray(base_size=15) +
      ggplot2::ggtitle(effect_name) +
      ggplot2::ylab("Power") +
      ggplot2::xlab("Total Number of Participants") +
      ggplot2::ylim(0,1) +
      ggplot2::geom_line(ggplot2::aes(x=n_total,y=power),size=2) +
      ggplot2::geom_point(ggplot2::aes(x=n_total,y=power),size=5)
    if (! is.null(crit_power)) {
      p <- p +
        ggplot2::geom_hline(yintercept = crit_power, linetype="longdash")
      crit_n <- dat$n_total[dat$power >= crit_power]
      if (length(crit_n) > 0) {
        crit_n <- min(crit_n)
        p <- p +
          ggplot2::geom_vline(xintercept = crit_n, linetype="dotdash") +
          ggplot2::geom_text(ggplot2::aes(x=crit_n, y=0.5, label=paste0("\nCritical N: ",crit_n)), colour="blue", angle=90, size=6)
      }
    }

    print(p) # print required for ggplot in for loop

    if (! is.null(plot_dir)) {grDevices::dev.off()}
  }
}
