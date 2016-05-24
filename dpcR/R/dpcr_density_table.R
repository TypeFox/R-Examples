#' Calculate Density of Multiple dPCR runs
#' 
#' Calculates the density of the number of positive molecules or the average number 
#' of molecules per partition of \code{\link{dpcr}} objects.
#' 
#' @param input an object of class \code{\linkS4class{dpcr}}.
#' @inheritParams dpcr_density
#' @return A list (with the length equal to the number of runs in \code{input}) of 
#' data frames containing densities and borders of confidence intervals.
#' @author Michal Burdukiewicz, Stefan Roediger.
#' @seealso \code{\link{dpcr_density}} for easy analysis and plots of single runs.
#' @export
#' @keywords dplot hplot
#' @examples
#' dens <- dpcr_density_table(six_panels)
#' 
#' # create plot using ggplot2
#' library(ggplot2)
#' 
#' ggplot(dens[["Experiment2.2"]], aes(x = x, y = y)) + 
#'   geom_line() + 
#'   geom_area(aes(fill = !(conf_up | conf_low))) +
#'   scale_y_continuous("Density") +
#'   scale_fill_discrete("0.95 CI")
#' 
dpcr_density_table <- function(input, average = FALSE, methods = "wilson", conf.level = 0.95) {
  res <- lapply(1L:ncol(input), function(run_id) {
    single_run <- extract_dpcr(input, run_id)
    kn <- unlist(summary(single_run, print = FALSE)[["summary"]][1, c("k", "n")])
    
    conf <- dpcr_density(k = kn["k"], n = kn["n"], average = average, 
                         methods = methods, 
                         conf.level = conf.level, plot = FALSE)
    
    dens <- data.frame(dpcr_calculator(kn["k"], kn["n"], average = average))
    colnames(dens) <- c("x", "y")
    
    dens[["conf_low"]] <- dens[["x"]] <= conf[["lower"]] 
    dens[["conf_up"]] <- dens[["x"]] >= conf[["upper"]]
    dens
  })
  names(res) <- colnames(input)
  res
}