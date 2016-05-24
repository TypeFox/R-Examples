#' Plot results of GlobalSens function
#' @description Plot results of of \code{\link{CalculateGlobalSens}} function.
#' @param global.out output from \code{\link{CalculateGlobalSens}} function.
#' @param legend.label string with name for the legend.
#' @param x.label string with the name of x axis.
#' @param y.label string with the name of y axis.
#' @param mm.label string with the name of the envelope calculated using the minimum and maximum ranges.
#' @param sd.label string with the name of the envelope calculated using the mean +- standard deviation ranges.
#' @details Font size of saved plots is usually different to the font size seen in graphic browsers. Before changing font sizes, see the final result in saved (or preview) plots.
#'  
#' Other details of the plot can be modifyed using appropriate functions from \code{ggplot2} package.
#' @references \url{http://oswaldosantos.github.io/capm}
#' @seealso \link[deSolve]{plot.deSolve}.
#' @export
#' @examples 
#' #####################
#' ## SolveIASA model ##
#' #####################
#' 
#' ## Parameters and intial conditions.
#' pars.solve.iasa = c(
#'    b1 = 21871, b2 = 4374,
#'    df1 = 0.104, dm1 = 0.098, df2 = 0.125, dm2 = 0.118,
#'    sf1 = 0.069, sf2 = 0.05, sm1 = 0.028, sm2 = 0.05,
#'    k1 = 98050, k2 = 8055, h1 = 1, h2 = 0.5,
#'    a = 0.054, alpha = 0.1, v = 0.2, z = 0.1)
#'    
#' init.solve.iasa = c(
#'    f1 = 33425, fs1 = 10865,
#'    m1 = 38039, ms1 = 6808,
#'    f2 = 3343, fs2 = 109,
#'    m2 = 3804, ms2 = 68)
#'    
#' 
#' # Solve for point estimates.
#' solve.iasa.pt <- SolveIASA(pars = pars.solve.iasa, 
#'                           init = init.solve.iasa, 
#'                           time = 0:15, method = 'rk4')
#' 
#' ## Set ranges 10 % greater and lesser than the
#' ## point estimates.
#' rg.solve.iasa <- SetRanges(pars = pars.solve.iasa)
#'
#' ## Calculate golobal sensitivity of combined parameters.
#' ## To calculate global sensitivity to each parameter, set
#' ## all as FALSE.
#' glob.all.solve.iasa <- CalculateGlobalSens(
#'   model.out = solve.iasa.pt,
#'   ranges = rg.solve.iasa, 
#'   sensv = 'n2', all = TRUE)
#' 
#' ### Plot the sensitivities of combined parameters.
#' PlotGlobalSens(glob.all.solve.iasa)
#'
PlotGlobalSens <- function(global.out = NULL, x.label = 'Time', y.label = 'Population', legend.label = 'Sensitivity range', mm.label = 'min - max', sd.label = 'mean +- sd   ') {
  x <- Mean <- Min <- Max <- Sd <- NULL
  if (colnames(global.out)[length(global.out)] == 'param') {
    ggplot(global.out, aes(x = x, y = Mean)) +
      geom_ribbon(aes(ymin = Min, ymax = Max, fill = 'red'), 
                  alpha = .6) +
      geom_ribbon(aes(ymin = Mean - Sd, ymax = Mean + Sd, fill = 'blue'), 
                  alpha = .6) +
      geom_line() + facet_wrap( ~ param) +
      xlab(x.label) + ylab(y.label) +
      scale_fill_manual(
        name = legend.label,
        values = c('blue', 'red'),
        labels = c(sd.label, mm.label)) +
      theme(legend.position = 'top') +
      guides(fill = guide_legend(
        override.aes = list(alpha = 0.5)))
  } else {
    ggplot(global.out, aes(x = x, y = Mean)) +
      geom_ribbon(aes(ymin = Min, ymax = Max, fill = 'red'), 
                  alpha = .6) +
      geom_ribbon(aes(ymin = Mean - Sd, ymax = Mean + Sd, fill = 'blue'), 
                  alpha = .6) +
      geom_line() +
      xlab(x.label) + ylab(y.label) +
      ylim(0, max(global.out$Max)) +
      scale_fill_manual(
        name = legend.label,
        values = c('blue', 'red'),
        labels = c(sd.label, mm.label)) +
      theme(legend.position = 'top') +
      guides(fill = guide_legend(
        override.aes = list(alpha = 0.4)))    
  }
}
