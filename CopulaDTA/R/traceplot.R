#' Trace plot using ggplot2.
#'
#' @param x A cdtafit object from \link{fit}.
#' @param ... additional options. See \link[rstan]{stan_trace} for more details.
#' @return A ggplot trace plot of the parameters of the models mean structure.
#' @examples
#' \dontrun{
#' fit1 <- fit(model1,
#'                 SID='ID',
#'                 data=telomerase,
#'                 iter=2000,
#'                 warmup=1000,
#'                 thin=1,
#'                 seed=3)
#'
#' traceplot(fit1)
#' }
#'@export
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}

traceplot.cdtafit <- function(x, ...){
    #open new window
    g <- rstan::stan_trace(x@fit, pars=c('MUse', 'MUsp'), ...)
    if (grDevices::dev.interactive()) grDevices::dev.new()
    print(g)	}
