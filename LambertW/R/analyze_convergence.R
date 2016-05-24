#' @title Analyze convergence of Lambert W estimators
#' @export
#'
#' @description Analyzes the feasibility of a Lambert W x F distribution for a
#'     given dataset based on bootstrapping.  In particular it checks whether
#'     parameter estimates support the hypothesis that the data indeed follows a
#'     Lambert W x F distribution with finite mean and variance of the input
#'     distribution, which is an implicit assumption of Lambert W x F random
#'     variables in Goerg (2011).
#' 
#' See Goerg (2016) for an alternative definition that does not rely on fnite
#' second order moments (set \code{use.mean.variance = FALSE} to use that type
#' of Lambert W \eqn{\times} F distributions).
#'
#' @details
#' Stehlik and Hermann (2015) show that when researchers use the IGMM algorithm
#'     outlined in Goerg (2011) erroneously on data that does not have finite
#'     input variance (and hence mean), the algorithm estimates do not converge.
#'
#' In practice, researchers should of course first check if a given model is
#' appropriate for their data-generating process.  Since original Lambert W x F
#' distributions assume that mean and variance are finite, it is not a given
#' that for a specific dataset the Lambert W x F setting makes sense.
#'
#' The bootstrap analysis reverses Stehlik and Hermann's argument and checks
#'     whether the IGMM estimates \eqn{\lbrace \hat{\tau}^{(n)} \rbrace_{n}}
#'     converge for increasing (bootstrapped) sample size \eqn{n}: if they do,
#'     then modeling the data with a Lambert W x F distribution is appropriate;
#'     if estimates do not converge, then this indicates that the input data is
#'     too heavy tailed for a classic skewed location-scale Lambert W x F
#'     framework. In this case, take a look at (double-)heavy tailed Lambert W x
#'     F distributions (\code{type = 'hh'}) or unrestricted location-scale
#'     Lambert W x F distributions (\code{use.mean.variance = FALSE}). For
#'     details see Goerg (2016).
#'
#' @references
#' Stehlik and Hermann (2015). ``Letter to the Editor''. Ann. Appl. Stat. 9
#'    2051. doi:10.1214/15-AOAS864 -- \url{http://projecteuclid.org/euclid.aoas/1453994190}
#' 
#' @param LambertW_fit,object,x an object of class \code{"LambertW_fit"} with an
#'     \code{IGMM} or \code{MLE_LambertW} estimate.
#' @param sample.sizes sample sizes for several steps of the convergence
#'     analysis.  By default, one of them equals the length of the original
#'     data, which leads to improved plots (see
#'     \code{\link{plot.convergence_LambertW_fit}}); it is not necessary,
#'     though.
#' @param ... additional arguments passed to \code{\link{bootstrap}} or
#'     \code{\link[boot]{boot.ci}} in \pkg{boot} package.
#' @examples
#' \dontrun{
#' 
#' sim.data <- list("Lambert W x Gaussian" = 
#'                     rLambertW(n = 100, distname = "normal", 
#'                               theta = list(gamma = 0.1, beta = c(1, 2))),
#'                  "Cauchy" = rcauchy(n = 100))
#' # do not use lapply() as it does not work well with match.call() in
#' # bootstrap()
#' igmm.ests <- list()
#' conv.analyses <- list()
#' for (nn in names(sim.data)) {
#'   igmm.ests[[nn]] <- IGMM(sim.data[[nn]], type = "s")
#'   conv.analyses[[nn]] <- analyze_convergence(igmm.ests[[nn]])
#' }
#' plot.lists <- lapply(conv.analyses, plot)
#' for (nn in names(plot.lists)) {
#'   plot.lists[[nn]] <- lapply(plot.lists[[nn]], "+", ggtitle(nn))
#' }
#' 
#' require(gridExtra)
#' for (jj in seq_along(plot.lists[[1]])) {
#'   grid.arrange(plot.lists[[1]][[jj]], plot.lists[[2]][[jj]], ncol = 2)
#' }
#' }
#' 


analyze_convergence <- function(LambertW_fit, 
                                sample.sizes = 
                                  round(seq(0.2, 1, length = 5) * 
                                          length(LambertW_fit$data)),
                                ...) {
  stopifnot(inherits(LambertW_fit, "LambertW_fit"),
            is.numeric(sample.sizes),
            sample.sizes > 0,
            !any(duplicated(sample.sizes)))
  
  sample.sizes <- sort(sample.sizes, decreasing = FALSE)
  
  boots <- list()
  for (ss in sample.sizes) {
    boots[[as.character(ss)]] <- bootstrap(LambertW_fit, sample.size = ss,
                                           ...)
  }
  
  out <- list(boots = boots,
              original.sample.size = length(LambertW_fit$data))
  class(out) <- c(class(out), "convergence_LambertW_fit")
  return(out)
}

#' @rdname analyze_convergence
#' @param type type of confidence interval from bootstrap estimates. Passes this
#' argument along to \code{\link[boot]{boot.ci}}.  However, contrary to 
#' the \code{type} argument in \code{\link[boot]{boot.ci}}, the \code{summary}
#' function can only take one of \code{c("basic", "norm", "perc", "bca")}.
#' See \code{\link[boot]{boot.ci}} for details.
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#' @export

summary.convergence_LambertW_fit <- function(object, 
                                             type = c("basic", "norm", "perc", "bca"),
                                             ...) {
  stopifnot(inherits(object, "convergence_LambertW_fit"),
            is.character(type))
  
  type <- match.arg(type)
  
  num.params <- length(object$boots[[1]]$t0)
  param.names <- names(object$boots[[1]]$t0)

  dot.args <- list(...)
  dot.args[["type"]] <- type
  
  requireNamespace("boot")
  .extractCI <- function(boot.out, index) {
    bootci.out <- boot::boot.ci(boot.out, index = index, type = type, ...)
    # select the entry of the list with the name of 'type', which contains
    # the estimates of confidence intervals returned by boot.ci
    mat.intervals <- tail(bootci.out, 1)[[1]]
    # last two entries of that row vector are lower and upper cis
    lower.upper <- tail(mat.intervals[1, ], 2)
    return(lower.upper)
  }
  
  boots <- object$boots
  cis <- lapply(boots, function(x) {
                out <- sapply(seq_len(num.params), 
                              function(ii) {
                                .extractCI(boot.out = x, index =  ii)
                               })
                out <- rbind(out,
                             colMeans(x$t), apply(x$t, 2, sd))
                colnames(out) <- param.names
                rownames(out) <- c("lower", "upper", "mean", "sd")
                return(out)
                })
  cis <- lapply(cis, reshape2::melt, varnames = c("ci", "parameter"))
  for (ss in names(cis)) {
    cis[[ss]][["sample.size"]] <- as.numeric(ss)
  }
  cis <- do.call("rbind", cis)
  cis <- reshape2::dcast("parameter  + sample.size ~ ci", data = cis, 
                         direction = "wide")
  
  boots.dt <- list()
  for (hh in names(boots)) {
    boots.dt[[hh]] <- data.frame(boots[[hh]]$t)
    boots.dt[[hh]] <- cbind(boots.dt[[hh]], 
                            "sample.size" = as.numeric(hh))
  }
  boot.dt <- do.call("rbind", boots.dt)
  colnames(boot.dt) <- c(param.names, "sample.size")
  rownames(boot.dt) <- NULL
  

  # use sample size which equals total (ie, the last one in the list)
  orig.param.est <- 
    object$boots[[as.character(object$original.sample.size)]]$t0
  
  out <- list(boots = boot.dt,
              cis = cis,
              estimate = orig.param.est)
  class(out) <- c(class(out), "summary.convergence_LambertW_fit")
  return(out)
}

#' @rdname analyze_convergence
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export

plot.convergence_LambertW_fit <- function(x, ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is not available.  Please install it",
         "to use this function.") 
  }
  
  object <- summary(x)
  out <- list()
  
  # generate colours for parameters
  num.params <- length(unique(object$cis$parameter))
  param.cols <- RColorBrewer::brewer.pal(n = max(3, min(5, num.params)), 
                                         name = "Set1")
  if (num.params > 5) {
    param.cols <- grDevices::colorRampPalette(param.cols)(num.params)
  }
  param.cols <- head(param.cols, num.params)

  out[["cis"]] <- 
    with(object$cis, 
         ggplot(object$cis, 
                aes(fill = parameter, colour = parameter)) +
    geom_ribbon(aes(x = sample.size, ymin = lower, ymax = upper),
                alpha = 0.4, size = 1) +
    geom_point(aes(x = sample.size, y = lower), pch = 19, size = 4) + 
    geom_point(aes(x = sample.size, y = upper), pch = 19, size = 4)) + 
    ylab("")
  
  out[["sd"]] <- 
    with(object$cis, 
         ggplot(object$cis, 
                aes(x = sample.size, 
                    y = sd, 
                    fill = parameter, 
                    colour = parameter)))
  
  out[["sd.sqrt.n"]] <- 
    with(object$cis,
         ggplot(object$cis, 
           aes(x = sample.size, 
                      y = sd * sample.size^(0.5), 
                      fill = parameter, 
                      colour = parameter)))
  
  boots.long <- reshape2::melt(object$boots, variable.name = "parameter",
                               id.vars = "sample.size")

  # add same colors and legend on bottom
  for (vv in c("cis", "sd", "sd.sqrt.n")) {
    out[[vv]] <- out[[vv]] + 
      theme(legend.position = "bottom") +
      scale_fill_manual("Parameter", values = param.cols) +
      scale_colour_manual("Parameter", values = param.cols) +
      facet_grid(parameter~., scales = "free") +
      xlab("Sample size n")
  }
  # add points and lines
  for (vv in c("sd", "sd.sqrt.n")) {
    out[[vv]] <- out[[vv]] +
      geom_line(size = 2, alpha = 0.6) +
      geom_point(pch = 19, size = 4)
  }
  
  out[["sd.sqrt.n"]] <- 
    out[["sd.sqrt.n"]] + ylab(expression(hat(sigma) %.% sqrt(n)))
  
  out[["sd"]] <- out[["sd"]] + ylab(expression(hat(sigma)))
  
  boots.long$sample.size <- factor(boots.long$sample.size)
  # generate colours for sample sizes
  num.sample.sizes <- length(unique(boots.long$sample.size))
  
  sample.size.cols <- 
      RColorBrewer::brewer.pal(n = max(3, min(7, num.sample.sizes)),
                               name = "YlOrRd")
  if (num.sample.sizes > 7) {
    sample.size.cols <- 
        grDevices::colorRampPalette(sample.size.cols)(num.sample.sizes)
  }
  sample.size.cols <- head(sample.size.cols, num.sample.sizes)
  
  out[["density"]] <- 
    with(boots.long,
         ggplot(boots.long,
                aes(value, fill = sample.size, group = sample.size)) +
    geom_density(alpha = 0.4) +
    facet_wrap(~parameter, scales = "free") +
    theme(legend.position = "bottom") + 
    scale_fill_manual("Sample size", values = sample.size.cols) + 
    scale_colour_manual("", values = sample.size.cols))
  
  # add original parameter estimate if available
  if (!is.null(object$estimate)) {
    param.dt <- data.frame(parameter = names(object$estimate),
                           estimate = object$estimate)

    # add vertical line at overall estimate
    out[["density"]] <- out[["density"]] + 
      with(param.dt,
           geom_vline(data = param.dt,
                      aes(xintercept = estimate), colour = "black",
                      linetype = "dashed", size = 2, alpha = 0.75))
    out[["cis"]] <- out[["cis"]] +
      with(param.dt,
           geom_hline(data = param.dt,
                      aes(yintercept = estimate, colour = parameter),
                      linetype = "dashed", size = 2))
  }
  class(out) <- c(class(out), "plot.convergence_LambertW_fit")
  return(out)
}
