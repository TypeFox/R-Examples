#' Calculate Density of Single dPCR Run
#' 
#' Calculates and plots the density of the number of positive
#' molecules or the average number of molecules per partition. Can be used for
#' both array digital PCR and droplet digital PCR.
#' 
#' @param k Total number of positive molecules.
#' @param n Total number of partitions.
#' @param average If \code{TRUE}, calculates density of the average number of
#' molecules per partition. If \code{FALSE}, instead performs calculations for
#' the total number of positive molecules.
#' @param methods Method for calculating the confidence interval. 
#' Possible values are: \code{"wilson"}, \code{"agresti-coull"},
#' \code{"exact"}, \code{"prop.test"}, \code{"profile"}, \code{"lrt"},
#' \code{"asymptotic"}, \code{"bayes"}, \code{"cloglog"}, \code{"logit"},
#' \code{"probit"}. Default value is \code{"wilson"}. See Details.
#' @param conf.level The level of confidence to be used in the confidence
#' interval. Values from 0 to 1 and -1 to 0 are acceptable.
#' @param plot If \code{TRUE}, plots density plot.
#' @param bars plot on density plot bars for discrete values of lambda.
#' @param ... Additional arguments send to \code{plot} function.
#' @return A data frame with one row containing bounds of the confidence
#' intervals and a name of the method used to calculate them.
#' @author Michal Burdukiewicz, Stefan Roediger.
#' @seealso Computation of confidence intervals: \link[binom]{binom.confint}, 
#' 
#' The browser-based graphical user interface for this function: \link{dpcr_density_gui}.
#' @references Brown, Lawrence D., T. Tony Cai, and Anirban DasGupta.
#' \emph{Confidence Intervals for a Binomial Proportion and Asymptotic
#' Expansions.} The Annals of Statistics 30, no. 1 (February 2002): 160--201.
#' @keywords dplot hplot
#' @examples
#' 
#' # Calculate the average number of molecules per partition and show the area
#' # of the confidence interval (left plot) and the area within the 
#' # confidence interval
#' par(mfrow = c(1,2))
#' dpcr_density(k = 25, n = 55, average = TRUE, methods = "wilson", 
#' 	     conf.level = 0.95)
#' dpcr_density(k = 25, n = 55, average = TRUE, methods = "wilson", 
#' 	     conf.level = -0.95)
#' par(mfrow = c(1,1))
#' 
#' # By setting average to FALSE the total number of positive molecules is 
#' # calculated
#' dpcr_density(k = 25, n = 55, average = FALSE, methods = "wilson", 
#' 	     conf.level = 0.95)
#' 
#' @export dpcr_density
dpcr_density <- function(k, n, average = FALSE, methods = "wilson", 
                         conf.level = 0.95, plot = TRUE, 
                         bars = FALSE, ...) {
  dat <- dpcr_calculator(k, n, average)
  conf <- binom.confint(k, n, methods = methods, conf.level = conf.level)
  if (average) {
    conf[, c(4:6)] <- - log(1 -  conf[, c(4:6)])
    names(conf)[4] <- "lambda"
    xlab = "Molecules/partition"
    main = "Number of molecules per partition"
  } else {
    conf[, 4:6] <- conf[, 4:6] * n
    xlab = "Positive partitions"
    main = "Number of positive partitions"
  }
  names(conf)[2] <- "k"
  if (plot) {
    plot_distr(dat, ylab = "Density", xlab = xlab,
               main = main, bars = bars, ...)
    plot_conf_int(conf[1, 4:6], dat, "left", conf_int_col = adjustcolor("cyan4", 
                                                                        alpha.f = 0.15), 
                  conf_int_border = adjustcolor("cyan4", alpha.f = 0.15))
    plot_conf_int(conf[1, 4:6], dat, "right", conf_int_col = adjustcolor("cyan4", 
                                                                         alpha.f = 0.15), 
                  conf_int_border = adjustcolor("cyan4", alpha.f = 0.15))
  }
  conf
}

plot_conf_int <- function(conf_int, data, side, 
                          conf_int_col = adjustcolor("cyan4", alpha.f = 0.15), 
                          conf_int_border = adjustcolor("cyan4", alpha.f = 0.15), 
                          ...) {
  #get position of conf between values in data
  y_val <- y_val_conf(conf_int, data, side)
  if (!is.null(y_val)) {
    if (side == "left") {
      id1 <- 1
      id2 <- y_val[1]
      dat_plot <- rbind(data[id1:id2,], y_val[3:2])   
    }
    if (side == "right") {
      id1 <- y_val[1]
      id2 <- nrow(data)
      dat_plot <- rbind(y_val[3:2], data[id1:id2,])
    }
    dat_plot <- rbind(c(dat_plot[1,1], 0), dat_plot, c(dat_plot[nrow(dat_plot),1], 0))
    polygon(dat_plot, border = NA, col = conf_int_col, ...)
    lines(dat_plot, col = conf_int_border)
  }
}
