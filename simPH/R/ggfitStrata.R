#' Graph fitted stratified survival curves from Cox Proportional Hazards models
#'
#' \code{ggfitStrata} graphs fitted survival curves created with
#' \code{\link{survfit}} using ggplot2.
#'
#' @param obj a \code{survfit} object.
#' @param byStrata logical, whether or not you want to include all of the
#' stratified survival curves on one plot or separate them into a grid arranged
#' plot.
#' @param xlab a label for the plot's x-axis
#' @param ylab a label of the plot's y-axis
#' @param title plot title.
#' @param lcolour line color. Currently only works if \code{byStrata = TRUE}.
#' The default it \code{lcolour = "#2C7FB8"} (a bluish hexadecimal colour)
#' @param rcolour confidence bounds ribbon color. Either a single color or a
#' vector of colours. The default it \code{lcolour = "#2C7FB8"}
#' (a bluish hexadecimal colour)
#'
#' @description This function largely improves \code{\link{plot.survfit}}. It
#' plots the curves using ggplot2 rather than base R graphics. One major
#' advantage is the ability to split the survival curves into multiple plots and
#' arrange them in a grid. This makes it easier to examine many strata at once.
#' Otherwise they can be very bunched up.
#'
#' Note: the strata legend labels need to be changed manually (see
#' \code{revalue}) in the \code{survfit} object with the \code{strata}
#' component.
#'
#' @examples
#' # Load packages
#' library(survival)
#' library(simPH)
#' 
#' # Subset data
#' bladder1 <- bladder[bladder$enum < 5, ]
#' 
#' # Estimate coxph model (note that this model is for code illustration only)
#' M1 <- coxph(Surv(stop, event) ~ (rx + size + number) + strata(enum) +
#'                 cluster(id), bladder1)
#' 
#' # Find predicted values
#' M1Fit <- survfit(M1)
#' 
#' # Plot strata in a grid
#' ggfitStrata(M1Fit, byStrata = TRUE)
#' 
#' # Plot all in one
#' ggfitStrata(M1Fit, byStrata = FALSE)
#' 
#' @seealso \code{\link{survfit}}, \code{ggplot2} and
#' \code{\link{strata}}
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @export

ggfitStrata <- function(obj, byStrata = FALSE, xlab = "", ylab = "", title = "",
                        lcolour = "#2C7FB8", rcolour = "#2C7FB8")
{
    Strata <- Lower <- Upper <- StrataC <- Time <- Survival <- NULL
    sFit <- obj
    time <- sFit$time
    lower <- sFit$lower
    upper <- sFit$upper
    S <- sFit$surv
    strata <- sFit$strata
    strata <- factor(rep(names(strata), strata), levels = names(strata))
    TempData <- data.frame(Time = time, Lower = lower,
                            Upper = upper, Survival = S,
                            Strata = strata)

    if (byStrata == FALSE){
        ggplot(data = TempData, aes(x = Time,
                                     y = Survival,
                                     color = Strata,
                                     fill = Strata)) +
        geom_line() +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = I(0.1)) +
        xlab(xlab) + ylab(ylab) + ggtitle(title) +
        theme_bw()

    } else if (byStrata == TRUE){
        TempData$StrataC <- gsub("=", "", TempData$Strata)
        TempData$StrataC <- gsub(" ", "", TempData$StrataC)
        eachStrata <- unique(TempData$StrataC)
        p <- list()
        for (i in eachStrata){
        SubData <- subset(TempData, StrataC == i)
        p[[i]] <- ggplot(data = SubData, aes(x = Time,
                                           y = Survival)) +
                                   geom_line(colour = lcolour) +
                                   geom_ribbon(aes(ymin = Lower,
                                                ymax = Upper),
                                                alpha = I(0.1),
                                                colour = NA,
                                                fill = rcolour) +
                                   xlab("") + ylab("") +
                                   ggtitle(paste(i, "\n")) +
                                   theme_bw()
        }
        Grid <- do.call(grid.arrange, c(p, top = title, bottom = xlab,
                        left = ylab))
    }
}
