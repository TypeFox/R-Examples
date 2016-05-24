#' Dynamically simulate one scenario
#'
#' \code{OneScen} is an internal function to dynamically simulate one scenario.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats coef quantile vcov
#'
#' @keywords internals
#' @noRd

OneScen <- function(obj, ldv, n, scen, sig, num, shocks){
    # CRAN requirements
    times <- sigma.sqr <- alpha.sqr <- NULL

    # Create lower and upper percentile bounds of the confidence interval
    Bottom <- (1 - sig)/2
    Top <- 1 - Bottom

    # Create data frame to fill in with simulation summaries
    SimSum <- data.frame()
    ShockVals <- data.frame()

    # Change data frame for shock values
    for (i in 1:n) {
        if (is.null(shocks))  {
            scenTemp <- scen
        }
        else if (!is.null(shocks)) {
            if (i %in% shocks[, "times"]) {
                scenTemp <- scen
                for (x in names(shocks)[-1]) {
                    shocksTemp <- subset(shocks, times == i)
                    scenTemp[, x] <- shocksTemp[1, x]
                }
            }
            else if (!(i %in% shocks)) {
                scenTemp <- scen
            }
        }

        # Parameter estimates & Variance/Covariance matrix
        Coef <- coef(obj)
        VC <- vcov(obj)

        # Draw covariate estimates from the multivariate normal distribution
        Drawn <- mvrnorm(n = num, mu = Coef, Sigma = VC)
        DrawnDF <- data.frame(Drawn)

        ##### The intercept
        #names(DrawnDF) <- names(scenTemp)

        # Find predicted value
        qiDF <- data.frame(DrawnDF[, 1]) # keep the intercept
        for (u in names(scenTemp)) {
            qiDF[, u] <- DrawnDF[, u] * scenTemp[, u]
        }
        PV <- rowSums(qiDF)

        # Create summary data frame
        time <- i
        ldvMean <- mean(PV)

        ldvLower <- quantile(PV, prob = Bottom, names = FALSE)
        ldvUpper <- quantile(PV, prob = Top, names = FALSE)
        ldvLower50 <- quantile(PV, prob = 0.25, names = FALSE)
        ldvUpper50 <- quantile(PV, prob = 0.75, names = FALSE)

        # Shock variable values
        if (!is.null(shocks)) {
            ShockNames <- names(shocks)[-1]
            TempShock <- scenTemp[, ShockNames]
            ShockVals <- rbind(ShockVals, TempShock)
        }
        # Combine
            TempOut <- cbind(time, ldvMean, ldvLower, ldvUpper, ldvLower50,
                            ldvUpper50)
            SimSum <- rbind(SimSum, TempOut)

            # Change lag variable for the next simulation
            scen[, ldv] <- ldvMean
    }
    # Clean up shocks
    if (!is.null(shocks)) {
        CleanNames <- paste0("shock.", ShockNames)
        names(ShockVals) <- CleanNames
        SimSum <- cbind(ShockVals, SimSum)
    }
    col_idx <- grep("time", names(SimSum))
    SimSum <- SimSum[, c(col_idx, (1:ncol(SimSum))[-col_idx])]
    return(SimSum)
}

#' Extract legend from ggplot2 object
#'
#' @source Hadley Wickham
#' @import ggplot2
#' @keywords internals
#' @noRd

gLegend <- function(Plot){
    tmp <- ggplot_gtable(ggplot_build(Plot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}
