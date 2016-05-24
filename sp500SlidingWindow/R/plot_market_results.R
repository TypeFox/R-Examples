#' Plot The Daily Chart For One Window
#'
#' @author George Fisher
#'
#' #@description
#'
#' #@details
#'
#' #@examples
#'
#' @export
#'
#'
#' @return No explicit returned values, merely the side-effect of a graphs
#' sent to the output_path folder
#'
#' @param filtered_data data.frame of the window's daily data
#' @param statistics summary statistics for this window
#' @param output_path file path to write the graph to
#'
plot_market_results <- function(filtered_data,
                                statistics,
                                output_path) {

    start_year <- statistics[['start_year']]
    end_year   <- statistics[['end_year']]
    ending_bal <- fmt(statistics[['ending_bal']])
    IRR        <- paste0(round(statistics[['IRR']]*100, 1),"%")

    plot_name <- paste0(start_year, "-", end_year, ".png")
    png(filename=paste0(output_path, plot_name))

    plot(filtered_data$Date, filtered_data$bal,
         pch=20, type="l", xlab=NA, ylab=NA, yaxt="n",
         main=paste0(start_year, " - ", end_year), cex.main=1.25,
         sub=paste0("Ending Bal: ", ending_bal, ";  IRR: ", IRR))

    if (filtered_data$bal[nrow(filtered_data)] <= 0) abline(h=0, col="red", lty=2)

    axis(2, las=2, at=axTicks(2), labels=fmt(axTicks(2)), cex.axis=0.72)
    grid(nx=length(axTicks(1))+3, ny=NULL)

    dev.off()
}
