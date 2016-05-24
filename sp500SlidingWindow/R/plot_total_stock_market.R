#' Plot the Entire Stock Market History from 1950
#'
#' @author George Fisher
#'
#' #@description
#'
#' #@details
#'
#' #@examples
#'
#' @return No explicit returned values, creates a plot
#' paste0("SP500_", start_year, "-", end_year, ".png") in the
#' output_path
#'
#' @param raw_data output of SP500TR_1950()
#' @param output_path folder to write the file and the plot
#'
#'@export
plot_total_stock_market <- function(raw_data, output_path) {

    fractional_years <- length(unique(raw_data$Year))-1 +
                        lubridate::month(lubridate::now()) / 12

    cagr <- CAGR(raw_data$Adj.Close[1],
                 raw_data$Adj.Close[nrow(raw_data)],
                 fractional_years, type = "geometric")

    start_year <- lubridate::year(raw_data$Date[1])
    end_year   <- lubridate::year(raw_data$Date[nrow(raw_data)])

    plot_name = paste0("SP500_", start_year, "-", end_year, ".png")
    png(filename=paste0(output_path, plot_name))

        plot(raw_data$Date, raw_data$Adj.Close, type="l", ylab=NA, yaxt="n",
             main=paste("Performance of the S&P 500 with Dividends\n",
                        start_year, " - ", end_year),
             xlab=paste0("Compound Annual Growth Rate ", round(cagr*100,1),"%"))
        axis(2, las=2, at=axTicks(2), labels=fmt(axTicks(2)))
        grid()

    dev.off()
}
