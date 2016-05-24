#' Plot the Effect of Fees
#'
#' @author George Fisher
#'
#' #@description
#'
#' #@details
#'
#' #@examples
#'
#' @return No explicit returned values, just a plot to output_path
#'
#' @param raw_data s&p500 daily data
#' @param annual_fee annual fee
#' @param output_path file path of plot
#'
#'@export
plot_effect_of_fees <- function(raw_data, annual_fee, output_path) {

    daily_fee <- annual_fee / mean(table(lubridate::year(raw_data$Date)))

    temp_df  <- data.frame(Date     = raw_data$Date,
                          Adj.Close = raw_data$Adj.Close,
                          chg       = calc_chg(raw_data$Adj.Close))

    temp_df$haircut <- raw_data$Adj.Close
    for (i in 2:nrow(temp_df)) {
        temp_df$haircut[i] <- temp_df$haircut[i-1] * temp_df$chg[i] * (1 - daily_fee)
    }

    start_year <- lubridate::year(temp_df$Date[1])
    end_year   <- lubridate::year(temp_df$Date[nrow(temp_df)])
    num_years  <- end_year - start_year + 1

    png(filename=paste0(output_path, "effect_of_fees.png"))

        plot(temp_df$Date, temp_df$Adj.Close,
             pch=20, type="l", xlab=NA, ylab=NA, yaxt="n",
             main=paste0("Effect of an ", annual_fee,
                         " Annual Fee over ", num_years, " years"))
        lines(temp_df$haircut, pch=20, col="blue")
        axis(2, las=2, at=axTicks(2), labels=fmt(axTicks(2)))
        grid(nx=length(axTicks(1))+3, ny=NULL)
        legend("topleft", legend=c("S&P 500 TR","S&P 500 TR minus fees"),
               lty=c(1,1), col=c("black","blue"))

    dev.off()

}
