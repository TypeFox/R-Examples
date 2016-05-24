#' Plot the Ending Balances and Withdrawals
#'
#' @author George Fisher
#'
#' #@description
#'
#' #@details
#'
#' #@examples
#'
#' @return No explicit returned values, creates the plots
#' 'ending_bal.png' and 'withdrawals.png' in the output_path
#'
#' @param window_df data.frame output of sp500SlidingWindow
#' @param output_path folder to write the file and the plot
#' @param lump_sum the amount initiall invested
#' @param annual_withdrawal the amount of hoped-for withdrawals
#' @param annual_inflation each year the withdrawal goes up this much
#' @param window_width the number of years in the sliding window
#'
#'@export
plot_ending_bal_withdrawals <- function(window_df, output_path, lump_sum,
                                        annual_withdrawal, annual_inflation,
                                        window_width) {

    # ================ collect stats =================

    #
    total_hoped_for_wd <- 0
    last <- annual_withdrawal
    for (i in 1:window_width) {
        last = last * (1 + annual_inflation)
        total_hoped_for_wd <- total_hoped_for_wd + last
    }
    ylim = c(0, total_hoped_for_wd)

    # ================ create a file with text stats =================

    # fileName <- paste0(output_path,"distribution.txt")
    # fileCon <- file(fileName, "at") #  opened for writing text
    #
    #     cat(file=fileCon, sep = "", "Worst Window ", wrst_window,
    #         " Best Window ", best_window, "\n")
    #
    # close(fileCon)

    # ================ plot the charts
    png(filename=paste0(output_path, "ending_bal.png"))

        barplot(window_df$ending_bal, ylab=NA, yaxt="n",
                names.arg=window_df$start_year,
                main=paste0("ending balances, ", window_width,"-year windows"),
                sub=paste0(annual_withdrawal))
        axis(2, las=2, at=axTicks(2), labels=fmt(axTicks(2)), cex.axis=0.72)


    dev.off()

    png(filename=paste0(output_path, "total withdrawals.png"))

        barplot(window_df$withdrawals, ylab=NA, yaxt="n",ylim=ylim,
                names.arg=window_df$start_year,
                main=paste0("total withdrawals, ", window_width,"-year windows"),
                sub=paste0(annual_withdrawal))
        axis(2, las=2, at=axTicks(2), labels=fmt(axTicks(2)), cex.axis=0.72)
        abline(h=total_hoped_for_wd, col="red", lty=2)

    dev.off()
}
