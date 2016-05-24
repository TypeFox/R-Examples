#' Plot the Distribution of Ending Balances
#'
#' @author George Fisher
#'
#' #@description
#'
#' #@details
#'
#' #@examples
#'
#' @return No explicit returned values, (1) writes 'distribution.txt' with
#' statistics about the windows (2) creates a plot 'distribution.png' in the
#' output_path showing the distribution of outcomes
#'
#' @param window_df data.frame output of sp500SlidingWindow
#' @param output_path folder to write the file and the plot, if NULL
#' then the image is plotted to the screen
#' @param window_width number of years in each window
#'
#'@export
plot_ending_bal_distribution <- function(window_df, output_path, window_width) {

    # ================ collect stats =================

    min_end_bal     <- min(window_df$ending_bal)
    max_end_bal     <- max(window_df$ending_bal)

    avg_end_bal     <- mean(window_df$ending_bal, na.rm = FALSE)
    bal_hist_obj    <- hist(window_df$ending_bal, plot=FALSE)
    mod_end_bal     <- bal_hist_obj$mids[which.max(bal_hist_obj$counts)]

    mod_window_idx  <- which.min(abs(window_df$ending_bal - mod_end_bal))
    best_window_idx <- which(window_df$ending_bal==max(window_df$ending_bal))
    wrst_window_idx <- which(window_df$ending_bal==min(window_df$ending_bal))

    # find the year with the ending balance closest to hist's modal value
    closest_mod_bal <- window_df$ending_bal[mod_window_idx]

    avg_irr         <- mean(window_df$IRR, na.rm = TRUE)
    min_irr         <- min(window_df$IRR)
    max_irr         <- max(window_df$IRR, na.rm = TRUE)
    closest_mod_irr <- window_df$IRR[mod_window_idx]

    best_window <- paste0(window_df$start_year[best_window_idx], "-",
                         window_df$end_year[best_window_idx])
    wrst_window <- paste0(window_df$start_year[wrst_window_idx], "-",
                         window_df$end_year[wrst_window_idx])

    # ================ create a file with text stats =================

    fileName <- paste0(output_path,"distribution.txt")
    fileCon <- file(fileName, "wt") #  opened for writing text

        cat(file=fileCon, sep = "", "Worst Window ", wrst_window,
            " Best Window ", best_window, "\n")
        cat(file=fileCon, sep = "", "Ending Balance Range [", fmt(min_end_bal),
            " - ", fmt(max_end_bal), "]\n")
        cat(file=fileCon, sep = "", "IRR Range [", round(min_irr*100,2), "% - ",
            round(max_irr*100,2), "%]\n")
        cat(file=fileCon, sep = "", "Mean IRR ",  round(avg_irr*100, 1), "%",
            "; Modal IRR ", round(closest_mod_irr*100, 1), "%\n")
        cat(file=fileCon, sep = "", "Mean Ending Balance: ",  fmt(avg_end_bal),
            "; Modal Ending Balance: ", fmt(closest_mod_bal), "\n")

    close(fileCon)

    # ================ plot the distribution =================

    mod_freq <- max(bal_hist_obj$counts)

    xlab <- paste0 (
        paste0("Mean Ending Balance: ",  fmt(avg_end_bal),
               "; Modal Ending Balance: ", fmt(closest_mod_bal)), "\n",
        paste0("Mean IRR ",  round(avg_irr*100, 1), "%",
               "; Modal IRR ", round(closest_mod_irr*100, 1), "%")
    )

    if (!is.null(output_path))
        png(filename=paste0(output_path, "distribution.png"))

        hist(window_df$ending_bal,
             main=paste0("Distribution of Ending Balances\n",
                         window_df$start_year[1], " to ",
                         window_df$end_year[nrow(window_df)],
                         ", ", window_width, "-Year Sliding Window"),
             xaxt="n",xlab=NA, sub=xlab)
        axis(1, las=0, at=axTicks(1), labels=fmt(axTicks(1)))
        abline(v=as.numeric(sub(",", "", avg_end_bal, fixed = TRUE)),lty=3)
        text(as.numeric(sub(",", "", avg_end_bal, fixed = TRUE)), mod_freq-4, "mean")
        text(mod_end_bal, mod_freq-2, "mode")

        if (!is.null(output_path)) dev.off()
}
