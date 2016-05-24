#' Sliding Window Investment Analysis
#'
#' @author George Fisher
#'
#' @description Uses the S&P500 daily data on a series of windows and generates
#' graphs & statistics on the performance.
#'
#' @details The daily market data for the S&P 500 from 1950 to the present
#' is broken into a series of periods or windows of equal length (except for
#' the last period). The investment and expense data provided is analyzed to
#' see how it would fare in each widow. Graphs and statistics are produced.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' investment_vector <- seq(1,30)*10
#' withdrawal_vector <- c(investment_vector[1:10]  * 0.15,
#'                    investment_vector[11:20] * 0.35,
#'                    investment_vector[21:30] * runif(10, min=0.01, max=0.90))
#' window_df <- sp500SlidingWindow(investment_vector,
#'                                withdrawal_vector, output_path='~/Downloads/')
#' }
#'
#' @return data.frame with summary statistics for each window plus
#'  the side-effects of graphs written to the output_path given
#'
#' @param investment_vector a vector of annual investments,
#' (positive values represent investment into the account), length must be window_width
#' @param withdrawal_vector a vector of annual withdrawals,
#' (positive values represent withdrawal from the account), length must be window_width
#' @param window_width the number of years in each window
#' @param annual_fee the total annual percent removed by the investment managers
#' @param output_path file path to a folder in which graphs & statistics will be saved
#'
#' @importFrom graphics axTicks
#' @importFrom graphics axis
#' @importFrom graphics grid
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @importFrom graphics barplot
#' @importFrom graphics points
#' @importFrom graphics text
#' @importFrom graphics abline
#' @importFrom graphics hist
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom magrittr %>%
#'
sp500SlidingWindow <- function(investment_vector,
                               withdrawal_vector,
                               window_width = 30,
                               annual_fee   = 0.0125,
                               output_path  = "~/") {
    Date <- NULL

    # ======================= error checks =======================

    # make sure things are the right length
    if (length(investment_vector) != length(withdrawal_vector) ||
        length(investment_vector) != window_width)
        stop('input vector length must equal window_width')

    # make sure output_path ends in a slash
    if (substr(output_path, nchar(output_path), nchar(output_path)) != '/')
        output_path <- paste0(output_path, '/')

    # ================ plot the investments and withdrawals ================
    plot_investments_expense(investment_vector, withdrawal_vector, output_path)

    # ================ read in and sort the stock market data ================
    raw_data <- SP500TR_1950() %>% dplyr::arrange(Date)

    raw_data$Year <- lubridate::year(raw_data$Date)

    # ================ plot the total stock market ================
    plot_total_stock_market(raw_data, output_path)

    # ================ loop through each window ================
    net_investment <- investment_vector - withdrawal_vector
    start_year     <- lubridate::year(raw_data$Date[1])
    window_df      <- NULL
    windows_procd  <- 0

    while (TRUE) {
        end_year <- start_year + (window_width-1)
        if (end_year > lubridate::year(lubridate::now())) break

        # select the range of years for this window
        # =========================================
        filtered_data <- raw_data %>%
            dplyr::filter(lubridate::year(Date) >= start_year &
                          lubridate::year(Date) <= end_year)

        # calculate the change in the Adj.Close values
        # ============================================
        filtered_data$chg <- calc_chg(filtered_data$Adj.Close)

        # add investments and withdrawals at the start of each year
        # ========================================================
        # indexes of the first of each year
        year_idx <- c(1, 1+which(diff(filtered_data$Year)!=0))

        filtered_data$inv <- 0
        filtered_data$inv[year_idx] <- investment_vector

        filtered_data$wdr <- 0
        filtered_data$wdr[year_idx] <- withdrawal_vector

        filtered_data$net <- 0
        filtered_data$net[year_idx] <- net_investment

        # what's each day's fee?
        # ======================
        daily_fee <- annual_fee / mean(table(lubridate::year(filtered_data$Date)))

        # calculate the daily account balance
        # ===================================
        filtered_data$bal    <- NA
        filtered_data$bal[1] <- filtered_data$inv[1] -
                                filtered_data$wdr[1]

        for(i in 2:nrow(filtered_data)) {
            # previous balance
            #   times the change in today's value
            #      haircut by the fee
            # plus any new net investment
            prev_bal <- filtered_data$bal[i-1]
            filtered_data$bal[i] <-
                ((prev_bal * filtered_data$chg[i]) * (1-daily_fee))

            filtered_data$bal[i] <- filtered_data$bal[i] +
                                    filtered_data$inv[i]

            # if the withdrawal amount is less than the
            # current balance, set withdrawal equal to balance
            # and set the balance equal to zero
            withdrawal <- filtered_data$wdr[i]
            if (withdrawal >= filtered_data$bal[i]) {
                withdrawal           <- filtered_data$bal[i]
                filtered_data$bal[i] <- 0
            } else {
                filtered_data$bal[i] <- filtered_data$bal[i] - withdrawal
            }

            filtered_data$wdr[i] <- withdrawal
        }

        # ============== append to data.frame ================
        ending_bal <- filtered_data$bal[nrow(filtered_data)]

        cash_flow                    <- filtered_data$wdr[year_idx] -
                                        filtered_data$inv[year_idx]
        cash_flow[length(cash_flow)] <- cash_flow[length(cash_flow)] + ending_bal

        irr_err_msg <- paste("The FinCal::irr funtion in the sp500SlidingWindow function",
                      "has raised an exception. You should look at the",
                      "investment and expense vectors to be sure you",
                      "are not setting up an impossible situation.",
                      "The image already created at",
                             paste0(output_path, 'investments_expense.png'),
                             "may help.")
        IRR <- tryCatch(
            {
                FinCal::irr(cf = cash_flow)
            },
            error=function(e) {
                #message(irr_err_msg)
                #stop(paste("IRR Error", e))
                return(NA)
            },
            warning=function(w) {
                #message(irr_err_msg)
                #stop(paste("IRR Warning", w))
                return(NA)
            },
            finally={}
            )

        window_df <- rbind(window_df,
                          data.frame(start_year = start_year,
                                     end_year   = end_year,
                                     ending_bal = ending_bal,
                                     inv        = sum(filtered_data$inv),
                                     wdr        = sum(filtered_data$wdr),
                                     IRR        = IRR))

        # ============== plot the market results for this window ================
        plot_market_results(filtered_data, window_df[nrow(window_df),], output_path)

        # ============== do it again ================
        start_year    <- start_year + 1
        windows_procd <- windows_procd + 1
    }

    # ================= create stats file and plot distribution of ending balances
    plot_ending_bal_distribution(window_df, output_path, window_width)
    fileName <- paste0(output_path, "distribution.txt")
    fileCon <- file(fileName, "at") #  opened for appending text

        cat(file=fileCon, sep = "", "Total Invested ",  fmt(sum(investment_vector)), "\n")
        cat(file=fileCon, sep = "", "Total Withdrawn ", fmt(sum(withdrawal_vector)), "\n")
        cat(file=fileCon, sep = "", "Total Net Investment ", fmt(sum(net_investment)), "\n")
        cat(file=fileCon, sep = "", "Mean Net Investment ",  fmt(mean(net_investment)), "\n")
        cat(file=fileCon, sep = "", "Total windows processed ", windows_procd, "\n")

    close(fileCon)

    # ================= plot the effect of fees
    plot_effect_of_fees(raw_data, annual_fee, output_path)

    return(window_df)
}
