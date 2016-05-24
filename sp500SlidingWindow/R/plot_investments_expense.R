#' Plot the Investments and Expenses
#'
#' @author George Fisher
#'
#' @description Produces a plot of the investments, expenses and net investment.
#'
#' #@details
#'
#' @export
#'
#' @examples
#' \dontrun{
#' investment_vector <- seq(1,30)*10
#' expense_vector <- c(investment_vector[1:10]  * 0.15,
#'                    investment_vector[11:20] * 0.35,
#'                    investment_vector[21:30] * runif(10, min=0.01, max=0.90))
#' plot_investments_expense(investment_vector,
#'                          expense_vector,
#'                          output_path=NULL)
#' }
#'
#' @return No explicit returned values, merely the side-effects of a graph
#' named 'investments_expense.png' written to output_path
#'
#' @param investment_vector a vector of annual investments
#' @param expense_vector a vector of annual withdrawals
#' @param output_path file path to a folder in which the graph will be saved,
#' NULL indicates plot, don't save
#'
plot_investments_expense <- function(investment_vector,
                                     expense_vector,
                                     output_path) {

    net_investment <- investment_vector - expense_vector

    ymax <- max(investment_vector,
               expense_vector,
               net_investment)*1.10
    ymin <- min(investment_vector,
               expense_vector,
               net_investment)*0.90

    if (!is.null(output_path)) {
        # make sure output_path ends in a slash
        if (substr(output_path, nchar(output_path), nchar(output_path)) != '/')
            output_path <- paste0(output_path, '/')

        png(filename=paste0(output_path, 'investments_expense.png'))
    }

    plot(investment_vector, pch=20, col="blue",
         ylim=c(ymin, ymax),
         xlab="Year", ylab=NA, yaxt="n",
         main="Investments and Withdrawals")
    axis(2, las=2, at=axTicks(2), labels=fmt(axTicks(2)), cex.axis=0.72)
    grid()

    points(expense_vector, pch=20, col="red")
    points(net_investment, pch=10, col="black", cex=0.75)

    legend("center", legend=c("Investments","Withdrawals","Net Investment"),
           pch=c(20,20, 10), col=c("blue","red","black"))

    if (!is.null(output_path))
        dev.off()

}
