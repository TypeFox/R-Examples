#' @title
#' Print Method for Seroincidence Summary Object
#'
#' @description
#' Custom \code{\link{print}} function to show output of the seroincidence summary \code{\link{summary.seroincidence}}.
#'
#' @param
#' x A list containing output of function \code{\link{summary.seroincidence}}.
#' @param
#' ... Additional arguments affecting the summary produced.
#'
#' @examples
#'
#' \dontrun{
#'
#' # estimate seroincidence
#' seroincidence <- estimateSeroincidence(...)
#'
#' # calculate summary statistics for the seroincidence object
#' seroincidenceSummary <- summary(seroincidence)
#'
#' # print the summary of seroincidence object to the console
#' print(seroincidenceSummary)
#'
#' # or simply type (appropriate print method will be invoked automatically)
#' seroincidenceSummary
#' }
#'
#' @export
print.summary.seroincidence <- function(x, ...)
{
    cat("Seroincidence estimated given the following setup:\n")
    cat(paste("a) Antibodies:    ", paste(x[["Antibodies"]], collapse = ", ")), "\n")
    cat(paste("b) Strata:        ", paste(x[["Strata"]], collapse = ", ")), "\n")
    censorLimits <- x[["CensorLimits"]]
    cat(paste("c) Censor limits: ", paste(sapply(names(censorLimits), FUN = function(name) {paste(name, censorLimits[name], sep = " = ")}), collapse = ", "), "\n"))
    cat(paste("d) Quantiles:     ", paste(x[["Quantiles"]], collapse = ", ")), "\n")

    cat("\n Seroincidence estimates:\n")
    print(x[["Results"]])
}
