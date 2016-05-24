#' @title
#' Print Method for Seroincidence Object
#'
#' @description
#' Custom \code{\link{print}} function to show output of the seroincidence calculator \code{\link{estimateSeroincidence}}.
#'
#' @param
#' x A list containing output of function \code{\link{estimateSeroincidence}}.
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
#' # print the seroincidence object to the console
#' print(seroincidence)
#'
#' # or simply type (appropriate print method will be invoked automatically)
#' seroincidence
#' }
#'
#' @export
print.seroincidence <- function(x, ...)
{
    cat("Seroincidence object estimated given the following setup:\n")
    cat(paste("a) Antibodies:    ", paste(x[["Antibodies"]], collapse = ", ")), "\n")
    cat(paste("b) Strata:        ", paste(x[["Strata"]], collapse = ", ")), "\n")
    censorLimits <- x[["CensorLimits"]]
    cat(paste("c) Censor limits: ", paste(sapply(names(censorLimits), FUN = function(name) {paste(name, censorLimits[name], sep = " = ")}), collapse = ", "), "\n"))

    cat("\n")
    cat("This object is a list containing the following items:\n")
    cat("Fits - List of outputs of optim function per stratum.\n")
    cat("Antibodies - Input parameter antibodies of function estimateSeroincidence.\n")
    cat("Strata - Input parameter strata of function estimateSeroincidence.\n")
    cat("CensorLimits - Input parameter censorLimits of function estimateSeroincidence.\n")

    cat("\n")
    cat("Call summary function to obtain output results.\n")
}
