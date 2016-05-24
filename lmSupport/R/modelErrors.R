modelErrors <- function(Model)
{
    cat(sprintf('Sum of Squared Errors (SSE):  %.2f\n\n', sum(residuals(Model)^2)), sep='')
    print(residuals(Model))
    invisible(residuals(Model))
}