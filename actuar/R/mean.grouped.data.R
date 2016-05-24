### ===== actuar: An R Package for Actuarial Science =====
###
### Mean (TODO: and summaries) of grouped data objects
###
### See Klugman, Panjer & Willmot, Loss Models, Wiley, 1998.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

mean.grouped.data <- function(x, ...)
{
    ## Get group boundaries
    cj <- eval(expression(cj), envir = environment(x))

    ## Compute group midpoints
    midpoints <- cj[-length(cj)] + diff(cj)/2

    ## Drop the boundaries column and convert to matrix for use in
    ## crossprod()
    x <- as.matrix(x[-1])

    ## Compute mean per column
    drop(crossprod(x, midpoints))/colSums(x)
}
