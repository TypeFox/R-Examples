### ===== actuar: An R Package for Actuarial Science =====
###
### Normal and Normal Power Approximation of the total amount of
### claims distribution
###
### See Dayken, Pentikanen and Pesonen, Practical Risk Theory for
### Actuaries, Chapman & Hall, 1994.
###
### AUTHORS:  Vincent Goulet <vincent.goulet@act.ulaval.ca>
### and Louis-Philippe Pouliot

normal <- function(mean, variance)
{
    ## Approximate the total amount of claims distribution using the first
    ## two moments.
    FUN <- function(x) pnorm(x, mean = mean, sd = sqrt(variance))

    environment(FUN) <- new.env()
    assign("mean", mean, envir = environment(FUN))
    assign("variance", variance, envir = environment(FUN))
    attr(FUN, "source") <- "function(x) pnorm(x, mean = mean, sd = sqrt(variance))"
    FUN
}

npower <- function(mean, variance, skewness)
{
    ## Approximate the total amount of claims distribution using the first
    ## three moments.
    FUN <- function(x)
        ifelse(x <= mean, NA,
               pnorm(sqrt(1 + 9/skewness^2 + 6 * (x - mean)/(sqrt(variance) *
                                                             skewness)) -
                     3/skewness))

    environment(FUN) <- new.env()
    assign("mean", mean, envir = environment(FUN))
    assign("variance", variance, envir = environment(FUN))
    assign("skewness", skewness, envir = environment(FUN))
    attr(FUN, "source") <- "function(x) ifelse(x <= mean, NA, pnorm(sqrt(1 + 9/skewness^2 + 6 * (x - mean)/(sqrt(variance) * skewness)) - 3/skewness))"
    FUN
}
