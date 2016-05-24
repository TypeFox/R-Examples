#' Calculates the required sample size for a one-way ANOVA.
#'
#' @param iv List of data for the treatment to be tested.
#' @param sig.level Desired significance level (default is 0.05).
#' @param power Desired level of power (default is 0.80).
#' @return Returns the recommended sample size given the conditions to achieve the desired \code{power}.
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences} (2nd ed.). Hillsdale, N.J.: Lawrence Erlbaum Associates.
#' @examples
#' # Exercise 8.10, p.391 from Cohen (1988)
#' main.eff <- list(name = "Teaching", levels = 4, eta.sq = 0.0588)
#' # Running the function with default settings
#' n.oneway(iv = main.eff)
#' @export
n.oneway <- function(iv=iv, sig.level = 0.05, power = 0.80) {
    name.iv <- iv$name
    effect.levels <- iv$levels

    # Obtain user effect size if needed
    if(length(iv$eta.sq) != 0) {
        effect.size <- iv$eta.sq
    }
    else {
        effect.size <- 0.01
    }

    # Checking user input for errors
    if(iv$levels < 2) {
        stop("There must be atleast 2 levels in a treatment.")
    }
    es.list <- list("small", "med", "large")
    if(is.character(iv$eta.sq) == TRUE) {
        if(!(iv$eta.sq %in% es.list)) {
            stop("Accepted string values are: 'small', 'med', 'large'.")
        }
    }
    else if(is.numeric(iv$eta.sq)) {
        if(iv$eta.sq <= 0) {
            stop("The effect size must be greater than 0.")
        }
    }
    else if(power <= 0) {
        stop("Power must be greater than 0.")
    }
    else if(sig.level <= 0) {
        stop("The significance level must be greater than 0.")
    }

    # Convert the eta-squared value to f
    f.value <- f.oneway(effect.size = effect.size)

    # Sample size calculations
    sample.size <- pwr::pwr.anova.test(k = effect.levels, f = f.value, sig.level = sig.level, power = power)

    out.data <- get.oneway(sample.size = sample.size, name.iv = name.iv)
    cat(sprintf("\nThe following is the sample size recommendation to \nachieve a power of %1.2f at a significance level of %1.2f.\n", power, sig.level))
    return(out.data)
}
