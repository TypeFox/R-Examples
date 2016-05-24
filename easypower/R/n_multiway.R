#' Sample size calculations for factorial ANOVAs
#'
#' @param iv1 The list of data for treatment 1.
#' @param iv2 The list of data for treatment 2.
#' @param iv3 (optional) The list of data for treatment 3.
#' @param iv4 (optional) The list of data for treatment 4.
#' @param interaction.eta2 (optional) Either a character string or numeric value of the desired eta squared. Default is set to "small".
#' @param sig.level (optional) Desired significance level. Default value is 0.05.
#' @param power (optional) Desired level of power. Default value is 0.80.
#' @param result The amount of data that will be output to the user (\emph{default = "all"}).
#' The following are the three output options the user may specify:
#' \itemize{
#'  \item result = "all" - Outputs the sample size recommendations for all treatments and all possible interactions.
#'  \item result = "highest" - Outputs the highest recommended sample size.
#'  \item result = "select" - Outputs specific results to the user. If there has been previous research on an
#'  effect, the user may input a numeric value for the effect size. The output will consist of the highest recommended sample size
#'  and the recommendations where the user input a numeric value for the effect size of a treatment.
#'}
#' @param ... Extra interactions to pass in. In order to change the effect size of a specific interaction an
#' interaction effect may be added to the function. It must take the form: \emph{int# = int.eff.#}.
#' @details
#' Acceptable effect size character string values and their numeric equivalents are: "small" (0.01), "med" (0.06), and "large" (0.14).
#' @note
#' Sample size recommendations are rounded up to the nearest integer. More detailed examples on n.multiway can be viewed in the vignette.
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences} (2nd ed.). Hillsdale, N.J.: Lawrence Erlbaum Associates.
#' @examples
#' # Exercise 8.15, p.400 from Cohen (1988)
#' # Defining the treatments
#' main.eff.1 <- list(name = "R", levels = 2, eta.sq = 0.123)
#' main.eff.2 <- list(name = "C", levels = 4, eta.sq = 0.215)
#' # Running n.multiway
#' n.multiway(iv1=main.eff.1, iv2=main.eff.2, interaction.eta2 = 0.079)
#' # To just view highest
#' n.multiway(iv1=main.eff.1, iv2=main.eff.2, interaction.eta2 = 0.079, result = "highest")
#'
#' # Exercise 8.14, p.397 from Cohen (1988)
#' # Defining the treatments and interaction
#' main.eff.1 <- list(name = "Sex", levels = 2, eta.sq = 0.0099)
#' main.eff.2 <- list(name = "Age", levels = 3, eta.sq = 0.0588)
#' main.eff.3 <- list(name = "Conditions", levels = 4, eta.sq = 0.1506)
#' # Running n.multiway
#' n.multiway(iv1=main.eff.1, iv2=main.eff.2, iv3=main.eff.3, interaction.eta2 = 0.0588)
#'@export
n.multiway <- function(iv1=NULL, iv2=NULL, iv3=NULL, iv4=NULL, interaction.eta2="small", sig.level=0.05, power=0.80, result="all",...) {
    # Extract user information from input
    name.ivs <- c(iv1$name, iv2$name, iv3$name, iv4$name)
    effect.lev <- c(iv1$levels, iv2$levels, iv3$levels, iv4$levels)
    ivs <- name.ivs
    num.ivs <- length(name.ivs)

    # Add default effect sizes and insert user's IV effect sizes
    all.ivs <- list(iv1, iv2, iv3, iv4)
    es.strings <- c(rep("small", num.ivs))
    for(i in 1:num.ivs) {
        if(length(all.ivs[[i]]$eta.sq) != 0) {
            es.strings[i] <- all.ivs[[i]]$eta.sq
        }
    }

    # Check user input for errors in main effects
    for(i in 1:num.ivs) {
        if(effect.lev[i] < 2) {
            stop("The number of levels in a treatment must be 2 or greater.")
        }
        if(is.character(all.ivs[[i]]$eta.sq) == TRUE) {
            if(!(all.ivs[[i]]$eta.sq %in% list("small", "med", "large"))) {
                stop("Accepted string values are: 'small', 'med', 'large'.")
            }
        }
        else if(is.numeric(all.ivs[[i]]$eta.sq) == TRUE) {
            if(all.ivs[[i]]$eta.sq <= 0) {
                stop("The effect size must be greater than 0.")
            }
        }
        else if(power <= 0) {
            stop("Power must be greater than 0.")
        }
        else if(sig.level <= 0) {
            stop("The significance level must be greater than 0.")
        }
    }

    # Get total number of interactions
    count.variables <- sapply(1:num.ivs, FUN = function(x) {
        choose(num.ivs, x)
    }
    )

    total.variables <- sum(count.variables)
    total.interactions <- total.variables - num.ivs

    df.from.lev <- sapply(1:num.ivs, FUN = function(x) {
        c(effect.lev[x] - 1)
    }
    )

    # Create vectors with variable names, interactions, and df
    variable.names <- var.names(num.ivs, ivs)
    df <- df.vector(num.ivs, df.from.lev)

    # Make a vector of the effect size strings for the interactions
    es.strings.ints <- c(rep(interaction.eta2, total.interactions))
    all.es.strings <- c(es.strings, es.strings.ints)

    my.frame <- data.frame("variable" = variable.names, "df" = df, "eta.sq" = all.es.strings)

    # Get interaction data from user
    extra.ints <- list(...)
    count.extra.ints <- length(extra.ints)

    # Inserting the new effect sizes for the desired interactions
    if(count.extra.ints > 0) {
        # Check interaction data for errors
        for(i in 1:count.extra.ints) {
            if(is.character(extra.ints[[i]]$eta.sq) == TRUE & !(extra.ints[[i]]$eta.sq %in% list("small", "med", "large"))) {
                stop("Accepted string values are: 'small', 'med', 'large'.")
            }
            else if(is.numeric(extra.ints[[i]]$eta.sq) == TRUE & extra.ints[[i]]$eta.sq <= 0) {
                stop("The effect size for the interaction must be greater than 0.")
            }
        }
        my.frame$variable <- as.character(my.frame$variable)
        my.frame$eta.sq <- as.character(my.frame$eta.sq)
        for(i in 1:count.extra.ints) {
            if(extra.ints[[i]]$name %in% my.frame$variable) {
                my.frame$eta.sq[my.frame$variable == extra.ints[[i]]$name]  <- extra.ints[[i]]$eta.sq
            }
            else {
                stop(sprintf("The interaction term %s does not exist.", extra.ints[[i]]$name))
            }
        }
        my.frame$variable <- as.factor(my.frame$variable)
        my.frame$eta.sq <- as.factor(my.frame$eta.sq)
    }

    # Create f2 vector
    f2 <- f2.vector(eta.sq = as.vector(my.frame$eta.sq))
    my.frame <- data.frame(my.frame, "f2" = f2)

    # Create vector with levels for all treatments and interactions
    total.lev <- lev.vector(num.ivs, effect.lev)
    my.frame <- data.frame(my.frame, "lev" = total.lev)

    # Power calculations
    # Sample size calculations based on user data
    final.sample <- sample.size(my.frame = my.frame, sig.level = sig.level, power = power)

    final.frame <- data.frame(my.frame, "sample.size" = final.sample)
    sort.frame <- final.frame[order(-final.frame$sample.size), ]

    # Output results to user
    if(result == "highest") {
        out <- get.highest(sort.frame=sort.frame)
        cat(sprintf("\nThe following is the largest recommended total sample size.\n\n"))
        cat(sprintf("Desired power: %1.2f\nSignificance level: %1.2f\n", power, sig.level))
        cat(sprintf("Effect size used in calculations: Cohen's f-squared\n"))
        cat(sprintf("Cutoffs: small = 0.01, med = 0.06, large = 0.14\n"))
        return(out)
    }
    else if(result == "all") {
        out <- get.all(final.frame)
        cat(sprintf("\nThe following sample size recommendations are for each treatment and all possible interactions.\n"))
        cat(sprintf("Sample sizes are calculated independently using the estimated effect size to achieve \nthe desired power level.\n\n"))
        cat(sprintf("Desired power: %1.2f\nSignificance level: %1.2f\n", power, sig.level))
        cat(sprintf("Effect size used in calculations: Cohen's f-squared\n"))
        cat(sprintf("Cutoffs: small = 0.01, med = 0.06, large = 0.14\n\n"))
        names(out) <- gsub("\\.", " ", names(out))
        print(out, row.names = FALSE)
        cat(sprintf("\n"))
    }
    else if(result == "select") {
        out <- get.select(sort.frame = sort.frame)
        cat(sprintf("\nThe following is the highest sample size required and the sample size \nrecommendations where a numeric value for effect size was entered.\n"))
        cat(sprintf("Sample sizes are calculated independently using the estimated \neffect size to achieve the desired power level.\n\n"))
        cat(sprintf("Desired power: %1.2f\nSignificance level: %1.2f\n", power, sig.level))
        cat(sprintf("Effect size used in calculations: Cohen's f-squared\n"))
        cat(sprintf("Cutoffs: small = 0.01, med = 0.06, large = 0.14\n\n"))
        names(out) <- gsub("\\.", " ", names(out))
        print(out, row.names = FALSE)
        cat(sprintf("\n"))
    }
}
