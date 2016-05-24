#############################################################################
# utils.R
#############################################################################

utils.symbols.replace <- function(x,
                                  replacements)
{
    index <- x %in% names(replacements)

    ifelse(index, replacements[x], x)
}

#########
# format functions
#########
utils.get.format <- function(x) UseMethod("utils.get.format")

utils.get.format.default <- function(x)
{
    return("object")
}

utils.get.format.character <- function(x)
{
    return("character")
}

utils.get.format.numeric <- function(x)
{
    return("numeric")
}

utils.get.format.htest <- function(x)
{
    switch(EXPR = x$method,
           "One Sample t-test" = "t.test",
           "Welch Two Sample t-test" = "t.test",
           "Fisher's Exact Test for Count Data" = "fisher",
           "Pearson's product-moment correlation" = "pearson",
           "Pearson's Chi-squared test" = "chisquared",
           "Shapiro-Wilk normality test" = "shapiro",
           "Bartlett test of homogeneity of variances" = "bartlett",
           "Friedman rank sum test" = "bartlett",
           "Two-sample Kolmogorov-Smirnov test" = "ks.test",
           "object")
}

utils.get.format.summary.aovlist <- function(x)
{
    return("summary.aovlist")
}

utils.get.format.summary.lm <- function(x)
{
    return("summary.lm")
}
