#####
## Tests for summary tables
require("papeR")

if (require("nlme")) {
    ## Use dataset Orthodont
    data(Orthodont, package = "nlme")

    ## Get summary for continuous variables
    (tab1 <- summarize(Orthodont, type = "numeric"))
    get_option(tab1, "sep")
    summary(tab1)
    ## check handling of digits and sep
    (tab1a <- summarize(Orthodont, type = "numeric", digits = 3, sep = TRUE))
    get_option(tab1a, "sep")

    ## Change statistics to display
    summarize(Orthodont, quantiles = FALSE, type = "numeric")
    summarize(Orthodont, quantiles = FALSE, count = FALSE, type = "numeric")
    (tmp <- summarize(Orthodont, mean_sd = FALSE, type = "numeric"))

    ## Get summary for categorical variables
    (tab2 <- summarize(Orthodont, type = "fac"))
    get_option(tab2, "sep")
    summary(tab2)
    ## check handling of digits and sep
    (tab2a <- summarize(Orthodont, type = "fac", digits = 4, sep = FALSE))
    get_option(tab2a, "sep")

    ## use fraction instead of percentage
    summarize(Orthodont, percent = FALSE, type = "fac")

    ## try using the tables with Markdown
    if (require("knitr")) {
        kable(tab1)
        kable(tab2)
    }

    if (require("xtable")) {
        ans <- xtable(tab1)
        print(ans)
        ## grouped
        xtable(summarize(Orthodont, group = "Sex"))
        print(xtable(tab2))
    }
}

