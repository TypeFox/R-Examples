
library(gridGraphics)

## Use the Berkeley admission data as in Friendly (1995).
x <- aperm(UCBAdmissions, c(2, 1, 3))
dimnames(x)[[2]] <- c("Yes", "No")
names(dimnames(x)) <- c("Sex", "Admit?", "Department")

fourfold1 <- function() {
    ## Fourfold display of data aggregated over departments, with
    ## frequencies standardized to equate the margins for admission
    ## and sex.
    ## Figure 1 in Friendly (1994).
    fourfoldplot(margin.table(x, c(1, 2)))
}

fourfold2 <- function() {
    ## Fourfold display of x, with frequencies in each table
    ## standardized to equate the margins for admission and sex.
    ## Figure 2 in Friendly (1994).
    fourfoldplot(x)
}

fourfold3 <- function() {
    ## Fourfold display of x, with frequencies in each table
    ## standardized to equate the margins for admission. but not
    ## for sex.
    ## Figure 3 in Friendly (1994).
    fourfoldplot(x, margin = 2)
}

plotdiff(expression(fourfold1()), "fourfold-1")
plotdiff(expression(fourfold2()), "fourfold-2")
plotdiff(expression(fourfold3()), "fourfold-3")

plotdiffResult()
