
library(gridGraphics)

plot.table1 <- function() {
    ## 1-d tables
    (Poiss.tab <- table(N = stats::rpois(200, lambda = 5)))
    plot(Poiss.tab, main = "plot(table(rpois(200, lambda = 5)))")
}

plot.table2 <- function() {
    plot(table(state.division))
}

plot.table3 <- function() {
    ## 4-D :
    plot(Titanic, main ="plot(Titanic, main= *)")
}

plotdiff(expression(plot.table1()), "plot.table-1")
plotdiff(expression(plot.table2()), "plot.table-2", width=15)
plotdiff(expression(plot.table3()), "plot.table-3")

plotdiffResult()
