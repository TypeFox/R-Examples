
library(gridGraphics)

assocplot1 <- function() {
    ## Aggregate over sex:
    x <- margin.table(HairEyeColor, c(1, 2))
    assocplot(x, main = "Relation between hair and eye color")
}

plotdiff(expression(assocplot1()), "assocplot-1")

plotdiffResult()
