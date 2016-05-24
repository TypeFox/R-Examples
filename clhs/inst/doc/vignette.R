### R code from vignette source 'vignette.Rnw'

###################################################
### code chunk number 1: make_things_reproducible
###################################################
set.seed(42)


###################################################
### code chunk number 2: load_diamonds
###################################################
data(diamonds, package = 'ggplot2')
diamonds <- data.frame(diamonds)
head(diamonds)
nrow(diamonds)


###################################################
### code chunk number 3: simple_clhs
###################################################
library(clhs)
res <- clhs(diamonds, size = 100, progress = FALSE, iter = 1000)
str(res)


###################################################
### code chunk number 4: cost_clhs
###################################################
diamonds$cost <- runif(nrow(diamonds))
res_cost <- clhs(diamonds, size = 100, progress = FALSE, iter = 1000, cost = 'cost')


###################################################
### code chunk number 5: plot_clhs_1
###################################################
res <- clhs(diamonds, size = 100, simple = FALSE, progress = FALSE, iter = 1000)
plot(res)


###################################################
### code chunk number 6: plot_clhs_2
###################################################
plot(res, c('obj', 'box'))


###################################################
### code chunk number 7: plot_clhs_3
###################################################
res_cost <- clhs(diamonds, size = 100, progress = FALSE, iter = 1000, cost = 'cost', simple = FALSE)
plot(res_cost, c('obj', 'cost'))


###################################################
### code chunk number 8: plot_clhs_4
###################################################
plot(res_cost, c('obj', 'cost', 'box'))


