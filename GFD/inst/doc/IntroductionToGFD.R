## ------------------------------------------------------------------------
library(GFD)
data(pizza)

## ------------------------------------------------------------------------
head(pizza)

## ------------------------------------------------------------------------
set.seed(1234)
model1 <- GFD(Delivery ~ Crust * Coke * Bread, data = pizza, nperm = 10000, alpha = 0.05)
summary(model1)

## ------------------------------------------------------------------------
data("curdies")
set.seed(987)
nested <- GFD(dugesia ~ season + season:site, data = curdies)
summary(nested)

## ------------------------------------------------------------------------
plot(model1, factor = "Crust:Coke:Bread", legendpos = "center", main = "Delivery time of pizza", xlab = "Bread")
plot(model1, factor = "Crust:Coke", legendpos = "topleft", main = "Two-way interaction", xlab = "Coke", col = 3:5, pch = 17)
plot(nested, factor = "season:site", xlab = "site")

## ------------------------------------------------------------------------
calculateGUI()

