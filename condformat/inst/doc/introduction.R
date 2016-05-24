## ------------------------------------------------------------------------
data(iris)
library(condformat)
condformat(iris[c(1:5,70:75, 120:125),]) +
  rule_fill_discrete(Species) + 
  rule_fill_discrete(Sepal.Width, Sepal.Length,
                     expression = Sepal.Width > Sepal.Length - 2.25,
                     colours = c("TRUE" = "#7D00FF")) + 
  rule_fill_gradient2(Petal.Length)

