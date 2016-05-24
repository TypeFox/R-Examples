## ---- include=FALSE------------------------------------------------------
library(exploreR)


## ------------------------------------------------------------------------

reset()


## ------------------------------------------------------------------------

sampleData <- iris

head(sampleData)


## ------------------------------------------------------------------------
regressResults <- masslm(sampleData, "Sepal.Length", ignore = "Species")

regressResults


## ------------------------------------------------------------------------
massregplot(sampleData, "Sepal.Length", ignore = "Species")


## ------------------------------------------------------------------------
stand.Petals <- standardize(sampleData, c("Petal.Width", "Petal.Length"))

head(stand.Petals)


