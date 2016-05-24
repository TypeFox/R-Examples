## load library
require("GMD")

## load data
data(iris)

## create common bins
n <- 30                                    # the number of bins
breaks <- gbreaks(iris[,"Sepal.Length"],n) # the boundary of bins 

## create a list of histograms
Sepal.Length <-
  list(setosa=ghist(iris[iris$Species=="setosa","Sepal.Length"],breaks=breaks),
       versicolor=ghist(iris[iris$Species=="versicolor","Sepal.Length"],breaks=breaks),
       virginica=ghist(iris[iris$Species=="virginica","Sepal.Length"],breaks=breaks)
       )

## convert to a `hist' object
x <- as.mhist(Sepal.Length)

## get bin-wise summary statistics
summary(x)

## visualize the histograms
plot(x,beside=FALSE,
     main="Histogram of Sepal.Length of iris",xlab="Sepal.Length (cm)")
