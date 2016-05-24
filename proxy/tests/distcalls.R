#################################
## Test for dist calls
##################################

library(proxy)

set.seed(20140107)

## get all measures
proxies = pr_DB$get_entry_names()

## remove special cases
proxies = setdiff(proxies, c("Mahalanobis", "Minkowski", "Stiles", "Levenshtein", "fJaccard"))

## create test data
x = matrix(1:100, 10)

## test function: checks if dist(x) == dist(x,x) for all measures,
## and if diag(dist(x, x)) == diag(x, x, pairwise = TRUE)
prtest <- function(...) {
    CD <- dist(x, x, ...)
    all(as.matrix(dist(x, ...)) == CD) &&
      all(diag(CD) == dist(x, x, pairwise = TRUE, ...))
}

## loop over all measures (except special cases)
for (i in proxies)
    {cat(i); prtest(i); cat(": OK.\n")}

## Minkowski
writeLines("Minkowski:")
for (j in c(0, 0.5, 1, 2, 3, Inf))
    {cat("p =", j); prtest(i = "Minkowski", p = j); cat(": OK.\n")}

## Mahalanobis (need non-singular matrix)
x = as.matrix(iris[1:50,-5])
prtest("Mahalanobis")

## fJaccard (needs values in unit interval)
x = as.matrix(1:100/100, 10)
prtest("fJaccard")

## produce binary matrix
x = matrix(rbinom(100,1,0.7), 10)

## Stiles (gives a lot of warnings due to log)
tmp = dist(x, "Stiles")
tmp = dist(x, x, "Stiles")

## try again (almost) all measures, this time with binary data to check
## conversions
for (i in proxies)
    {cat(i); prtest(i); cat(": OK.\n")}
for (j in c(0, 0.5, 1, 2, 3, Inf))
    {cat("p =", j); prtest(i = "Minkowski", p = j); cat(": OK.\n")}

## Levenshtein distance
s <- c("A", "quick", "brown", "fox", "jumps", "over", "the", "lazy", "dog")
all(as.matrix(dist(s, "Levenshtein")) == dist(s, s, "Levenshtein"))

## Test auto-conversion
x = iris[,-5]
prtest()


