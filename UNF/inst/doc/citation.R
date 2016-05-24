## ------------------------------------------------------------------------
library("UNF")

## ------------------------------------------------------------------------
library("datasets")
library("digest")
data(iris)

## ------------------------------------------------------------------------
digest(iris, "md5")
digest(head(iris), "md5")

## ------------------------------------------------------------------------
digest(iris, "md5")
digest(rev(iris), "md5")

## ------------------------------------------------------------------------
iris2 <- iris
names(iris2) <- gsub("\\.","",names(iris2))
names(iris)
names(iris2)
digest(iris, "md5")
digest(iris2, "md5")

## ------------------------------------------------------------------------
library("tools")
save(iris, file = "iris.RData")
write.csv(iris, file = "iris.csv", row.names = FALSE)
md5sum("iris.RData")
md5sum("iris.csv")

## ------------------------------------------------------------------------
# read CSV back into memory
d <- read.csv("iris.csv")
identical(iris, d)
digest(iris, "md5")
digest(d, "md5")

## ----, echo = FALSE------------------------------------------------------
# cleanup files
unlink("iris.RData")
unlink("iris.csv")

## ------------------------------------------------------------------------
x1 <- 1:20
x2 <- x1 + 1e-7
identical(digest(x1), digest(x2))

## ------------------------------------------------------------------------
write.csv(iris, file = "iris.csv", row.names = FALSE)
iris2 <- read.csv("iris.csv")
identical(iris, iris2)
identical(unf(iris), unf(iris2))

## ----, echo = FALSE------------------------------------------------------
# cleanup files
unlink("iris.csv")

## ------------------------------------------------------------------------
x1 <- 1:20
x2 <- x1 + 1e-7
identical(unf(x1), unf(x2))
x3 <- x1 + 1e-3
identical(unf(x1), unf(x3))
identical(unf(x1, digits = 3), unf(x2, digits = 3))

## ------------------------------------------------------------------------
x1 <- 1:20
x2 <- as.character(x1)
identical(unf(x1), unf(x2))

## ------------------------------------------------------------------------
unf(iris)

## ----, echo=FALSE--------------------------------------------------------
cat("> Anderson, Edgar (1935). The irises of the Gaspe Peninsula, *Bulletin of the American Iris Society*, 59, 2-5.", unf(iris)$formatted, "\n")

## ------------------------------------------------------------------------
unf(iris)

## ------------------------------------------------------------------------
unf(iris[iris$Species == 'setosa',])
unf(iris[iris$Species == 'versicolor',])
unf(iris[iris$Species == 'virginica',])

## ------------------------------------------------------------------------
set.seed(123)
unf(iris[sample(20, 1:nrow(iris), FALSE),])

