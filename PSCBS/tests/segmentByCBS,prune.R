library("PSCBS")

## Compare segments
assertMatchingSegments <- function(fitM, fit) {
  chrs <- getChromosomes(fitM)
  segsM <- lapply(chrs, FUN=function(chr) {
    getSegments(extractChromosome(fitM, chr))
  })
  segs <- lapply(fit[chrs], FUN=getSegments)
  stopifnot(all.equal(segsM, segs, check.attributes=FALSE))
}

## Simulate data
set.seed(0xBEEF)
J <- 1000
mu <- double(J)
mu[200:300] <- mu[200:300] + 1
mu[350:400] <- NA
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- sort(runif(length(y), max=length(y))) * 1e5

data <- list()
for (chr in 1:2) {
  data[[chr]] <- data.frame(chromosome=chr, x=x, y=y)
}
data$M <- Reduce(rbind, data)

## Segment
message("*** segmentByCBS()")
fit <- lapply(data, FUN=segmentByCBS)
print(fit)
assertMatchingSegments(fit$M, fit)

## Join segments
message("*** joinSegments()")
fitj <- lapply(fit, FUN=joinSegments)
print(fitj)
assertMatchingSegments(fitj$M, fitj)

## Reset segments
message("*** resetSegments()")
fitj <- lapply(fit, FUN=resetSegments)
print(fitj)
assertMatchingSegments(fitj$M, fitj)

## Prune by SD undo
message("*** pruneBySdUndo()")
fitp <- lapply(fit, FUN=pruneBySdUndo)
print(fitp)
assertMatchingSegments(fitp$M, fitp)

## Prune by hierarchical clustering
message("*** pruneByHClust()")
fitp <- lapply(fit, FUN=pruneByHClust, k=1L)
print(fitp)
assertMatchingSegments(fitp$M, fitp)
