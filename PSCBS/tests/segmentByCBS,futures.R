library("PSCBS")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xBEEF)

# Number of loci
J <- 1000

mu <- double(J)
mu[200:300] <- mu[200:300] + 1
mu[350:400] <- NA # centromere
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- sort(runif(length(y), max=length(y))) * 1e5
w <- runif(J)
w[650:800] <- 0.001

## Create multiple chromosomes
data <- knownSegments <- list()
for (cc in 1:3) {
  data[[cc]] <- data.frame(chromosome=cc, y=y, x=x)
  knownSegments[[cc]] <- data.frame(
    chromosome=c(   cc,  cc,  cc),
    start     =x[c(  1, 350, 401)],
    end       =x[c(349, 400,   J)]
  )
}
data <- Reduce(rbind, data)
str(data)
knownSegments <- Reduce(rbind, knownSegments)
str(knownSegments)

message("*** segmentByCBS() via futures ...")

message("*** segmentByCBS() via 'lazy' futures without attaching 'future' ...")
future::plan("lazy")
print(future::plan)
fitL <- segmentByCBS(data, seed=0xBEEF, verbose=TRUE)
print(fitL)


message("*** segmentByCBS() via futures with 'future' attached ...")
library("future")
oplan <- plan()

strategies <- c("eager", "lazy")
if (supportsMulticore()) strategies <- c(strategies, "multicore")

## Test 'async' futures?
pkg <- "async"
if (require(pkg, character.only=TRUE)) {
  backend("local")
  strategies <- c(strategies, "batchjobs")
}

message("Future strategies to test: ", paste(sQuote(strategies), collapse=", "))

fits <- list()
for (strategy in strategies) {
  message(sprintf("- segmentByCBS() using '%s' futures ...", strategy))
  plan(strategy)
  fit <- segmentByCBS(data, seed=0xBEEF, verbose=TRUE)
  fits[[strategy]] <- fit
  stopifnot(all.equal(fit, fits[[1]]))
  stopifnot(all.equal(fit, fitL))
}


message("*** segmentByCBS() via futures with known segments ...")
fits <- list()
dataT <- subset(data, chromosome == 1)
for (strategy in strategies) {
  message(sprintf("- segmentByCBS() w/ known segments using '%s' futures ...", strategy))
  plan(strategy)
  fit <- segmentByCBS(dataT, knownSegments=knownSegments, seed=0xBEEF, verbose=TRUE)
  fits[[strategy]] <- fit
  stopifnot(all.equal(fit, fits[[1]]))
}

message("*** segmentByCBS() via futures ... DONE")


## Cleanup
plan(oplan)
rm(list=c("fits", "dataT", "data", "fit"))

