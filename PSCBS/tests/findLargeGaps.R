library("PSCBS")

# Simulating copy-number data
set.seed(0xBEEF)

# Simulate CN data
J <- 1000
mu <- double(J)
mu[200:300] <- mu[200:300] + 1
mu[350:400] <- NA # centromere
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- seq(from=1, to=100e6, length.out=J)

data <- data.frame(chromosome=0L, x=x)

gaps <- findLargeGaps(x=x, minLength=1e6)
print(gaps)
stopifnot(is.data.frame(gaps))
stopifnot(nrow(gaps) == 0L)
segs <- gapsToSegments(gaps)
print(segs)
stopifnot(is.data.frame(segs))
stopifnot(nrow(segs) == 1L)


gaps <- findLargeGaps(data, minLength=1e6)
print(gaps)
stopifnot(is.data.frame(gaps))
stopifnot(nrow(gaps) == 0L)
segs <- gapsToSegments(gaps)
print(segs)
stopifnot(is.data.frame(segs))
stopifnot(nrow(segs) == 1L)


## Add missing values
data2 <- data
data$x[30e6 < x & x < 50e6] <- NA
gaps <- findLargeGaps(data, minLength=1e6)
print(gaps)
stopifnot(is.data.frame(gaps))
stopifnot(nrow(gaps) == 1L)
segs <- gapsToSegments(gaps)
print(segs)
stopifnot(is.data.frame(segs))
stopifnot(nrow(segs) == 3L)



# BUG FIX: Issue #6
gaps <- findLargeGaps(chromosome=rep(1,10), x=1:10, minLength=2)
print(gaps)
stopifnot(is.data.frame(gaps))
stopifnot(nrow(gaps) == 0L)
# BUG FIX: Issue #9
segs <- gapsToSegments(gaps)
print(segs)
stopifnot(is.data.frame(segs))
stopifnot(nrow(segs) == 1L)


# BUG FIX: PSCBS GitHub Issue #8
gaps <- try({
  findLargeGaps(chromosome=rep(1,3), x=as.numeric(1:3), minLength=1)
})
stopifnot(inherits(gaps, "try-error"))
