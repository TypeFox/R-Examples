# This function uses a resampling approach for computing
# the probability that at least two people will have
# the same birthday in a group of "n" people.
pbday <- function(n) {
  ntests <- 1000
  pop <- 1:365
  anydup <- function(i) any(duplicated(sample(pop, n, replace=TRUE)))
  sum(sapply(seq(ntests), anydup)) / ntests
}

# Execute the pbday function for group sizes 1 to 100
x <- 1:100

# Sequential version
if (FALSE) {
  prob <- sapply(x, pbday)
}

# Parallel version
if (TRUE) {
  library(nws)
  s <- sleigh()
  prob <- unlist(eachElem(s, pbday, x))
}

# Display a plot of the probability vs. group size
plot(x, prob, main='Birthday Paradox', xlab='Group Size', ylab='Probability')

# Compare with the analytical results, using the
# standard R "pbirthday" function.
# The results are only expected to be close.
print(all.equal(prob, sapply(x, pbirthday)))
