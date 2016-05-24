# This example takes the sum of the square roots of the numbers from
# one to a billion.  Note that this can't be done with the expression:
#
#   sum(sqrt(1:1000000000))
#
# since that attempts to allocate a 3.7 Gb vector, which generates
# an error in versions of R including 2.10.1.

library(itertools)
library(foreach)

# Size of input vector
n <- 1000000000

# The best value for chunkSize depends on how much memory you have.
# Generally, if chunkSize is too small, the code becomes inefficient.
# If chunkSize is too big, you either run out of memory or suffer from
# virtual memory thrashing.  I think that a value of 5 million should
# avoid memory thrashing on most modern computers without being too
# inefficient, but if you have 2 Gigabytes of memory or more, you
# might want to increase this value.
chunkSize <- 5000000

r <- foreach(x=isplitIndices(n, chunkSize=chunkSize), .combine='sum') %do% {
  sqrt(x)
}

cat(sprintf('sum(sqrt(1:%d)) = %e\n', n, r))
