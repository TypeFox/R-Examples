# This example shows how to do resampling using the "times"
# convenience function.  I then show a slightly fancier version
# using foreach with the ".final" argument to perform the
# post processing.
#
# Thanks to Daniel Kaplan for this basic method of resampling, and for
# inspiring the "times" function, as well as the foreach package itself.

library(foreach)
library(iterators)

# Define a simple function for resampling a data frame
resample <- function(d, n=nrow(d)) d[sample(nrow(d), n, replace=TRUE),]

# Set the seed so we can reproduce our results
set.seed(123)

# Resample and compute the mean a thousand times
x <- times(1000) %do% with(resample(quakes), mean(mag))

# Print the standard deviation of the mean values
print(sd(x))

# Now let's do the same thing with foreach and use the .final argument
# setting the seed the same as before
set.seed(123)
xsd <- foreach(icount(1000), .combine='c', .final=sd) %do%
  with(resample(quakes), mean(mag))

# Print the standard deviation of the mean values returned by foreach/%do%
print(xsd)
