preprocess.data <- function (data, lower, upper) 
{
#
# Preprocess data to produce epiparameters.
# data: numeric vector of data
# lower: numeric lower bound, optional
# upper: numeric upper bound, optional
#
out <- list ()

mean.x <- mean (data)
sd.x <- sd(data)
data <- sort (data)
#
# Lower bound: if NULL, use min(data) - 2 sd; if -Inf, use
# mean (data) - 10 sd; otherwise use specified value
if (missing (lower) || is.null (lower))
    lower <- out$m0 <- data[1] - 2 * sd.x
else if (lower == -Inf)
    lower <- out$m0 <- mean.x - 10 * sd.x
else
    out$m0 <- lower

# Upper bound, likewise, with +2 or +10 sds.
if (missing (upper) || is.null (upper))
    upper <- out$mN <- data[length(data)] + 2 * sd.x
else if (upper == Inf)
    upper <- out$mN <- mean.x + 10 * sd.x
else
    out$mN <- upper
#
# Keep only data within these ranges
#
out$data <- data[data >= lower & data <= upper]
return (out)
}
