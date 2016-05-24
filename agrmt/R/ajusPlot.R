ajusPlot <- function(V, tolerance=0.1, variant="modified", ...) {
# plot AJUS with additional information
# additional arguments can now be passed on (e.g col="red", lwd=2)
# removing NA to deal with real-life data
# (e.g. time series with missing data)
a <- ajus(na.omit(V), tolerance=tolerance, variant=variant) # do not consider missing values
l <- length(V)          # for plotting
m <- max(V, na.rm=TRUE) # for plotting
i <- min(V, na.rm=TRUE) # for difference d
d <- m - i              # for plotting
#plot(V, type="b", axes=FALSE, xlab="", ylab="")
miss <- !is.na(V)   # for which years do data exist?
plot(which(miss), na.omit(V), type="b", axes=FALSE, xlab="", ylab="", xlim=c(1,l), ...) # plot where data exist
axis(1, at=1:l) # control axes
axis(2)
if (m==0) m <- 1 # avoid plotting type on line if "F" at 0
if (d==0) d <- m # avoit plotting type on line if "F"
text(1+l/5, m-d/5, paste("type", a$type, sep=" "), cex=2)
text(1+l/5, m-d/3, paste("skew", a$skew, sep=" "))
}
