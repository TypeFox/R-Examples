"plot.knncat" <- 
function (x, ...) 
{
#
# Plot.knncat: plot method for knncat objects
#
# Arguments: x: knncat object, from knncat().
#
y <- unlist (x$phi)
maxx <- length(x$phi)
#
# Set up plot without axes; add y axix
#
ten.pct <- maxx * .10
plot (c(1-ten.pct, maxx+ten.pct), range(y), xlab = "Variable", ylab = "Phi", 
      type = "n", axes=FALSE)
box ()
axis (2)
#
# Add x axis with proper labels
#
axis (1, at = 1:maxx, labels=names(x$vars))
#
# Now plot the phi's, using the levels names.
#
for (i in 1:length(x$vars))
{
    vec <- x$phi[[i]]
    text (i, vec, names(vec))
}
}
