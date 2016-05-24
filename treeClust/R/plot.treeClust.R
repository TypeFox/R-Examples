plot.treeClust <- function (x, extended = TRUE, ...) 
{
#
# Plot method for objects of class treeClust.
#        x: Object of class treeClust
# extended: if TRUE, include variables whose trees were dropped
#      ...: Other argument to plot()
#
if (extended)
    tbl <- x$extended.tbl
else
    tbl <- x$tbl
dev.rat <- tbl[,"DevRat"] / max (tbl[,"DevRat"])
plot (dev.rat, type = "n", ylim = c(0, 1),
 xlim = c(1,length (dev.rat)), xlab = "Variable Number", ylab = "Scaled Deviance Ratio",
  ...)

text (1:length (dev.rat), dev.rat, tbl[,"Size"])
}

