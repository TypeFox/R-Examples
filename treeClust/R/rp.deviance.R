rp.deviance <- function (x, ...) 
{
#
# Compute deviance at every node of an rpart tree. It has already
# been computed for regression trees, but for classification trees
# we have to compute it ourselves.
#
if (x$numresp == 1) # regression tree
    return (x$frame$dev)
#
# If the response for this tree had k levels, the yval2 element of
# the frame has 2k + 2 columns. We drop the first (plurality class label)
# and last (node probability). Among those that remain, the first k
# are within-node class counts and the second k, within-node class probs.
#
yv <- x$frame$yval2
k <- (ncol(yv) - 2)/2
N <- yv[,2:(k+1)] # upper-case since debugger uses "n"
p <- yv[,(k+2):(2 * k + 1)]
dev.mat <- -2 * N * log (p)
dev.mat[is.na (dev.mat)] <- 0
return (rowSums (dev.mat))
}

