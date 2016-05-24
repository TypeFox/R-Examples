treeClust <- function(dfx, d.num = 1, col.range = 1:ncol (dfx), verbose = F, 
         final.algorithm, k, control = treeClust.control(), rcontrol = rpart.control (), ...)
{
if(!is.data.frame(dfx))
    stop("This function requires a data frame")

leaf.matrix <- as.data.frame(matrix(0., nrow(dfx), ncol(dfx)))
#
# Fail if any columns have names with embedded spaces.
#
dimnames(leaf.matrix) <- dimnames(dfx)
nm <- names(dfx)
if (length (grep (" ", nm) > 0))
    stop ("Some columns have embedded spaces in their names, and that's not good.")
#
# Final.algorithm has to be agnes, pam, clara, or kmeans (or missing). For everything except
# agnes, a k is required.
#
additional.args <- list (...)
if (!missing (final.algorithm)) {
    if (!is.element (final.algorithm, c("agnes", "pam", "clara", "kmeans")))
        stop ("Unrecognized final algorithm")
#
# For pam, clara or k-means, "k" must be present. For agnes, set it to -1 if 
# it's missing.
#
    if (missing (k))
        if (is.element (final.algorithm, c("kmeans", "pam", "clara")))
            stop ("With kmeans, pam or clara, specify the number of clusters 'k'")
        else k <- -1
        if (k == -1 & control$cluster.only == T) stop ("Cluster.only TRUE requires k")

}
#
# Set up the two-column "result" matrix.
#
results <- matrix(0., nrow = ncol(dfx), ncol = 2.)
dimnames(results) <- list(dimnames(dfx)[[2.]], c("DevRat", "Size"))
#
# Set up the big list of trees if asked -- except that we will need
# that list no matter what, when d.num == 4 or we're asked for "newdata."
# Save that onto "control" for later use.
#
control$keep.trees <- FALSE
if (control$return.trees || control$return.newdata || d.num == 4) {
    control$keep.trees <- TRUE
    big.list.of.trees <- vector("list", ncol(dfx))
}

#
# If the user asked for a cluster, try to create it. "Parproc" will be
# TRUE if we're processing in parallel.
#
parproc <- FALSE
if (is.numeric (control$parallelnodes) && control$parallelnodes > 1) {
    cl <- parallel::makePSOCKcluster (round (control$parallelnodes))
    on.exit (parallel::stopCluster (cl))
    parproc <- TRUE
}
#
# Loop over columns, possibly using the cluster. For each column, call treeClust.rpart() to build
# a tree. Fill the "results" matrix with the DevRat and Size entries. Keep the trees if (a) asked or
# (b) we want to compute "newdata" or (c) we want to compute d4.
#
if (parproc) {
    caout <- parallel::clusterApplyLB (cl, col.range, treeClust.rpart,
                                dfx = dfx, d.num = d.num, control = control,
                                rcontrol = rcontrol)
    keepers <- sapply (caout, function (x) x$Size > 1)
    leaf.matrix <- as.data.frame (sapply (caout[keepers], function (x) x$leaf.where))
    names (leaf.matrix) <- names (dfx)[keepers]
    results <- as.data.frame (t(sapply (caout, function (x) unlist (x[1:2]))))
    if (control$keep.trees)
        big.list.of.trees <- lapply (caout[keepers], function (x) x$tree)
} else {
    for(i in col.range) {
        if(verbose > 0)
            cat("Creating rpart tree with column", i, "\n")
        out.i <- treeClust.rpart(i, dfx, d.num, control, rcontrol)
        results[i,] <- c(out.i$DevRat, out.i$Size)
        if (results[i, "Size"] > 1) {
            leaf.matrix[,i] <- out.i$leaf.where
            if (control$keep.trees)
                big.list.of.trees[[i]] <- out.i$tree
        }
    }
    if (control$keep.trees)
        big.list.of.trees <- big.list.of.trees[results[,"Size"] > 1]
    leaf.matrix <- leaf.matrix[,results[,"Size"] > 1]
}

if(!any(results[, "Size"] > 1))
    stop("No tree produced anything! Panic!")
#
# Save original results. For the moment we will keep both the full set
# of results and the one that shows only trees that were kept. This
# is a little redundant, but we want (# row in result) to be equal to
# (# of trees).
#
original.results <- results
results <- results[results[,"Size"] > 1,, drop=F]
#
# If there's a final algorithm, prepare for it now. Agnes and pam use
# the inter-point distances, plus any additional arguments, so we'll
# compute those. Also compute them if we've asked for them explicitly.
#
if (control$return.dists == TRUE ||
    (!missing (final.algorithm) && is.element (final.algorithm, c("pam", "agnes"))))
{
# Compute dists if needed. d1 and d2 don't use the trees.
    dists <- tcdist (tbl = results, mat = leaf.matrix, 
            trees = big.list.of.trees, d.num = d.num)
}
#
# K-means and clara use "newdata". Compute that thing if asked.
#
if (control$return.newdata == TRUE ||
    (!missing (final.algorithm) && is.element (final.algorithm, c("kmeans", "clara"))))
{
    newdata <- tcnewdata (tbl = results, mat = leaf.matrix, 
                              trees = big.list.of.trees, d.num = d.num)
}
#
# Now, if there's a final algorithm, call it.
#
if (missing (final.algorithm)) {
    final.algorithm <- "None"
    final.clust <- NULL
} else {
#
# We call "agnes" or "pam" by "do.call", which saves a copy of the dists
# in the call element. That thing is huge and unnecessary, so we remove it.
#
    if (final.algorithm == "agnes") {
        final.clust <- do.call (final.algorithm, list (x = dists, ...))
        final.clust$call$x <- "deleted"
        if (control$cluster.only == TRUE)
            final.clust <- cutree (final.clust, k = k)
    }
    if (final.algorithm == "pam") {
        final.clust <- do.call (final.algorithm, list (x = dists, k = k, ...))  
        final.clust$call$x <- "deleted"
        if (control$cluster.only == TRUE)
            final.clust <- final.clust$clustering
    }
#
# k-means are clara use the "newdata" data, which has p-1 columns for each
# tree with p leaves.
#
    if (is.element (final.algorithm, c("clara", "kmeans"))) {
        if (final.algorithm == "kmeans")
            final.clust <- kmeans (x = newdata, centers = k)
        else 
            final.clust <- clara (x = newdata, k = k)
        if (control$cluster.only == TRUE)
                final.clust <- final.clust$cluster
    }
} # end of "final algorithm" stuff

#
# If we were only asked for the clustering, return that.
#
if (control$cluster.only)
    return (final.clust)
#
# Set up return value; add requested stuff.
#
return.val <- list(call = match.call(), d.num = d.num, tbl = results, extended.tbl = original.results,
                   final.algorithm = final.algorithm, 
                   final.clust = final.clust, additional.args = additional.args)
if (control$return.trees)
    return.val$trees <- big.list.of.trees
if (control$return.dists)
    return.val$dists <- dists
if (control$return.newdata)
    return.val$newdata <- newdata
if (control$return.mat)
    if (!missing (final.algorithm) && is.element (final.algorithm, c("clara", "kmeans")))
        return.val$mat <- newdata
    else
        return.val$mat <- leaf.matrix
class(return.val) <- "treeClust"

return(return.val)
}

