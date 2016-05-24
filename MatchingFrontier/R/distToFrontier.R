# On 1.18.2015, fixed serious bug that
# was calculating the wrong Y values for
# the frontier, occasionally making it appear
# non-monotonic.

distToFrontier <-
function(distance.mat){
    N <- sum(dim(distance.mat))
    cat("Calculating theoretical frontier...\n")

    # Get matches for T and C
    row.mins <- apply(distance.mat, 1, function(x) min(x))
    col.mins <- apply(distance.mat, 2, function(x) min(x))
    
    row.mins.inds <- apply(distance.mat, 1, function(x) as.integer(names(which.min(x))))
    col.mins.inds <- apply(distance.mat, 2, function(x) as.integer(names(which.min(x))))
    
    matched.to <- c(row.mins.inds, col.mins.inds)[order(as.integer(names(c(row.mins.inds, col.mins.inds))))]
    
    minimums <- c(row.mins, col.mins)
    
    sorted.minimums <- sort(unique(c(row.mins, col.mins)), decreasing = TRUE)
    drop.order <- lapply(sorted.minimums, function(x) as.integer(names(minimums[minimums == x])))
    
    cat("Calculating information for plotting the frontier...\n")
    weighted.vals <- unlist(lapply(drop.order, function(x) length(x))) * sorted.minimums
    Xs <- c(0, cumsum(lapply(drop.order, function(x) length(x))))
    Xs <- Xs[1:(length(Xs) - 1)]

    Ys <- rev(cumsum(rev(weighted.vals))) / (N - Xs)

    # Checks to confirm monotonically decreasing. Since
    # that's theoretically impossible, if the condition is
    # met, there's a serious bug somewhere in the code. 
    if(any(diff(Ys) > 0 )){
        stop('Something is very wrong. Email clucas@fas.harvard.edu.')
    }
    return(list(drop.order = drop.order, Xs = Xs, Ys = Ys, matched.to = matched.to, distance.mat = distance.mat))
}

