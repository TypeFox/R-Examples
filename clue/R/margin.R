cl_margin <-
function(x)
{
    if(is.cl_hard_partition(x))
        out <- rep.int(1, n_of_objects(x))
    else if(is.cl_partition(x)) {
        x <- cl_membership(x)
        i <- seq_len(nrow(x))
        j <- cbind(i, max.col(x))
        out <- x[j]
        x[j] <- 0
        out <- out - x[cbind(i, max.col(x))]
    }
    else
        stop("Argument 'x' must be a partition.")

    out
}
