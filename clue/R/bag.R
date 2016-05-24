cl_bag <-
function(x, B, k = NULL,
         algorithm = "kmeans", parameters = NULL,
         method = "DFBC1", control = NULL)
{
    ## Currently, method 'DFBC1' (Dudoit-Fridlyand BagClust1) is the
    ## only one available, and argument 'control' is ignored.
    
    ## Construct reference partition.
    algorithm <- match.fun(algorithm)
    reference <- eval(as.call(c(list(algorithm, x),
                                if(!is.null(k)) list(k),
                                parameters)))

    ## Construct bootstrap ensemble.
    clusterings <-
        cl_boot(x, B, k, algorithm, parameters, resample = TRUE)

    ## Construct Dudoit-Fridlyand BagClust1 consensus partitions,
    ## suitably generalized ...
    ## <NOTE>
    ## In principle, this could be turned into a "constructive" method
    ## for cl_consensus(), also allowing for weights (straightforward).
    ## E.g.,
    ##    .cl_consensus_partition_DFBC1(clusterings, weights, control)
    ## where either 'control specifies a reference partition, or the
    ## first element of 'clusterings' is taken as such.
    ## </NOTE>
    k <- max(sapply(c(clusterings, reference), n_of_classes))
    M_ref <- cl_membership(reference, k)
    M <- matrix(0, NROW(M_ref), k)
    for(b in seq_len(B)) {
        mem <- cl_membership(clusterings[[b]], k)
        ## Match classes to reference partition.
        ind <- solve_LSAP(crossprod(M_ref, mem), maximum = TRUE)
        M <- M + mem[, ind]
    }
    as.cl_partition(cl_membership(as.cl_membership(M / B), k))
}
