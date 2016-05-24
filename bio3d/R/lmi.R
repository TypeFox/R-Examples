lmi <- function (trj, grpby = NULL, ncore=1) {
    ncore <- setup.ncore(ncore)

# rm:r-value matrix
    cm <- var(trj)

# mclapply or lapply
    if (ncore > 1) {
        lmiapply = mclapply
    } else {
        lmiapply = lapply
    }  
    
    rm <- cov2dccm(cm, method = "lmi", ncore = ncore)

# group by or not
    if (!is.null(grpby)) {
        if (ncol(trj) != (length(grpby) * 3)) 
            stop("dimension miss-match in 'trj' and 'grpby', check lengths")
        inds <- bounds(grpby, dup.inds = TRUE)
        l <- dim(inds)[1]
        m <- matrix(, ncol = l, nrow = l)
        ij <- pairwise(l)
# list3: lmi 
        list3 <- lmiapply(1:nrow(ij), function(k) max(rm[(inds[ij[k, 1], "start"]:inds[ij[k, 1], "end"]), (inds[ij[k, 2], "start"]:inds[ij[k, 2], "end"])], na.rm = TRUE))
        list3 <- unlist(list3)

        for (k in 1:nrow(ij)) {
            m[ij[k, 1], ij[k, 2]] <- list3[k]
        }
        m[lower.tri(m)] = t(m)[lower.tri(m)]
        diag(m) <- 1; rm=m
    }
    class(rm) = c("dccm", "matrix")
    return(rm)
}


