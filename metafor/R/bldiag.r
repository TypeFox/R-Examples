bldiag <-
function (...) 
{
    mlist <- list(...)
    if (length(mlist) == 1) 
        mlist <- unlist(mlist, recursive = FALSE)
    csdim <- rbind(c(0, 0), apply(sapply(mlist, dim), 1, cumsum))
    out <- array(0, dim = csdim[length(mlist) + 1, ])
    add1 <- matrix(rep(1:0, 2), ncol = 2)
    for (i in seq(along = mlist)) {
        indx <- apply(csdim[i:(i + 1), ] + add1, 2, function(x) x[1]:x[2])
        if (is.null(dim(indx))) {
            out[indx[[1]], indx[[2]]] <- mlist[[i]]
        }
        else {
            out[indx[, 1], indx[, 2]] <- mlist[[i]]
        }
    }
    return(out)
}
