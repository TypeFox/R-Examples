hierar <-
function (W, x, type = c("person", "relation")) 
{
    if (isTRUE(attr(W, "class") == "Rel.Box") == FALSE) 
        stop("\"W\" must be a \"Rel.Box\" class.")
    if (isTRUE(is.null(W$lbs) == TRUE) == FALSE) {
        if (isTRUE(is.numeric(x) == TRUE) == TRUE) {
            if (isTRUE(length(W$lbs) < x) == TRUE) 
                stop("\"x\" exceeds the order of the specified network")
        }
        if (isTRUE(is.character(x) == TRUE) == TRUE) {
            if (isTRUE(x %in% W$lbs) == FALSE) 
                stop("\"x\" is not an actor in the specified network")
        }
    }
    else if (isTRUE(is.null(W$lbs) == TRUE) == TRUE) {
        if (isTRUE(dim(W$W)[1] < x) == TRUE) 
            stop("\"x\" exceeds the order of the specified network")
    }
    ifelse(isTRUE(is.numeric(x) == TRUE), X <- x, X <- which(W$lbs == 
        x))
    rele <- unique(t(W$W[X, , ]))
    switch(match.arg(type), person = {
        Ph <- as.data.frame(array(0, dim = c(nrow(t(rele)), nrow(t(rele)))))
    }, relation = {
        Rh <- as.data.frame(array(0, dim = c(nrow(rele), nrow(rele))))
    })
    switch(match.arg(type), person = {
        for (j in 1:nrow(t(rele))) {
            for (i in 1:nrow(t(rele))) {
                if ((as.numeric(any(rele[, i] < rele[, j])) == 
                  1 && as.numeric(any(rele[, j] < rele[, i])) == 
                  0) | as.numeric(all(rele[, i] == rele[, j])) == 
                  1) Ph[i, j] <- 1
            }
        }
        rm(j)
        for (i in 1:nrow(t(rele))) {
            if (sum(rele[, i]) == 0) Ph[i, ] <- 0
        }
        rm(i)
    }, relation = {
        for (j in 1:nrow(rele)) {
            for (i in 1:nrow(rele)) {
                if ((as.numeric(any(rele[i, ] < rele[j, ])) == 
                  1 && as.numeric(any(rele[j, ] < rele[i, ])) == 
                  0) | as.numeric(all(rele[i, ] == rele[j, ])) == 
                  1) Rh[i, j] <- 1
            }
        }
        rm(j)
    })
    switch(match.arg(type), person = {
        ph <- as.matrix(Ph)
    }, relation = {
        rh <- as.matrix(Rh)
    })
    switch(match.arg(type), person = {
        if (isTRUE(is.null(W$lbs) == FALSE)) dimnames(ph)[[1]] <- dimnames(ph)[[2]] <- W$lbs
        ph <- list(ph = ph, pers = x)
        return(ph)
    }, relation = {
        if (isTRUE(is.null(dimnames(W$W)[[3]]) == FALSE)) attributes(rh)$dimnames[[2]] <- attributes(rh)$dimnames[[1]] <- dimnames(rele)[[1]]
        rh <- list(rh = rh, pers = x)
        return(rh)
    })
}
