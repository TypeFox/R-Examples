iinc <-
function (inc, PO, prsep = ", ", print.eqs = FALSE) 
{
    pi <- (transf(inc, "listmat", ord = dim(PO)[1], prsep = prsep) + 
        PO)
    ls <- list()
    k <- 1
    for (i in 1:(nrow(pi) - 1)) {
        for (j in (i + 1):nrow(pi)) {
            if (isTRUE(all.equal(pi[i, ], pi[j, ])) == TRUE) {
                if (i != j) {
                  ls[[k]] <- c(i, j)
                  k <- k + 1
                  if (isTRUE(print.eqs == TRUE) == TRUE) 
                    print(paste(i, j, sep = " = "))
                }
            }
        }
        rm(j)
    }
    rm(i)
    if (length(ls) == 0) 
        return(1:dim(PO)[1])
    if (sum(transf(inc, "listmat", ord = dim(PO)[1]), PO) == 
        dim(PO)[1] * dim(PO)[2]) 
        return(rep(1, each = dim(PO)[1]))
    j <- 1
    attr(ls[[1]], "names")[1] <- attr(ls[[1]], "names")[2] <- j
    if (length(ls) > 1) {
        for (i in (j + 1):length(ls)) {
            if (any(ls[[i]] %in% ls[[j]]) | any(ls[[j]] %in% 
                ls[[i]])) 
                attr(ls[[i]], "names")[1] <- attr(ls[[i]], "names")[2] <- j
        }
        rm(i)
        for (j in 2:length(ls)) {
            if (isTRUE(is.null(attr(ls[[j]], "names")) == TRUE) == 
                TRUE) {
                for (i in 1:(j - 1)) {
                  ifelse((any(ls[[i]] %in% ls[[j]]) | any(ls[[j]] %in% 
                    ls[[i]])), attr(ls[[j]], "names")[1] <- attr(ls[[j]], 
                    "names")[2] <- as.integer(attr(ls[[i]], "names")[1]), 
                    NA)
                }
                rm(i)
            }
        }
        rm(j)
        for (j in 2:length(ls)) if (isTRUE(is.null(attr(ls[[j]], 
            "names")) == TRUE) == TRUE) {
            for (i in j:length(ls)) {
                if (any(ls[[i]] %in% ls[[j]])) 
                  attr(ls[[i]], "names")[1] <- attr(ls[[i]], 
                    "names")[2] <- j
            }
            rm(i)
            if (j < length(ls)) {
                for (k in (j + 1):length(ls)) {
                  if (is.null(attr(ls[[k]], "names")) == TRUE) {
                    for (i in 1:length(ls)) {
                      if (is.null(attr(ls[[i]], "names")) == 
                        FALSE) {
                        if (any(ls[[k]] %in% ls[[i]])) 
                          attr(ls[[k]], "names")[1] <- attr(ls[[k]], 
                            "names")[2] <- as.numeric(attr(ls[[i]], 
                            "names")[1])
                      }
                    }
                    rm(i)
                  }
                }
                rm(k)
            }
        }
        rm(j)
    }
    clu <- vector()
    for (i in 1:length(ls)) {
        clu[i] <- as.numeric(attr(ls[[i]], "names")[1])
    }
    rm(i)
    f <- nlevels(factor(clu))
    x <- clu
    cls <- vector()
    length(cls) <- length(clu)
    for (i in 1:f) {
        cls[which(x == as.numeric(levels(factor(clu)))[i])] <- i
    }
    rm(i)
    nls <- ls
    for (i in 1:length(nls)) {
        attr(nls[[i]], "names")[1] <- attr(nls[[i]], "names")[2] <- cls[i]
    }
    rm(i)
    cl <- integer()
    length(cl) <- dim(PO)[1]
    for (k in 1:f) {
        for (i in 1:length(nls)) {
            if (attr(nls[[i]], "names")[1] == as.character(k)) {
                cl[nls[[i]][1]] <- k
                cl[nls[[i]][2]] <- k
            }
        }
        rm(i)
    }
    rm(k)
    n <- as.numeric(attr(stats::na.omit(cl), "na.action"))
    for (i in 1:length(n)) {
        cl[n[i]] <- f + i
    }
    rm(i)
    rm(f)
    return(cl)
}
