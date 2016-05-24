reducs <-
function (x, clu, labels = NULL) 
{
    lngt <- nlevels(factor(clu))
    cls <- list()
    for (i in 1:lngt) {
        cls[[i]] <- which(clu == i)
    }
    rm(i)
    if (isTRUE(attr(x, "class")[1] == "Semigroup") == FALSE) {
        stop("\"x\" should be an object of a \"Semigroup\" class.")
    }
    tmp <- x
    lb <- x$st
    ifelse(isTRUE(attr(tmp, "class")[2] == "symbolic") == TRUE, 
        x <- convert(tmp), x <- tmp$S)
    ifelse(isTRUE(is.array(x)) == TRUE, xa <- x, xa <- array(x, 
        dimnames = lb))
    px <- perm(xa, clu, rev = FALSE)
    tab <- tabulate(clu)
    bm <- array(dim = c(lngt, lngt))
    y <- h <- j <- 0
    for (i in 1:lngt) {
        k <- ((y + 1):(y + tabulate(clu)[i]))
        for (j in 1:lngt) {
            for (q in 1:length(cls)) {
                if (all(unique(as.numeric(levels(factor(as.matrix(x[which(clu == 
                  i), which(clu == j)]))))) %in% cls[[q]]) == 
                  TRUE) {
                  if (isTRUE(attr(tmp, "class")[2] == "symbolic") == 
                    TRUE) {
                    bm[i, j] <- lb[min(as.numeric(as.matrix(px[k, 
                      (h + 1):(h + tabulate(clu)[j])])))]
                  }
                  else if (isTRUE(attr(tmp, "class")[2] == "numerical") == 
                    TRUE) {
                    bm[i, j] <- min(as.numeric(as.matrix(px[k, 
                      (h + 1):(h + tabulate(clu)[j])])))
                  }
                }
            }
            if (isTRUE(j > 0) == TRUE) 
                h <- sum(tab[1:j])
        }
        rm(j)
        h <- 0
        y <- sum(tab[1:i])
    }
    rm(i)
    if (is.null(labels) == FALSE) {
        lbs <- labels
    }
    else if (is.null(labels) == TRUE) {
        lbs <- vector()
        for (i in 1:length(tabulate(clu))) {
            ifelse(isTRUE(attr(tmp, "class")[2] == "symbolic") == 
                TRUE, lbs[length(lbs) + 1] <- lb[which(clu == 
                i)[1]], lbs[length(lbs) + 1] <- as.numeric(dimnames(x)[[1]])[which(clu == 
                i)[1]])
        }
        rm(i)
    }
    for (i in 1:dim(bm)[1]) {
        for (j in 1:dim(bm)[1]) {
            if (isTRUE(bm[i, j] %in% lbs) == FALSE) {
                for (l in 1:length(lbs)) {
                  if (isTRUE(attr(tmp, "class")[2] == "symbolic") == 
                    TRUE) {
                    if (isTRUE(bm[i, j] %in% lbs[l]) == TRUE) 
                      bm[i, j] <- lbs[l]
                  }
                  else if (isTRUE(attr(tmp, "class")[2] == "numerical") == 
                    TRUE) {
                    if (isTRUE(clu[bm[i, j]] %in% as.numeric(lbs)[l]) == 
                      TRUE) 
                      bm[i, j] <- l
                  }
                }
                rm(l)
            }
        }
        rm(j)
    }
    rm(i)
    dimnames(bm)[[1]] <- dimnames(bm)[[2]] <- as.list(lbs)
    return(as.data.frame(bm))
}
