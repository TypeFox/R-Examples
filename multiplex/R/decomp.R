decomp <-
function (S, x, type = c("mc", "pi", "cc"), reduc = FALSE) 
{
    if (isTRUE(attr(S, "class")[1] == "Semigroup") == FALSE) 
        stop("\"S\" should be an object of a \"Semigroup\" class.")
    if (isTRUE(attr(x, "class") == "Pi.rels" || attr(x, "class")[1] == 
        "Congruence") == FALSE) 
        stop("\"x\" should be an object either of a \"Pi.rels\" or a \"Congruence\" class.")
    if (isTRUE(attr(x, "class") == "Pi.rels") == TRUE) {
        switch(match.arg(type), mc = poi <- x$mc, pi = poi <- x$pi, 
            cc = stop("Type options for a \"Pi.rels\" class object should be either \"pi\" or \"mc\""))
    }
    else if (isTRUE(attr(x, "class")[1] == "Congruence") == TRUE) {
        if (match.arg(type) != "cc") {
            warning("Other type options than \"cc\" are not suitable for the \"Congruence\" class object")
        }
        ifelse(isTRUE(attr(x, "class")[2] == "PO.Semigroup") == 
            TRUE, poi <- x$PO, NA)
    }
    if (isTRUE(attr(x, "class") == "Pi.rels") == TRUE) {
        clu <- list()
        lb <- list()
        length(lb) <- length(clu) <- dim(poi)[3]
        for (i in 1:dim(poi)[3]) {
            ifelse(isTRUE(all(as.vector(poi[, , i]) == 1L) == 
                TRUE) == TRUE, clu[[i]] <- rep(1, S$ord), clu[[i]] <- stats::cutree(stats::hclust(stats::dist(poi[, 
                , i])), k = length(cut(stats::as.dendrogram(stats::hclust(stats::dist(poi[, 
                , i]))), h = 0L)$lower)))
            attr(clu[[i]], "names") <- S$st
            lb[[i]] <- list()
            for (j in 1:length(tabulate(clu[[i]]))) {
                lb[[i]][[j]] <- noquote(attr(which(clu[[i]] == 
                  j), "names"))
            }
            rm(j)
        }
        rm(i)
    }
    else {
        clu <- x$clu
        lb <- list()
        length(lb) <- length(clu)
        for (i in 1:length(clu)) {
            attr(clu[[i]], "names") <- S$st
            lb[[i]] <- list()
            for (j in 1:length(tabulate(clu[[i]]))) {
                lb[[i]][[j]] <- noquote(attr(which(clu[[i]] == 
                  j), "names"))
            }
            rm(j)
        }
        rm(i)
    }
    if (reduc) {
        im <- list()
        po <- list()
        dm <- vector()
        length(po) <- length(im) <- length(clu)
        for (i in 1:length(clu)) {
            im[[i]] <- reducs(S, clu = as.vector(clu[[i]]))
            if (isTRUE(attr(x, "class") == "Pi.rels") == TRUE) {
                po[[i]] <- reduc(poi[, , i], clu = as.vector(clu[[i]]), 
                  labels = dimnames(im[[i]])[[1]])
            }
            else if (isTRUE(attr(x, "class")[2] == "PO.Semigroup") == 
                TRUE) {
                po[[i]] <- reduc(poi, clu = as.vector(clu[[i]]), 
                  labels = dimnames(im[[i]])[[1]])
            }
            else {
                NA
            }
            dm[length(dm) + 1] <- dim(im[[i]])[1]
        }
        rm(i)
        for (k in 1:length(clu)) {
            imm <- as.matrix(im[[k]])
            for (i in 1:dim(imm)[1]) {
                for (j in 1:dim(imm)[1]) {
                  if (isTRUE(imm[i, j] %in% dimnames(imm)[[1]]) == 
                    FALSE) {
                    for (l in 1:length(lb[[k]])) {
                      if (isTRUE(attr(S, "class")[2] == "symbolic") == 
                        TRUE) {
                        if (isTRUE(imm[i, j] %in% lb[[k]][[l]]) == 
                          TRUE) 
                          imm[i, j] <- lb[[k]][[l]][1]
                      }
                      else if (isTRUE(attr(S, "class")[2] == 
                        "numerical") == TRUE) {
                        if (isTRUE(attr(clu[[k]], "names")[imm[i, 
                          j]] %in% lb[[k]][[l]]) == TRUE) 
                          imm[i, j] <- which(attr(clu[[k]], "names") == 
                            lb[[k]][[l]][1])
                      }
                    }
                  }
                }
                rm(j)
            }
            rm(i)
            im[[k]] <- as.data.frame(imm)
        }
        rm(k)
        for (i in 1:length(clu)) attr(clu[[i]], "names") <- NULL
        ifelse(isTRUE(attr(x, "class") == "Pi.rels" || attr(x, 
            "class")[2] == "PO.Semigroup") == TRUE, lst <- list(clu = clu, 
            eq = lb, IM = im, PO = po, dims = dm), lst <- list(clu = clu, 
            eq = lb, IM = im, dims = dm))
        ifelse(isTRUE(attr(x, "class")[1] == "Congruence") == 
            TRUE, class(lst) <- c("Decomp", "Reduc", attr(x, 
            "class")[1], "cc"), class(lst) <- c("Decomp", "Reduc", 
            attr(x, "class")[1], match.arg(type)))
        return(lst)
    }
    else {
        for (i in 1:length(clu)) attr(clu[[i]], "names") <- NULL
        lst <- list(clu = clu, eq = lb)
        ifelse(isTRUE(attr(x, "class")[1] == "Congruence") == 
            TRUE, class(lst) <- c("Decomp", attr(x, "class")[1], 
            "cc"), class(lst) <- c("Decomp", attr(x, "class")[1], 
            match.arg(type)))
        return(lst)
    }
}
