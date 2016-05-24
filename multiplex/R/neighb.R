neighb <-
function (x, rs, type = c("und", "inn", "out"), inclx = FALSE, 
    k = 1, expand = FALSE) 
{
    if (isTRUE(attr(rs, "class") == "Rel.System") == FALSE) {
        if (is.array(rs) == FALSE) {
            stop("Data must be at least a stacked array of square matrices.")
        }
        warning("'rs' is transformed into a entire relational system of a \"Rel.System\" class.")
        rs <- rel.sys(rs, bonds = "entire")
    }
    if (isTRUE(k < 0L) == TRUE) 
        stop("'k' must not be negative.")
    if (isTRUE(k == 0L) == TRUE) 
        ifelse((inclx), return(x), return(character(0)))
    if (isTRUE(all(x %in% unique(unlist(dhc(as.character(rs$nodes)))))) == 
        TRUE) {
        if (isTRUE(length(rs$Ties) > 0L) == TRUE) {
            rst <- as.list(unlist(rs$Ties))
            srs <- list()
            for (i in 1:length(rst)) {
                tmp <- vector()
                if (length(rst[[i]]) > 0L) {
                  for (n in 1:length(x)) {
                    for (j in 1:length(rst[[i]])) {
                      if (x[n] %in% c(c(strsplit(rst[[i]][j], 
                        rs$prsep)[[1]][1], strsplit(rst[[i]][j], 
                        rs$prsep)[[1]][2]))) {
                        tmp <- append(tmp, rst[[i]][j])
                      }
                    }
                    rm(j)
                  }
                  rm(n)
                }
                srs[[i]] <- as.vector(unlist(tmp))
            }
            rm(i)
            attr(srs, "names") <- attr(rst, "names")
            nrs <- vector()
            for (i in 1:length(srs)) {
                if (isTRUE(length(srs[[i]]) > 0L) == TRUE) {
                  for (j in 1:length(srs[[i]])) {
                    switch(match.arg(type), und = nrs <- append(nrs, 
                      strsplit(srs[[i]][j], rs$prsep)[[1]][1]), 
                      inn = nrs <- append(nrs, strsplit(srs[[i]][j], 
                        rs$prsep)[[1]][1]), out = nrs <- append(nrs, 
                        (strsplit(srs[[i]][j], rs$prsep)[[1]][2])))
                    switch(match.arg(type), und = nrs <- append(nrs, 
                      strsplit(srs[[i]][j], rs$prsep)[[1]][2]), 
                      inn = NA, out = NA)
                  }
                  rm(j)
                }
            }
            rm(i)
            levels(factor(nrs))
            nb <- levels(factor(nrs))
            if (isTRUE(k > 1L) == TRUE) {
                if (isTRUE(expand == FALSE) == TRUE) {
                  for (K in 2:k) {
                    nb <- append(nb, ngbs(nb, rs, type = type))
                  }
                  rm(K)
                }
                else if (isTRUE(expand == TRUE) == TRUE) {
                  nb2 <- nb
                  nbk <- list()
                  if (!(inclx)) {
                    nbk[[1]] <- nb2[which(!(nb2 %in% x))]
                    ink <- 2L
                  }
                  else {
                    nbk[[1]] <- x
                    nbk[[2]] <- nb2[which(!(nb2 %in% x))]
                    ink <- 3
                    k <- k + 1L
                  }
                  for (K in ink:k) {
                    nb2 <- append(nb2, ngbs(nb2, rs, type = type))
                    nbk[[K]] <- nb2[which(!(nb2 %in% c(nb, unlist(nbk))))]
                  }
                  rm(K)
                  ifelse(!(inclx), attr(nbk, "names") <- paste0("k=", 
                    1:k), attr(nbk, "names") <- paste0("k=", 
                    seq_along((ink - 2L):k) - 1L))
                }
            }
            else {
                NA
            }
        }
        else if (isTRUE(length(rs$Ties) > 0L) == FALSE) {
            nb <- x
        }
        if (isTRUE(expand == FALSE) == TRUE) {
            ifelse(!(inclx), return(nb[which(!(nb %in% x))]), 
                return(nb))
        }
        else if (isTRUE(expand == TRUE) == TRUE) {
            return(nbk)
        }
    }
    else {
        warning("'x' is not part of the relational system provided.")
        x
    }
}
