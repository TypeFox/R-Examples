ngbs <-
function (x, rs, type = c("und", "inn", "out"), num = FALSE) 
{
    if (isTRUE(attr(rs, "class") == "Rel.System") == FALSE) 
        stop("'rs' must be a relational system of a \"Rel.System\" class.")
    ifelse(isTRUE(is.numeric(x)) == TRUE, x <- rs$nodes[x], NA)
    if (isTRUE(all(x %in% unique(unlist(dhc(as.character(rs$nodes)))))) == 
        TRUE) {
        rst <- as.list(unlist(rs$Ties))
        srs <- list()
        for (i in 1:length(rst)) {
            tmp <- vector()
            if (length(rst[[i]]) > 0) {
                for (n in 1:length(x)) {
                  for (j in 1:length(rst[[i]])) {
                    if (x[n] %in% c(c(strsplit(rst[[i]][j], rs$prsep)[[1]][1], 
                      strsplit(rst[[i]][j], rs$prsep)[[1]][2]))) {
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
            if (isTRUE(length(srs[[i]]) > 0) == TRUE) {
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
        nb <- levels(factor(nrs))
        if (num) {
            return(which(rs$nodes %in% nb[which(!(nb %in% x))]))
        }
        else {
            return(nb[which(!(nb %in% x))])
        }
    }
    else {
        x
    }
}
