jnt <-
function (xj, prsep = ", ") 
{
    if (isTRUE(length(xj) != 0) == TRUE) {
        Ltj <- FALSE
        ifelse(isTRUE(is.list(xj)) == TRUE, Ltj <- TRUE, NA)
        if (isTRUE(Ltj) == TRUE) {
            Xj <- list()
            jt <- list()
            length(Xj) <- length(jt) <- length(xj)
            for (i in 1:length(xj)) {
                if (isTRUE(length(xj[[i]]) != 0) == TRUE) {
                  tmpj <- as.list(xj[[i]])
                  for (j in 1:length(xj[[i]])) {
                    Xj[[i]] <- append(Xj[[i]], strsplit(tmpj[[j]], 
                      prsep)[[1]])
                  }
                  rm(j, tmpj)
                  Xj[[i]] <- unique(Xj[[i]])
                  if (length(Xj[[i]]) == 1) 
                    jt[[i]] <- Xj[[i]]
                  if (length(Xj[[i]]) > 1) 
                    jt[[i]] <- paste(Xj[[i]][1], Xj[[i]][2], 
                      sep = prsep)
                  if (length(Xj[[i]]) > 2) {
                    for (j in 3:length(Xj[[i]])) jt[[i]] <- paste(jt[[i]], 
                      Xj[[i]][j], sep = prsep)
                    rm(j)
                  }
                }
                else {
                  NA
                }
            }
            rm(i)
        }
        else if (isTRUE(Ltj) == FALSE) {
            vec <- vector()
            for (i in 1:length(xj)) {
                vec <- append(vec, strsplit(xj[i], prsep)[[1]])
            }
            rm(i)
            vec <- levels(factor(vec))
            if (length(vec) == 1) 
                jt <- vec
            if (length(vec) > 1) 
                jt <- paste(vec[1], vec[2], sep = prsep)
            if (length(vec) > 2) {
                for (i in 3:length(vec)) jt <- paste(jt, vec[i], 
                  sep = prsep)
                rm(i)
            }
        }
        ifelse(isTRUE(is.list(xj) == TRUE) == TRUE, attr(jt, 
            "names") <- attr(xj, "names"), NA)
        return(jt)
    }
    else {
        xj
    }
}
