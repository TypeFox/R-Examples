transl <-
function (lt, prsep = ", ") 
{
    llt <- levels(factor(lt))
    if (isTRUE(length(llt) == 1) == TRUE) {
        ifelse(isTRUE(length(jnt(llt)) == 1) == TRUE, return(jnt(llt)), 
            return(llt))
    }
    else {
        Ls <- as.list(llt)
        j <- 1
        tmp0 <- strsplit(Ls[[1]], prsep)[[1]]
        attr(Ls[[1]], "names") <- j
        for (i in 2:length(Ls)) {
            tmp2 <- strsplit(Ls[[i]], prsep)[[1]]
            ifelse((any(tmp2 %in% tmp0)) | (any(tmp0 %in% tmp2)), 
                attr(Ls[[i]], "names") <- as.integer(attr(Ls[[j]], 
                  "names")), NA)
        }
        rm(i)
        uno <- unlist(dhc(jnt(unlist(Ls)[which(attr(unlist(Ls), 
            "names") == j)])))
        for (i in which(attr(unlist(Ls), "names") == "")) {
            tmp2 <- strsplit(Ls[[i]], prsep)[[1]]
            ifelse((any(uno %in% tmp2)) | (any(tmp2 %in% uno)), 
                attr(Ls[[i]], "names") <- j, NA)
            uno <- unlist(dhc(jnt(unlist(Ls)[which(attr(unlist(Ls), 
                "names") == j)])))
        }
        rm(i)
        if (isTRUE(length(uno) > 1) == TRUE) {
            while (isTRUE(unlist(dhc(jnt(unlist(Ls)[which(attr(unlist(Ls), 
                "names") == j)]))) == uno) == TRUE) {
                for (i in which(attr(unlist(Ls), "names") == 
                  "")) {
                  tmp2 <- strsplit(Ls[[i]], prsep)[[1]]
                  ifelse((any(uno %in% tmp2)) | (any(tmp2 %in% 
                    uno)), attr(Ls[[i]], "names") <- j, NA)
                }
                rm(i)
            }
        }
        while (length(which(attr(unlist(Ls), "names") == "")) != 
            0) {
            j <- (j + 1)
            attr(Ls[[which(attr(unlist(Ls), "names") == "")[1]]], 
                "names") <- j
            if (length(which(attr(unlist(Ls), "names") == "")) != 
                0) {
                tmp0 <- unlist(dhc(jnt(unlist(Ls)[which(attr(unlist(Ls), 
                  "names") == j)])))
                for (i in which(attr(unlist(Ls), "names") == 
                  "")) {
                  tmp2 <- strsplit(Ls[[i]], prsep)[[1]]
                  ifelse((any(tmp2 %in% tmp0)) | (any(tmp0 %in% 
                    tmp2)), attr(Ls[[i]], "names") <- j, NA)
                }
                rm(i)
                tmp0 <- unlist(dhc(jnt(unlist(Ls)[which(attr(unlist(Ls), 
                  "names") == j)])))
                for (i in which(attr(unlist(Ls), "names") == 
                  "")) {
                  tmp2 <- strsplit(Ls[[i]], prsep)[[1]]
                  ifelse((any(tmp0 %in% tmp2)) | (any(tmp2 %in% 
                    tmp0)), attr(Ls[[i]], "names") <- j, NA)
                }
                rm(i)
            }
        }
        clu <- vector()
        for (i in 1:length(Ls)) clu[i] <- as.numeric(attr(Ls[[i]], 
            "names")[1])
        tls <- vector()
        for (i in 1:nlevels(factor(clu))) tls <- append(tls, 
            jnt(llt[clu == i]))
        return(levels(factor(tls)))
    }
}
