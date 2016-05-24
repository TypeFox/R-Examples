sprt <-
function (x, k, l) 
{
    if (isTRUE(is.matrix(x)) == FALSE) 
        x <- as.matrix(x)
    pr <- vector()
    pr <- append(pr, paste(x[k, k], x[l, k], sep = ", "))
    pr <- append(pr, paste(x[k, l], x[l, l], sep = ", "))
    pr <- transl(pr)
    if (length(pr) == 2 && (is.na(strsplit(pr[1], ", ")[[1]][2]) & 
        is.na(strsplit(pr[2], ", ")[[1]][2]))) {
        pr <- jnt(pr)
    }
    clu <- rep(0, nrow(x))
    for (i in 1:length(pr)) {
        clu[as.numeric(strsplit(levels(factor(pr[i])), ", ")[[1]])] <- i
    }
    rm(i)
    clu[which(clu == 0)] <- max(clu) + 1
    X <- as.semigroup(x)
    if (isTRUE(any(is.na(reducs(X, clu)))) == FALSE) {
        return(clu)
    }
    else {
        if (length(pr) == 1) {
            pr <- paste(pr, pr, sep = ", ")
        }
        for (i in 1:length(pr)) {
            ifelse(is.na(strsplit(pr[i], ", ")[[1]][2]), pr[i] <- paste(pr[i], 
                pr[i], sep = ", "), NA)
        }
        rm(i)
        for (i in 1:length(pr)) {
            pr <- append(pr, paste(x[(as.numeric(strsplit(pr[i], 
                ", ")[[1]][1])), k], x[(as.numeric(strsplit(pr[i], 
                ", ")[[1]][2])), k], sep = ", "))
            for (j in 1:length(pr)) {
                pr <- append(pr, paste(x[(as.numeric(strsplit(pr[i], 
                  ", ")[[1]][1])), (as.numeric(strsplit(pr[j], 
                  ", ")[[1]][1]))], x[(as.numeric(strsplit(pr[i], 
                  ", ")[[1]][1])), (as.numeric(strsplit(pr[j], 
                  ", ")[[1]][2]))], sep = ", "))
                pr <- append(pr, paste(x[(as.numeric(strsplit(pr[i], 
                  ", ")[[1]][2])), (as.numeric(strsplit(pr[j], 
                  ", ")[[1]][1]))], x[(as.numeric(strsplit(pr[i], 
                  ", ")[[1]][2])), (as.numeric(strsplit(pr[j], 
                  ", ")[[1]][2]))], sep = ", "))
            }
            rm(j)
            pr <- append(pr, paste(x[(as.numeric(strsplit(pr[i], 
                ", ")[[1]][1])), (as.numeric(strsplit(pr[i], 
                ", ")[[1]][2]))], x[(as.numeric(strsplit(pr[i], 
                ", ")[[1]][2])), (as.numeric(strsplit(pr[i], 
                ", ")[[1]][1]))], sep = ", "))
        }
        rm(i)
        tpr <- transl(pr)
        clu <- rep(0, nrow(x))
        for (i in 1:length(tpr)) {
            clu[as.numeric(strsplit(levels(factor(tpr[i])), ", ")[[1]])] <- i
        }
        rm(i)
        flt <- which(clu == 0)
        ifelse(isTRUE(length(flt) == 0) == FALSE, NA, return(clu))
        clu[which(clu == 0)] <- max(clu) + 1
        if (isTRUE(any(is.na(reducs(X, clu)))) == FALSE && isTRUE(sum(clu) == 
            ncol(x)) == FALSE) {
            return(clu)
        }
        else {
            xf <- x[flt, flt]
            tmp <- vector()
            for (i in 1:length(diag(xf))) {
                tmp <- append(tmp, paste(attr(diag(xf)[i], "names"), 
                  diag(xf)[i], sep = ", "))
            }
            rm(i)
            ttpr <- transl(c(tmp, tpr))
            clu <- rep(0, nrow(x))
            for (i in 1:length(ttpr)) {
                clu[as.numeric(strsplit(levels(factor(ttpr[i])), 
                  ", ")[[1]])] <- i
            }
            rm(i)
            clu[which(clu == 0)] <- max(clu) + 1
            Sfr <- reducs(X, clu)
            if (isTRUE(any(is.na(Sfr))) == FALSE) {
                return(clu)
            }
            else {
                ftt <- which(is.na(Sfr), arr.ind = TRUE)
                vect <- vector()
                for (i in 1:nrow(ftt)) {
                  vect <- append(vect, (c(ttpr, jnt(levels(factor(x[which(clu == 
                    ftt[i, 1]), which(clu == ftt[i, 2])]))))))
                }
                rm(i)
                vect <- transl(vect)
                clu <- rep(0, nrow(x))
                for (i in 1:length(vect)) {
                  clu[as.numeric(strsplit(levels(factor(vect[i])), 
                    ", ")[[1]])] <- i
                }
                rm(i)
                clu[which(clu == 0)] <- max(clu) + 1
                while (isTRUE(any(is.na(eval(Sfr2 <- reducs(X, 
                  clu))))) == TRUE) {
                  ftt2 <- which(is.na(Sfr2), arr.ind = TRUE)
                  if (isTRUE(nrow(ftt2) > 0) == TRUE) {
                    vect2 <- vector()
                    for (i in 1:nrow(ftt2)) {
                      vect2 <- append(vect2, (c(vect, jnt(levels(factor(x[which(clu == 
                        ftt2[i, 1]), which(clu == ftt2[i, 2])]))))))
                    }
                    rm(i)
                    vect2 <- transl(vect2)
                    clu <- rep(0, nrow(x))
                    for (i in 1:length(vect2)) {
                      clu[as.numeric(strsplit(levels(factor(vect2[i])), 
                        ", ")[[1]])] <- i
                    }
                    rm(i)
                    clu[which(clu == 0)] <- max(clu) + 1
                  }
                }
                return(clu)
            }
        }
    }
}
