confus <- function (clustering, fitted)
{
    clustering <- clustify(clustering)

    numplt <- length(clustering)
    numclu <- length(levels(clustering))
    pred <- apply(fitted,1,which.max)
    res <- matrix(0,nrow=numclu,ncol=numclu)
    for (i in 1:numplt) {
        res[clustering[i],pred[i]] <- res[clustering[i],pred[i]] + 1
    }
    correct <- sum(diag(res))
    percent <- correct/numplt
    rowsum <- apply(res,1,sum)
    colsum <- apply(res,2,sum)
    summar <- sum(rowsum*colsum)
    kappa <- ((numplt*correct) - summar) / (numplt^2 - summar)
    out <- list(confus=res,correct=correct,percent=percent,kappa=kappa)
    out
}

fuzconfus <- function (part, fitted, dis)
{
    clustering <- clustify(part)
    numplt<- length(clustering)
    numclu <- length(levels(clustering))
    pred <- apply(fitted,1,which.max)
    tmp <- part$ctc/diag(part$ctc)
    fuzerr <-  1- matrix(pmin(1,tmp),ncol=ncol(tmp))
    diag(fuzerr) <- 1

    res <- matrix(0,nrow=numclu,ncol=numclu)
    for (i in 1:numplt) {
        res[clustering[i],pred[i]] <- res[clustering[i],pred[i]] + 1
    }

    fuzres <- res

    for (i in 1:ncol(res)) {
        for (j in 1:ncol(res)) {
            if (i != j) {
                fuzres[i,j] <- res[i,j] * fuzerr[i,j]
                fuzres[i,i] <- fuzres[i,i] + res[i,j] * (1-fuzerr[i,j])
            }
        }
    }

    correct <- sum(diag(fuzres))
    percent <- correct/numplt
    rowsum <- apply(fuzres,1,sum)
    colsum <- apply(fuzres,2,sum)
    summar <- sum(rowsum*colsum)
    kappa <- ((numplt*correct) - summar) / (numplt^2 - summar)
    fuzres <- data.frame(fuzres)
    names(fuzres) <- as.character(names(table((clustering))))
    out <- list(confus=fuzres,correct=correct,percent=percent,kappa=kappa)
    out
}

