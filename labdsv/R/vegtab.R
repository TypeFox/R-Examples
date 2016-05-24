vegtab <- function(taxa,set,minval=1,pltord,spcord,pltlbl,trans=FALSE)
{
    if (missing(set)) {
        set <- seq(1:nrow(taxa))
    } else {
        set <- seq(1:nrow(taxa))[set]
        set <- set[!is.na(set)]
    }
    tmp <- taxa[set,]
    spcidx <- apply(tmp>0,2,sum)
    tmp <- tmp[,spcidx >= minval]

    if (missing(pltord)) {
        pltord <- seq(1:nrow(tmp))
    } else {
        pltord <- pltord[set]
    }
    if (missing(spcord)) {
        spcord <- -apply(tmp>0,2,sum)
    } else {
        spcord <- spcord[spcidx >= minval]
    }
    if (!missing(pltlbl)) {
        if (is.numeric(pltlbl)) {
            tmp <- cbind(pltlbl[set],tmp)
            dimnames(tmp)[[2]][1] <- deparse(substitute(pltlbl))
            spcord <- c(min(spcord)-1,spcord)
        } else {
            dimnames(tmp)[[1]] <- pltlbl
        }
    }
    tmp <- tmp[order(pltord),order(spcord)]
    if (trans==TRUE) {
        tmp <- t(tmp)
    }
    tmp
}
