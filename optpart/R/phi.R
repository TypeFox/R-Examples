phi <- function (taxa,clustering,minplt=10) 
{
    clustering <- clustify(clustering)
    names.out <- levels(clustering)
    clustering <- as.numeric(clustering)

    N <- nrow(taxa)
    N.p <- table(clustering)
    numplt <- apply(taxa>0,2,sum)
    taxa <- taxa[,numplt>=minplt]
    out <- matrix(NA,nrow=ncol(taxa),ncol=length(N.p))

    for (i in 1:ncol(taxa)) {
        n <- sum(taxa[,i]>0)
        n.p <- tapply(taxa[,i]>0,clustering,sum)
        numer <- (N*n.p) - (n*N.p)
        denom <- 
           sqrt(as.numeric(n*N.p)*as.numeric((N-n)*(N-N.p)))
        out[i,] <- numer/denom
     }
     out <- data.frame(out)
     row.names(out) <- names(taxa)
     names(out) <- names.out
     out
}
