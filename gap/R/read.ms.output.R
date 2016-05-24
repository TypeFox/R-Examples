read.ms.output <- function(msout, is.file=TRUE, xpose=TRUE, verbose=TRUE, outfile=NULL, outfileonly=FALSE)
{
    if (is.file) msout <- scan(file=msout, what=character(0), sep="\n", quiet=TRUE)
    if (is.na(msout[1])) stop("Usage: read.ms.output(msout), or read.ms.output(filename)")
    nsam <- as.integer(strsplit(msout[1], split=" ")[[1]][2])
    ndraws <- as.integer(strsplit(msout[1], split=" ")[[1]][3])
    result <- gamlist <- positions <- list()
    marker <- grep("prob",msout)
    probs <- sapply(strsplit(msout[marker], split=":"), function(vec) as.numeric(vec[2]))
    marker <- grep("time",msout)
    times <- sapply(strsplit(msout[marker], split="\t"), function(vec){ as.numeric(vec[2:3])})
    marker <- grep("segsites", msout)
    stopifnot(length(marker) == ndraws)
    segsites <- sapply(strsplit(msout[marker], split=" "), function(vec) as.integer(vec[2]))
    if (!is.null(outfile)) of <- file(outfile,"w")
    for (draw in seq(along=marker))
    {
        if (verbose) if (!(draw %% 1000)) cat(draw, " ")
        if (segsites[draw] > 0)
        {
            tpos <- strsplit(msout[marker[draw]+1], split=" ")
            positions[[draw]] <- as.numeric(tpos[[1]][2:(segsites[draw]+1)]) 
            haplotypes <- msout[(marker[draw] + 2):(marker[draw] + 2 + nsam - 1)]
            if (!is.null(outfile)) cat(paste(draw,1:nsam,haplotypes,sep="\t"),file=of,sep="\n")
            haplotypes <- strsplit(haplotypes, split="")
            h <- sapply(haplotypes, function(el) c(as.integer(el)))
            if(segsites[draw] == 1) h <- t(as.matrix(h))
        }
        else {
            h <- matrix(nrow=0, ncol=nsam)
            positions[[draw]]<- NA
        }
        if (xpose)
        {
           s <- h
           colnames(s) <- 1:nsam
        }
        else {
           s <- t(h)
           row.names(s) <- 1:nsam
        }
        if (outfileonly) gamlist[[draw]] <- NA else gamlist[[draw]] <- s
        stopifnot(all(dim(h) == c(segsites[draw],nsam))) 
    }
    if (verbose) cat("\n")
    if (!is.null(outfile)) close(of)
    z <- list(call=msout[1], seed=as.numeric(strsplit(msout[2]," ")[[1]]), nsam=nsam, nreps=ndraws,
              segsites=segsites, probs=probs, times=t(times), positions=positions, gametes=gamlist)
    invisible(z)
}
