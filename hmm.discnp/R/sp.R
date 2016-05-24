sp <- function (y, object = NULL, tpm, Rho, ispd=NULL, means=FALSE)
{
    if (!is.null(object)) {
        tpm <- object$tpm
        Rho <- object$Rho
        ispd <- object$ispd
    }
    if(is.null(ispd)) ispd <- revise.ispd(tpm)
    if(missing(y)) {
	y <- if(!is.null(object)) object$y else NULL
	if(is.null(y)) stop("No observation sequence supplied.\n")
    }
    y    <- charList(y)
    Rho  <- check.yval(y,Rho)
    lns  <- sapply(y,length)
    fy   <- ffun(y, Rho)
    rp   <- recurse(fy, tpm, ispd, lns)
    prbs <- rp$gamma
    if(means) {
	yval <- as.numeric(row.names(Rho))
	if(any(is.na(yval)))
		stop("Non-numeric y-values; means make no sense.\n")
        cmns <- apply(yval*Rho,2,sum)
	mns  <- apply(cmns*prbs,2,sum)
    }
    nseq <- length(lns)
    if (nseq == 1) {
	if(means) return(list(probs=prbs,means=mns))
	return(prbs)
    }
    xxx <- list()
    if(means) yyy <- list()
    istop <- 0
    for(i in 1:nseq) {
        istart <- istop+1
        istop  <- istop + lns[i]
        xxx[[i]] <- prbs[,istart:istop]
        if(means) yyy[[i]] <- mns[istart:istop]
    }
    if(means) return(list(probs=xxx,means=yyy))
    xxx
}
