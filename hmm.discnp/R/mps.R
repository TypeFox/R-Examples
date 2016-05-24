mps <- function(y,object=NULL,tpm,Rho,ispd=NULL) {
#
# Function mps: most probable states.
#

if(!is.null(object)) {
	tpm  <- object$tpm
	Rho  <- object$Rho
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
nseq <- length(y)
fy   <- ffun(y,Rho)
rp   <- recurse(fy, tpm, ispd, lns)
xxx  <- apply(rp$gamma, 2, which.max)
if(nseq==1) return(xxx)
rslt  <- list()
jstop <- 0
for(i in 1:nseq) {
        jstart <- jstop+1
        jstop  <- jstop + lns[i]
	rslt[[i]] <- xxx[jstart:jstop]
}
rslt
}
