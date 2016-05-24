pr <- function(s,y,object=NULL,tpm,Rho,ispd=NULL) {

if(!is.null(object)) {
	tpm <- object$tpm
	Rho <- object$Rho
	ispd <- object$ispd
}
if(missing(y)) {
	y <- if(!is.null(object)) object$y else NULL
	if(is.null(y)) stop("No observation sequence supplied.\n")
}
if(is.null(ispd)) ispd <- revise.ispd(tpm)
y <- charList(y)
Rho <- check.yval(y,Rho)

if(!is.list(s)) s <- list(s)
nseq <- length(s)
if(!length(y)%in%c(1,nseq))
	stop(paste("Mismatch between number of state sequences\n",
                   "and number of observation sequences.\n"))
rslt <- numeric(nseq)
for(i in 1:nseq) {
	ii <- if(length(y)==1) 1 else i
	tmp <- log(ispd[s[[i]][1]]) + log(Rho[y[[ii]][1],s[[i]][1]])
	n <- length(y[[ii]])
	for(j in 2:n)
        	tmp <- tmp + log(tpm[s[[i]][j-1],s[[i]][j]]) +
                             log(Rho[y[[ii]][j],s[[i]][j]])
	rslt[i] <- tmp
}

fy  <- ffun(y, Rho)
lns <- sapply(y,length)
rp  <- recurse(fy, tpm, ispd, lns)
ll  <- log(rp$llc)

jstop <- 0
for(i in 1:nseq) {
        jstart <- jstop+1
        jstop  <- jstop + lns[i]
	rslt[i] <- exp(rslt[i]-sum(ll[jstart:jstop]))
}

rslt
}
