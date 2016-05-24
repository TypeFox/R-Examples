`calcCond2Prob` <-
function(object,conditionalp=0.5) {

    if (!inherits(object, "randomLCA"))
        stop("Use only with 'randomLCA' objects.\n")
     if (!object$random)
        stop("Object must be random effects model.\n")
     if (!object$level2)
        stop("Object must be 2 level model.\n")
     if (missing(conditionalp)) conditionalp <- 0.5

	if (object$probit) outcomex <- qnorm(object$outcomep)
	else outcomex <- log(object$outcomep/(1-object$outcomep))
    nblocks <- dim(object$outcomep)[2]/object$level2size
# need to integrate over tau
	outcomep <- NULL
	for (ip in 1:length(conditionalp)) {
		offset <- qnorm(conditionalp[ip])
		for (i in 1:object$nclass) { 
			if (object$byclass) lambdacoef <- rep(as.vector(object$lambdacoef[i,]),
				times=dim(object$outcomep)[2]/object$level2size)
			else lambdacoef <- rep(as.vector(object$lambdacoef),
				times=dim(object$outcomep)[2]/object$level2size)
			if (object$byclass) taucoef <- object$taucoef[i]
			else taucoef <- object$taucoef
	
			if (object$probit) probs <- apply(as.matrix(norm.gauss.hermite(51)[,1]),1,function(x)
			  pnorm(outcomex[i,]+(offset+x*taucoef)*lambdacoef))
			else probs <- apply(as.matrix(norm.gauss.hermite(51)[,1]),1,function(x)
			  1/(1+exp(-outcomex[i,]-(offset+x*taucoef)*lambdacoef)))
			outcomep <- c(outcomep,apply(t(t(probs)*norm.gauss.hermite(51)[,2]),1,sum))
		}
	}
	outcome <- factor(rep(rep(1:object$level2size,times=nblocks),length(conditionalp)))
	class <- factor(rep(rep(1:object$nclass,each=object$level2size*nblocks),length(conditionalp)))
	block <- factor(rep(rep(1:nblocks,each=object$level2size,times=object$nclass),length(conditionalp)))
	perc <- factor(rep(conditionalp,each=length(as.vector(object$outcomep))))
	conddata <- data.frame(perc,class,block,outcome,outcomep)
	conddata
}
