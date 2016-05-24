`calcCondProb` <-
function(object,conditionalp=0.5) {
    if (!inherits(object, "randomLCA"))
        stop("Use only with 'randomLCA' objects.\n")
     if (missing(conditionalp)) conditionalp <- 0.5

    if (object$level2) blocksize <- object$level2size
    else blocksize <- object$blocksize
    
	nblocks <- dim(object$outcomep)[2]/blocksize

	if (object$probit) outcomex <- rep(as.vector(qnorm(object$outcomep)),each=length(conditionalp))
	else outcomex <- rep(log(as.vector(object$outcomep)/(1-as.vector(object$outcomep))),each=length(conditionalp))
	# get the offset corresponding to each percentile
	offset <- rep(qnorm(conditionalp),times=prod(dim(object$outcomep)))
	# get the variances corresponding to each 
	if (!is.null(object$lambdacoef)) {
		lambdacoef <- rep(as.vector(object$lambdacoef),each=ifelse(object$byclass,1,object$nclass)*length(conditionalp),
		times=dim(object$outcomep)[2]/blocksize)
		outcomex <- outcomex+offset*lambdacoef
	}
	if (object$probit) outcomep <- pnorm(outcomex)
	else outcomep <- 1/(1+exp(-outcomex))
	perc <- factor(rep(conditionalp,times=length(as.vector(object$outcomep))))
 
 	outcome <- factor(rep(1:blocksize,each=length(conditionalp)*object$nclass,times=nblocks))
	 class <- factor(rep(1:object$nclass,times=blocksize*nblocks,each=length(conditionalp)))
	 block <- factor(rep(1:nblocks,each=object$nclass*blocksize*length(conditionalp)))
	conddata <- data.frame(perc,class,block,outcome,outcomep)
	 conddata
 }

