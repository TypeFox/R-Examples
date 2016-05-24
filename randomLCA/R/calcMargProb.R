`calcMargProb` <-
function(object) {

    if (!inherits(object, "randomLCA"))
        stop("Use only with 'randomLCA' objects.\n")
     if (!object$random)
        stop("Object must be random effects model.\n")
	if (object$probit) outcomex <- qnorm(object$outcomep)
	else outcomex <- log(object$outcomep/(1-object$outcomep))
	if (object$level2) blocksize <- object$level2size
	else blocksize <- object$blocksize
	nblocks <- dim(object$outcomep)[2]/blocksize
  #browser()
    if (object$level2) {
    	outcomep <- NULL
    	for (i in 1:object$nclass) { 
    		if (object$byclass) lambdacoef <- object$lambdacoef[i,]
 			else lambdacoef <-object$lambdacoef
			if (object$byclass) taucoef <- object$taucoef[i]
			else taucoef <- object$taucoef
					
			# 1 and 2 are the points
			# 3 and 4 are the weights
			
			points <- cbind(expand.grid(norm.gauss.hermite(51)[,1],norm.gauss.hermite(51)[,1]),
				expand.grid(norm.gauss.hermite(51)[,2],norm.gauss.hermite(51)[,2]))
				
			if (object$probit) {
				wprobs <- t(apply(points,1, function(x) {
					x[3]*x[4]*pnorm(outcomex[i,]+(x[1]+taucoef*x[2])*rep(lambdacoef,nblocks))
					}))
			} else {
				wprobs <- t(apply(points,1, function(x) {
					x[3]*x[4]/(1+exp(-outcomex[i,]-(x[1]+taucoef*x[2])*rep(lambdacoef,nblocks)))
					}))
			}			
			outcomep <- c(outcomep,apply(wprobs,2,sum))	
		}
		outcome <- factor(rep(1:blocksize,times=nblocks))
		class <- factor(rep(1:object$nclass,each=blocksize*nblocks))
		block <- factor(rep(1:nblocks,each=blocksize,times=object$nclass))
    } else
    {
    	outcomep <- NULL
    	for (i in 1:object$nclass) { 
    		if (object$byclass) lambdacoef <- rep(as.vector(object$lambdacoef[i,]),
				times=dim(object$outcomep)[2]/blocksize)
			else lambdacoef <- rep(as.vector(object$lambdacoef),
				times=dim(object$outcomep)[2]/blocksize)
			if (object$probit) probs <- apply(as.matrix(norm.gauss.hermite(51)[,1]),1,function(x)
				pnorm(outcomex[i,]+x*lambdacoef))
			else probs <- apply(as.matrix(norm.gauss.hermite(51)[,1]),1,function(x)
				1/(1+exp(-outcomex[i,]-x*lambdacoef)))
			outcomep <- c(outcomep,apply(t(t(probs)*norm.gauss.hermite(51)[,2]),1,sum))
		}
		outcome <- factor(rep(1:blocksize,times=nblocks))
		class <- factor(rep(1:object$nclass,each=blocksize*nblocks))
		block <- factor(rep(1:nblocks,each=blocksize,times=object$nclass))
	}
	outcomep <- ifelse(outcomep>1-1.0e-14,1-1.0e-14,outcomep)
	outcomep <- ifelse(outcomep<1.0e-14,1.0e-14,outcomep)
	margdata <- data.frame(class,block,outcome,outcomep)
	margdata
}

