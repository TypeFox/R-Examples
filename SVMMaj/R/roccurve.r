roccurve <- function(q,y=attr(q,'y'),class=1,...){
	if(is.null(y)) error('Input y not found')
	if(!is.factor(y)) y <- factor(y)
	
	neg     <- sum(q<0)
	classes <- levels(y)
	obs     <- summary(y)

	ys    <- y[sort.int(q*(-1)^(class!=1),index.return=TRUE)$ix]
	c1    <- (ys==classes[-class])/obs[-class]
	c2	  <- (ys==classes[ class])/obs[ class]

	plot( c(0,cumsum(c1)), c(0,cumsum(c2)) , xlab = 'False Positive Rate',
		ylab = 'True Positive Rate',
		main=paste('ROC-curve\n','Class =',classes[class],
						', AUC ='  ,formatC( sum( c1 * cumsum(c2) ) ,digit=4)),
		type='l',...)
	lines(c(0,1),c(0,1),lty=2,col='grey')
	lines( mean( cumsum(c1)[neg:(neg+1)]),
		 mean( cumsum(c2)[neg:(neg+1)]),
		type='p',col='grey')
}

auc <- function(q,y=attr(q,'y')){
	if(is.null(y)) error('Input y not found')
	if(!is.factor(y)) y <- factor(y)
	neg     <- sum(q<0)
	classes <- levels(y)
  obs     <- summary(y)

	ys    <- y[sort.int(q,index.return=TRUE)$ix]
	c1    <- (ys==classes[2])/obs[2]
	c2	  <- (ys==classes[1])/obs[1]
	auc	  <- sum( c1 * cumsum(c2) )
	return(auc)
}

