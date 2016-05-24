"clt.examp" <-
function( n=1, reps=10000, nclass=16, norm.param=list(mean=0,sd=1),
          gamma.param=list(shape=1, rate=1/3), unif.param=list(min=0,max=1),
          beta.param=list(shape1=0.35, shape2=0.25)
         ) {
	# this function demonstrates the central limit theorem
	# by generating reps samples of size n from 4 different
	# distributions
	old.par <- par(oma=c(0,0,2,0), mfrow=c(2,2) )
	on.exit( par(old.par) )

	# Normal
        norm.param$n <- n*reps
	norm.mat <- matrix( do.call('rnorm',norm.param),
                             ncol=n )
	norm.mean <- rowMeans(norm.mat)

	x <- seq( min(norm.mean), max(norm.mean), length=50)
	normmax <- max( dnorm(x,mean(norm.mean),sd(norm.mean)) )
	tmp.hist <- hist( norm.mean, plot=FALSE , nclass=nclass)
	normmax <- max( tmp.hist$density, normmax )*1.05

	hist( norm.mean, main="Normal",xlab="x",col='skyblue'
             ,freq=FALSE,ylim=c(0,normmax), nclass=nclass)
	lines( x, dnorm(x,mean(norm.mean),sd(norm.mean)) )

	# gamma
        gamma.param$n <- n*reps
	exp.mat <- matrix( do.call('rgamma',gamma.param), ncol=n )
	exp.mean <- rowMeans(exp.mat)

	x <- seq( min(exp.mean), max(exp.mean), length=50)
	expmax <- max( dnorm(x,mean(exp.mean),sd(exp.mean)) )
	tmp.hist <- hist( exp.mean, plot=FALSE, nclass=nclass)
	expmax <- max( tmp.hist$density, expmax)*1.05

	hist( exp.mean, main="Gamma",xlab="x",col='skyblue',
             freq=FALSE,ylim=c(0,expmax), nclass=nclass)
	lines( x, dnorm(x,mean(exp.mean),sd(exp.mean)) )

	# Uniform
        unif.param$n <- n*reps
	unif.mat <- matrix( do.call('runif',unif.param),
                             ncol=n )
	unif.mean <- rowMeans(unif.mat)

	x <- seq( min(unif.mean), max(unif.mean), length=50)
	unimax <- max( dnorm(x,mean(unif.mean),sd(unif.mean)) )
	tmp.hist <- hist( unif.mean, plot=FALSE, nclass=nclass)
	unimax <- max( tmp.hist$density, unimax)*1.05

	hist( unif.mean, main="Uniform", xlab="x",col='skyblue',
             freq=FALSE,ylim=c(0,unimax), nclass=nclass)
	lines( x, dnorm(x,mean(unif.mean),sd(unif.mean)) )

	# Beta
        beta.param$n <- n*reps
	beta.mat <- matrix( do.call('rbeta',beta.param), ncol=n )
	beta.mean <- rowMeans(beta.mat)

	x <- seq( min(beta.mean), max(beta.mean), length=50)
	betamax <- max( dnorm(x,mean(beta.mean),sd(beta.mean)) )
	tmp.hist <- hist( beta.mean, plot=FALSE, nclass=nclass)
	betamax <- max( tmp.hist$density, betamax)

	hist( beta.mean, main="Beta", xlab="x",col='skyblue',
             freq=FALSE, ylim=c(0,betamax), nclass=nclass)
	lines( x, dnorm(x,mean(beta.mean),sd(beta.mean)) )

	mtext( paste("sample size =",n), outer=TRUE ,cex=2)
        invisible(NULL)
}

