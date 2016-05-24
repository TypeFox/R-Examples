power.examp <-
function(n=1, stdev=1, diff=1, alpha=0.05, xmin=-2, xmax=4)
{
	old.par <- par(mfrow=c(2,1), oma=c(0,0,3.1,0) )
	on.exit(par(old.par))


	n<-as.integer(n)
	stdev<-as.numeric(stdev)
	diff<-as.numeric(diff)
	alpha<-as.numeric(alpha)
	xmin<-as.numeric(xmin)
	xmax<-as.numeric(xmax)

	se <- stdev/sqrt(n)
	x <- seq( xmin, xmax, length=100 )

	# null hypothesis plots
	plot( x, dnorm(x,0,se), type="n", ylim=c(0, dnorm(0,0,se)*7/6), ylab="",
		main="Null Distribution")

	r <- qnorm(1-alpha,0,se)

	polygon( c(r, r, x[ x>r ]), c(0, dnorm(c(r,x[x>r]),0,se)), col='pink')

	abline(h=0)
	lines(x, dnorm(x,0,se), col='red' )

	abline(v=r)
	text(r,dnorm(0,0,se)*15/14, "--> rejection region", adj=0)
	axis(1,at=r, line=-0.75, cex=0.7)

	legend( par('usr')[2],par('usr')[4],xjust=1,bty='n',
	     fill='pink',legend=expression(alpha))

	# Alternative hypothesis plots
	plot( x, dnorm(x,0,se), type="n", ylim=c(0, dnorm(0,0,se)*7/6), ylab="",
		main="Alternative Distribution")

	polygon( c(r, r, x[ x>r ], max(x)), c(0, dnorm(c(r,x[x>r]),diff,se),0), col='lightblue')

	abline(h=0)

	lines(x, dnorm(x,diff,se), col='blue' )

	abline(v=r)
	text(r,dnorm(0,0,se)*15/14, "--> rejection region", adj=0)
	axis(1,at=r, line= -0.75, cex=0.7)

	legend( par('usr')[2],par('usr')[4],xjust=1,bty='n',
	     fill='lightblue',legend="Power")


	mtext(paste("se =",format(signif(se,3),nsmall=2),
                    "     z* =",format(signif(r,3),nsmall=2),
	"     power =", format(round( 1-pnorm(r,diff,se), 3 ),nsmall=2),
	"\n n =",format(n,width=3),"   sd =",format(stdev,nsmall=2),
	"   diff =",format(diff,nsmall=2),
                    "   alpha =",format(alpha,nsmall=3)
                    ), outer=TRUE, line=0, cex=1.5 )

	invisible( 1-pnorm(r,diff,se) )
}

