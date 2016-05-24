
#source( "~/Documents/elec_R/PPEB.R" )




## make contour plot of the optimization problem
## Note: alphas are multiplied by 100 to get in percents.
## Param     n = sample size
##           k = # errors of size d
##           d = mag of error (in percent)
##           e.max = mag of largest possible error (in percent)
##           xlim,ylim - area of plot to plot
trinomial.bound = function( n = 11, k = 2, d=40, e.max=100, 
				xlim=c(0.4,1), ylim=c(0,0.55),
				alpha.lvls=c(10),
				zero.threshold = 0.3,
				tick.lines=NULL,
				alpha.lwd = 2, 
				bold.first=FALSE,
				plot=TRUE,
				p.value.bound = NULL,
				grid.resolution=300,
				... ) {
					
					
					
## chance of getting something as or even less extreme
## than k errors of size d and no errors larger.
## (Step-Down Set)
prob.S = function( p_0, p_d = 1-p_0, n=11, k = 2 ) {
	
	if ( p_0 + p_d > 1 ) {
		NA
	} else {
		ks = 0:k
	
		sum( sapply( ks, function( x ) { 
			choose( n, x ) * ( p_0 ^ (n - x) ) * ( (p_d) ^ x ) 
		} ) )
	}	
}

## Derivative of above w.r.t p_0
d.prob.S.d.p_0 = function( p_0, p_d = 1-p_0, n=11, k = 2 ) {
	
	if ( p_0 + p_d > 1 ) {
		NA
	} else {
		ks = 0:k
	
		sum( sapply( ks, function( x ) { 
			choose( n, x ) * (n - x) * ( p_0 ^ (n - x - 1) ) * ( (p_d) ^ x ) 
		} ) )
	}
	
}

## Derivative of above w.r.t p_d
d.prob.S.d.p_d = function( p_0, p_d = 1-p_0, n=11, k = 2 ) {
	
	if ( p_0 + p_d > 1 ) {
		NA
	} else {
		ks = 1:k
	
		sum( sapply( ks, function( x ) { 
			choose( n, x ) * (x) * ( p_0 ^ (n - x) ) * ( (p_d) ^ (x-1) ) 
		} ) )
	}
}


## calc prob.S on vector of p_0s
prob.Ss = function( x, ... ) { sapply( x, prob.S, ... ) }

	
## the bound given the probability vector
## with, if desired, the lambda penalty for 
## the alpha constraint.
err.bound = function( p_0, p_d, n = 11, k = 2, d=40, 
					e.max=100, 
					penalize=FALSE, alpha=0.05, ... ) {
	if ( p_0 + p_d > 1 ) {
		NA
	} else {
		b = p_d * d + (1 - p_0 - p_d) * e.max
		if ( penalize ) {
			b = b - 100000 * (prob.S( p_0, p_d, n, k ) - alpha)^2
		}
		b
	}
}

## gradient of error bound function, with d/dp_0 and d/dp_d
grad.err.bound = function( x, d=40, e.max=100, 
						penalize=TRUE, alpha=0.05, ... ) {
	p_0 = x[1]
	p_d = x[2]
	
	if ( p_0 + p_d > 1 ) {
		c( NA, NA )
	} else {
		prS = prob.S( p_0, p_d )
		
		c( -e.max - 200000 * (prS - alpha) * d.prob.S.d.p_0(p_0, p_d, ... ),
			d - e.max - 200000 * (prS - alpha) * d.prob.S.d.p_d(p_0, p_d, ... ) )
	}
}



## the bound function that takes probs as a list
err.bound.pen = function( x, ... ) {
   err.bound( x[1], x[2], ... )
}

					
					
					
	p0s = seq( xlim[[1]], xlim[[2]], length.out = grid.resolution )
	pds = seq( ylim[[1]], ylim[[2]], length.out = grid.resolution )

	## PLOT THE ALPHA LINES
	vSS <- Vectorize(prob.S, c("p_0", "p_d"))

	# outer product, make 2D array applying vSS to all values of p0s and pds
	alphas = 100 * outer( p0s, pds, vSS, n=n, k=k )
	
	if ( plot ) {
		lwd = rep( alpha.lwd, length( alpha.lvls ) )
		if ( bold.first ) {
			lwd[1] = 2 * lwd[1]
		}
		contour( p0s, pds, alphas, levels=alpha.lvls, lwd=lwd, frame=FALSE, 
				xlab=expression(pi['0']), drawlabels=FALSE,
				ylab=bquote( pi['d'] ~ ", d =" ~ .(round(d,digits=2)) ),
				method="edge", ... )
	
		## BOUNDING BOX
		if ( ylim[2] < 1 ) {
			y0 = ylim[2]
		} else {
			y0 = 1
		}
		x0 = 1 - y0
		
		if ( ylim[1] > 0 ) {
			y1 = ylim[1]
		} else {
			y1 = 0	
		}
		x1 = 1 - y1

		if ( x0 < xlim[1] ) {
			x0 = xlim[1]
			y0 = 1 - x0
		}
		if ( x1 > xlim[2] ) {
			x1 = xlim[2]
			y1 = 1 - x1
		}
		
		segments( x0, y0, x1, y1 )
		segments( x0, 0, x1, 0 )
		
	}
	
	## THE BOUND CONTOUR LINES
	vBND = Vectorize( err.bound, c("p_0","p_d") )

	bounds = outer( p0s, pds, vBND, n=n, k=k, d=d)
	
	if ( plot ) {
		contour( p0s, pds, bounds, lty=3, add=TRUE )
		if ( !is.null( tick.lines ) ) {
			contour( p0s, pds, bounds, lty=5, levels=tick.lines, add=TRUE )
		}
	}
	
	## Finding the max along the main alpha line

	alpha = alpha.lvls[[1]]
	del = abs( alphas - alpha ) < zero.threshold
	zt = zero.threshold

	while( sum( del, na.rm=TRUE ) < 75 ) {
		zt = zero.threshold + zt
		warning( "zero threshold too small--increasing it to ", zt )
		del = abs( alphas - alpha ) < zt
	}

	del = !is.na(del) & del
	b = bounds
	b[ !del ] = 0
	loc = which.max(b) - 1
	loc.y = pds[ 1 + floor( loc / ncol(b) ) ]
	loc.x = p0s[ 1 + loc %% ncol(b) ]
	if ( plot ) {
		points( loc.x, loc.y, pch=19 )
	}
	
	## find the p.value corresponding to the passed bound
	if (!is.null( p.value.bound ) ) {
		
		del = abs( bounds - p.value.bound ) < zero.threshold / 5
		del = !is.na(del) & del
		a = alphas
		a[ !del ] = 0
		loc = which.max(a) - 1
		p.value = max(a)
		loc.y = pds[ 1 + floor( loc / ncol(a) ) ]
		loc.x = p0s[ 1 + loc %% ncol(a) ]
		if ( plot ) {
			points( loc.x, loc.y, pch=20, col="red" )
		}
	} else {
		p.value = NULL
	}

		
	list( n=n, k=k, d=d, e.max=e.max, max=max(b), 
			p=c(loc.x, loc.y, 1-loc.x-loc.y),
			p.value=p.value ) 
}






## SIMULATION FUNCTION
## Given a matrix of votes, calculate the weights for all precincts
## and do a sample.  
## Then, if sim.audit=TRUE, assume that p_d percent of the precincts
## (at random) have error, and the errors are due to vote miscounts of
## size 'swing'.
## Param  n: Sample size to draw.
##     sim.audit: pretend that auditing occurs at passed rates.
## Return:   List of taints found in such a circumstance OR precincts selected
##      with relevant attributes (including simulated errors, if asked) OR 
##      the number of non-zero taints and the size of largest taint.
tri.audit.sim = function( Z, n, p_d = 0.1, 
						swing = 5, 
						return.type = c("statistics","taints","precinct"),
						seed=NULL,
						PID = "PID",
						... )    {
	stopifnot( !is.null( Z$V[[PID]] ) )
	stopifnot( p_d > 0 && p_d < 1 )
	
	return.type = match.arg( return.type )	

	# make the truth.
	Z$V$swing.ep = ifelse( runif(Z$N) < p_d, 
			calc.pairwise.e_p( Z, Z$V, err.override=swing ), 0 )
	
	aud = tri.sample( Z, n, PID=PID, simplify=(return.type=="precinct"), ... )

	aud$taint =  aud$swing.ep / aud$e.max
	
	switch( return.type,
		"precinct" = aud,
		"statistics" = c( k=sum( aud$taint > 0 ), d=max( aud$taint ) ),
		"taints" = aud$taint )
}


## Generate sample.
##
## Param: n -- number of samples desired OR an audit.plan.tri object
## Return: List of precincts to audit (with weights--the number of times that
##        precinct was selected).
tri.sample = function( Z, n, seed=NULL, 
						print.trail=FALSE,
						simplify=TRUE,
						return.precincts = TRUE, 
						PID = "PID",
						known="known" ) {
	
	if ( is.audit.plan.tri( n ) ) {
		Z$V$e.max = switch( n$bound,
		 	"passed" = {
				if ( is.null( Z$V$e.max ) ) {
					stop( "No e.max column, and told to not calculate it" )
				}
				Z$V$e.max },
			"e.plus" = maximumMarginBound( Z ),
			"WPM" = 0.4 * Z$V[[Z$tot.votes.col]] / Z$margin
			)
		n = n$n 
	} else {
		if ( is.null(Z$V$e.max) ) {
			warning( "Computing e.max values to sample from" )
			Z$V$e.max = maximumMarginBound( Z )
		}
	}	
		
	stopifnot( n > 0 && n <= Z$N )
	
		
	if ( !is.null(seed) ) {
		cat( "setting seed to ", seed )
		set.seed( seed )	
	}
	
	stopifnot( !is.null( Z$V$e.max ) )
	
	## Note the replace=TRUE--technically, we should sample a SRS of vote units
	## but since we are looking at margins, "vote units" are abstract, arb small 
	## pieces of the margin.  So just sample with replacement.
	sample = sample( 1:nrow(Z$V), n, prob=Z$V$e.max, replace=TRUE )
	
	A = Z$V[sample,]
	if ( simplify ) {
		A = A[ order(A[[PID]]), ]
		A$count = table(A$PID)[ A$PID ]
		A = A[ !duplicated(A$PID), ]
	}
	if ( return.precincts ) {
		A
	} else {
		tri.sample.stats( A )
	}
}


tri.sample.stats = function( samp ) {
	
	c( sampled = nrow( samp ),
		votes = sum( samp$tot.votes ),
		mean.e.max = mean( samp$e.max ) )
		
}

is.audit.plan.tri = function( x ) {
	inherits( x, "audit.plan.tri" )
}








## Calculate estimated sample size to do a trinomial bound that would certify
## assuming a given estimate of low-error error rate
##
## Param: bound -- e.plus, WPM, or use the passed, previously computed, e.max values in the Z object.
tri.calc.sample = function( Z, beta=0.75, guess.N = 20, p_d = 0.1, swing = 5, 
							power=0.90,
							bound = c("e.plus", "WPM", "passed") ) {

	bound = match.arg(bound)
	
	if ( bound == "passed" ) {
		if ( is.null( Z$V$e.max ) ) {
			stop( "No e.max column, and told to not calculate it" )
		}
	} else if ( bound == "e.plus" ) {
		Z$V$e.max = maximumMarginBound( Z )
	} else {
		Z$V$e.max = 0.4 * Z$V[[Z$tot.votes.col]] / Z$margin
	}
	
	## Compute the taints that would be due to swings as large as 
	## passed.
	Z$V$swing.ep = calc.pairwise.e_p( Z, Z$V, err.override=swing )
	Z$V$taint = Z$V$swing.ep / Z$V$e.max
	
	V = Z$V[ order( Z$V$e.max, decreasing=TRUE ), ]
	cumMass = cumsum(V$e.max)
	
	T = sum( V$e.max )
	if ( p_d > 0 ) {
		P = (power^(1/guess.N) - 1 + p_d) / p_d
		cut = sum( cumMass <= P * T)
		d.plus=max(V$taint[1:cut])
	} else { 
		stopifnot( p_d == 0 )
		P = 1
		cut = nrow(V)
		d.plus=0
	}
	
	n = ceiling( log( 1 - beta ) / log( 1 + d.plus * p_d - 1/T ) )
	
	## calculated expected amount of work
	Ep = nrow(V) - sum( (1-V$e.max/T)^n )
	Evts =  sum( V[[Z$tot.votes.col]] * ( 1 - (1-V$e.max/T)^n ) )
	
	res <- list( beta=beta, stages=1, beta1=beta, 
				P=P, cut=cut, swing=swing, d.plus=d.plus, p_d = p_d, T=T, n=n,
				E.p = Ep, E.vts = Evts,
				bound=bound,
				Z = Z )
				
	class(res) <- "audit.plan.tri"
	res	
}




print.audit.plan.tri = function( x, ... ) {
	cat( "Audit plan: beta=", x$beta, "  stages=", x$stages, "  beta1=", x$beta1, 
			"\n\t\td^+=", x$d.plus, " (vote swing of ", x$swing, ")    p_d=", x$p_d, 
			"\n\t\tP=", x$P, "  cut=", x$cut, "  T=", x$T, "  1/T=", (1/x$T),
			"\n\t\tn=", x$n, "  met=", "  PPEB w/ ", x$bound, sep="",
			"\n\t\tE[# pcts audited]=", round( x$E.p, digits=1), "   E[votes audited]=", round(x$E.vts), "\n" )
#			"\n\t\tskipped=", x$skipped, " %mar lost =", (1-x$threshold), "\n" )
#	if ( length(x$stratas) > 1 ) {
#		print( rbind( x$stratas, x$ns ) )
#	}
}



# return the errors of a trinomial audit
trinomial.audit = function( Z, audit ) {

	audit$e.max = maximumMarginBound( Z, audit )	audit = audit.totals.to.OS(Z, audit)	errs = compute.audit.errors( Z, audit=audit, 
					w_p=weight.function("taint"), 
					bound.col="e.max"  )
			
	audit$MOS = errs$err * Z$margin
	audit$taint = errs$err.weighted	k = sum( (audit$taint > 0) * (audit$count) )	d.max = max(audit$taint)

	list( n=sum(audit$count), k=k, d.max=d.max, audit=audit )
}
