# Implementation of CAST system of Stark


## What does n choose k candidates look like?



# Make the cartoon example from the CAST paper as a voter data matrix.
# Param: vote.dist: reported votes for C1, C2, and C3 in order for all precincts.prompt
make.cartoon = function( n = 400, vote.dist=c(125,113,13), stratify=TRUE ) {
	n = n - n %% 4
	
	V = data.frame( strata = rep( c(rep(c("S1","S1.VBM"), 3),"S2","S2.VBM"), n/4 ),
					tot.votes = rep( 255, 2*n ),
					C1 = rep( vote.dist[[1]], 2*n ),
					C2 = rep( vote.dist[[2]], 2*n ),
					C3 = rep( vote.dist[[3]], 2*n ) )
	
				
	if ( !stratify ) {
		V$strata = "S1"
	}
	PID = paste( V$strata, rownames(V), sep="-" )
	names(PID) = "PID"
	V = cbind( PID, V )
	V$PID = as.character(V$PID)
	
	
	make.Z( V, c( "C1", "C2", "C3" ) )
}


## make the bad truth as descibed in Stark's paper.
make.truth.ex.bad = function( Z ) {
	n = nrow(Z$V)
	stopifnot( n %% 8 == 0 )
	C1 = c( rep( 80, n/8 ), rep( 124, 7*n/8 ) )
	C2 = c( rep( 160, n/8 ), rep( 113, 7*n/8 ) )
	C3 = c( rep( 13, n/8 ), rep( 15, 7*n/8 ) )
	reord = sample( 1:n, n )
	
	V = data.frame( strata = rep( c(rep(c("S1","S1.VBM"), 3),"S2","S2.VBM"), n/8 ),
					tot.votes = rep( 255, n ),
					C1=C1[reord], C2=C2[reord], C3=C3[reord] )
			
	PID = paste( V$strata, rownames(V), sep="-" )
	names(PID) = "PID"
	V = cbind( PID, V )
	V$PID = as.character(V$PID)
	
	make.Z( V, c( "C1", "C2", "C3" ) )
}

## Given a collection of counted votes, make a theoretical bad truth by packing all error 
## possible in the largest precincts
##
## Warning: if bound is WPM this error is made by simply adding the max amount of error
## to the first loser's total (so that total votes may in this case exceed the total votes
## of the precinct)--this could potentially cause trouble.  Be careful!
##
## Param t: an allowed backgound level of error for all precincts
##       bound.function: Can pass fractionOfVotesBound if you want to do 0.4WPM.
make.truth.opt.bad = function(Z, strata="strata", 
			bound=c("margin","WPM"), t=0 ) {
	bound = match.arg( bound )
	
	if ( bound=="WPM" ) {
		s = CAST.calc.sample(Z, t=t, bound.function=fractionOfVotesBound )
	} else {
		s = CAST.calc.sample(Z, t=t )
	}
	
	winner = Z$winner[[Z$f]]
	
	stopifnot( all( s$Z$V$PID == Z$V$PID ) )
	V = Z$V[ order(s$Z$V$e.max, decreasing=TRUE), ]
	baddies = 1:nrow(V) <= s$q
	V[winner] = pmax( 0,  V[[winner]] - t )
	
	tots = V$tot.votes[baddies]
	if ( bound=="margin" ) {
		newvotes = as.data.frame( matrix( 0, nrow=s$q, ncol=length(Z$C.names) ) )
		names(newvotes) = Z$C.names
		newvotes[[ Z$losers[[1]] ]] = tots
		V[ baddies, Z$C.names ] = newvotes 
	} else {
		V[ baddies, Z$losers[[1]] ] = V[ baddies, Z$losers[[1]] ] + round( 0.4 * V[ baddies, ]$tot.votes )
	}

	
	nZ = make.Z( V[ Z$V$PID, ], Z$C.names )
	stopifnot( nZ$winner != Z$winner )
	nZ$num.tweaked = s$q
	nZ
	
}

## make bad truth as described in Stark's paper (assuming fixed precinct size)
make.truth.opt.bad.strat = function(Z, strata="strata", t=3, shuffle.strata=FALSE) {
	s = CAST.calc.sample(Z, t=t)
	
	if ( shuffle.strata ) {
		pids = sample(Z$V$PID, s$q )
	} else {
		strats = split( Z$V$PID, Z$V[strata] )
		pids = lapply( names(s$stratas), function( ST ) {
			sample( strats[[ST]], floor( s$q * s$stratas[[ST]] / s$N ) )
		} )
		pids = unlist( pids )
		if ( ( ext <- ( s$q - length( pids ) ) ) > 0 ) {
			pids = c( pids, sample( setdiff( Z$V$PID, pids ), ext ) )
		}
	}
		
	Z$V["C1"] = Z$V["C1"] - t
	baddies = Z$V$PID %in% pids
	Z$V[ baddies, Z$C.names ] = cbind( rep( 0, s$q ), rep(255,s$q), rep(0,s$q) ) 
	make.Z( Z$V, Z$C.names )
}

make.ok.truth = function( Z, num.off = 8, amount.off = 5 ) {
	cut = floor(num.off / 2)
	off = sample( 1:nrow(Z$V), num.off )
	offUp = off[1:cut]
	offDown = off[(cut+1):num.off]
	CW = Z$winners[Z$f]
	CL = Z$losers[1]
	maxL = pmin( Z$V[ offUp, CL ], amount.off )
	Z$V[ offUp, CL ] = Z$V[ offUp, CL ] - maxL
	Z$V[ offUp, CW ] = Z$V[ offUp, CW ] + maxL
	maxL = pmin( Z$V[ offDown, CW ], amount.off )
	Z$V[ offDown, CW ] = Z$V[ offDown, CW ] - maxL
	Z$V[ offDown, CL ] = Z$V[ offDown, CL ] + maxL
	
	Z = countVotes(Z)
	
	Z
}


## Make a Z matrix with the desired charactaristics for simulation studies
## and debugging purposes.
## 
## Param M: margin
##       N: number of precincts
##       strata: number of strata
##       per.winner: percent of vote winner got (to get undervotes, 3rd cand, etc.)
## Return: Z voter matrix with desired charactaristics.
make.sample = function( M, N, strata=1, per.winner=NULL, 
						worst.e.max = NULL, 
						R = NULL, 
						tot.votes=100000 ) {
	if ( !is.null( per.winner ) ) {
		v = c( per.winner, per.winner - M, 1 - 2 * per.winner + M )
		names(v) = c("WNR", "LSR", "OTR" )
		stopifnot( min(v) >= 0 && max(v) < 1 )
		stopifnot( v[[1]] == max(v) )
	} else {
		stopifnot( M < 1 && M > 0 )
		v = c( (1+M) / 2, (1-M)/2 )
		names(v) = c("WNR", "LSR" )
	}
	
	st = paste( "ST", 1:strata, sep="-" )
	V = round( (tot.votes / N) * matrix( rep( v, N ), ncol=length(v), byrow=TRUE ) )
	V = data.frame( V )
	V = cbind( paste("P", 1:N, sep="-"), rep( st, length.out=N ), apply( V, 1, sum ), V )
	V[[1]] = as.character( V[[1]] )
	names(V) = c( "PID", "strata", "tot.votes", names(v) )

	if ( !is.null( worst.e.max ) ) {
		# calc how bad worst prec is, and then use p(x) = c/sqrt(x)
		# as density function of precincts.  Invert the CDF, and eval
		# at orders (1:N)/N to get the new precinct sizes.  Scale
		# and go!
		t1 = V[1, ]
		e.max = t1$tot.votes - t1$LSR + t1$WNR
		stopifnot( is.null( R ) )
		R = worst.e.max / e.max
	}
	if ( !is.null(R) ) {   # using the ratio of bad to baseline, go.
		stopifnot( R >= 0 )
		if ( R <= 2 ) {
			m = R / N
			wts = 1:N * m - N*m / 2 + 1
		} else {
         p = (R - 2)/(R - 1)
         C = (1 - p)/R^(1 - p)
        
       	 wts = (((1:N)/N)/(C * (R - 1)))^(R - 1)
		}

		V[names(v)] = round( V[names(v)] * wts )
		V["tot.votes"] = apply( V[names(v)], 1, sum )
	}
	
	make.Z( V, C.names=names(v) )
}


## Given a vector of precinct totals and the total votes for the winner
## and the loser, make a plausible precinct-by-precinct vote count that
## works. 
## Note: the margins of the precincts will all be the same as the margin
## of the overall race.
make.sample.from.totals = function( vote.W, vote.L, totals ) {
    GT = sum(totals)
    v = c(vote.W/GT, vote.L/GT)
    names(v) = c("WNR", "LSR")
    stopifnot(min(v) >= 0 && max(v) < 1)
    stopifnot(v[[1]] == max(v))
    N = length(totals)
    V = matrix(rep(v, N), ncol = length(v), byrow = TRUE)
    V = data.frame(V)
    V = cbind(rep("ST-1", length.out = N), totals, V)
    names(V) = c("strata", "tot.votes", names(v))
    V[names(v)] = round(V[names(v)] * totals)

totW = vote.W - sum( V$WNR )
totL = vote.L - sum( V$LSR )

# make total match (gets off due to rounding)
tweak = function( votes, flex, tot ) {
	cntr = 1
	while ( tot != 0 && cntr <= length(votes) ) {
		if ( tot > 0 ) {
			if ( flex[cntr] > 0 ) {
				votes[cntr] = votes[cntr] + 1
				tot = tot - 1
			}
		} else {
			if ( votes[cntr] > 0 ) {
				votes[cntr] = votes[cntr] - 1
				tot = tot + 1
			}
		}
		cntr = cntr + 1
	}
	votes
}

	
flex = with( V, tot.votes - WNR - LSR )
V$WNR = tweak( V$WNR, flex, totW )
flex = with( V, tot.votes - WNR - LSR )
V$LSR = tweak( V$LSR, flex, totL )

	Z =    make.Z(V, C.names = names(v))
	stopifnot( Z$total.votes == sum( totals ) )
	stopifnot( sum(Z$V$WNR) == vote.W )
	stopifnot( sum(Z$V$LSR) == vote.L )
	
	Z
}


	
	

make.sample.from.totals.margin = function( M,  totals, per.winner=NULL ) {
	GT = sum(totals)
	vote.W = round( GT* (0.5 + M/2) )
	vote.L = round( GT * (0.5 - M/2) )
	make.sample.from.totals( vote.W, vote.L, totals )
}



## Make an audit.plan given reported results for an election.  It gives back what to
## do for a single stage.  If stages is > 1, then it adjusts beta appropriately.
##
## param Z      : voter matrix
##       beta   : overall chance of correctly escalating a bad election to full recount
##       stages : number of auditing stages before full recount
##       t      : Threshold error for escalation -- if >= 1 then # votes, otherwise
##                fraction of margin.
##  small.cut   : Consider all precincts smaller than this to be as wrong as possible
##                and then remove them from the N to sample from.
##        drop  : Don't audit precincts in the drop column (of T/F) since they are known
##     as.taint : TRUE means interpret t as a taint in [0,1] by batch.  FALSE is a prop of
##                margin or number of votes (as desc. above)
CAST.calc.sample = function( Z, beta = 0.9, stages=1, t=3, as.taint=FALSE,
				small.cut = NULL,
				strata=NULL, drop=NULL, 
				method=c("select", "binomial","hypergeometric"),
				calc.e.max = TRUE,
				bound.function = maximumMarginBound ) {

	method=match.arg(method)
	
	if ( !is.null(strata) && is.null( Z$V[[strata]] ) ) {
		warning( "No strata column of name '", strata, "'--assuming single strata" )
		strata=NULL
	}
	
	if ( is.null( drop ) ) {
		Z$V$known = FALSE
		drop="known"
	}
	
	beta1 = beta^(1/stages)
	
	# convert t to fraction of margin
	if ( t >= 1 ) { t = t / Z$margin }
	
	# Ignore small precincts?  Even if not,
	# there is no point in auditing precincts smaller than
	# t.
	if ( as.taint ) {
		# do nothing---in the taint world, nothing is too small (unless
		# it is an empty precint)
		small.cut = 0
	} else if ( is.null( small.cut ) ) {
		small.cut = t
	} else if ( small.cut >= 1 ) {
		small.cut = small.cut / Z$margin
	}
	
	if ( calc.e.max ) {
		Z$V$e.max = bound.function( Z )
	}
	
	Z$V$small = Z$V$e.max <= small.cut
	
	# do not consider all precincts that are either known or too small to pay
	# attention to.
	skip = Z$V[[drop]] | Z$V$small
	Z$V$skip = skip
	
	if ( is.null(strata) ) {
		strat.size = sum( !skip )
	} else {
		strat.size = table( Z$V[ !skip, strata ] )
	}
	
	# Reduce the threshold by all non-known precincts that we are skipping,
	# as we must assume the error in them is as large as possible.
	# Also calc new N (number of remaining batches to select from).
	threshold = 1 - sum( Z$V$e.max[ !Z$V[[drop]] & Z$V$small ] )
	N = Z$N - sum(skip)
	
	if ( threshold <= 0 ) {
		# We are doomed.  Must audit everything.
		q = -1
		n = Z$N
		
	} else {	
		if ( as.taint ) {
			q = find.q( Z$V, drop=skip, t, "e.max", threshold=threshold,
							w_p = weight.function("taint") )
		} else {
			q = find.q( Z$V, drop=skip, t, Z$tot.votes.col, threshold=threshold )
		}
		
		stopifnot( N >= q )
	
	
		if ( method=="select") {
			method = ifelse( length(strat.size)==1, "hypergeometric","binomial")
		}
	
		if ( q == 0 ) {
			n = N
		} else if ( method=="hypergeometric" ) {
			n = ceiling( uniroot( function(x) { 
				lchoose( N - q, round(x) ) - lchoose( N, round(x) ) - log(1-beta1) }, 
								  c(0,N-q) )$root )
		} else {
			n = ceiling( log( 1-beta1 ) / log( ( N - q ) / N ) )
		}
	}
	
	if ( n > N ) {
		n = N
	}

	ns = ceiling( n * strat.size / N )
	
	# calculate estimated work
	V2 = Z$V[ !skip, ]
	if ( is.null( strata ) ) {
		E.votes = mean(V2$tot.votes) * ns
	} else {
		mns = tapply( V2$tot.votes, V2[[strata]], mean )
		E.votes = ns * mns[names(ns)]
	}
	
	res = list( beta=beta, beta1 = beta1, stages=stages, t=t, as.taint=as.taint, 
				q=q, N=N, n=n, 
				stratas=strat.size, strata=strata, ns=ns, method=method, Z=Z, 
				threshold=threshold, skipped = sum(skip), 
				bound.function=bound.function,
				E.votes=E.votes )
	class(res) = "audit.plan"
	res
}


CAST.calc.opt.cut = function( Z, beta = 0.9, stages=2, t=3, plot=FALSE, ... ) {

	cuts = t:max(2*Z$V$tot.votes)
	ns = rep( NA, length(cuts) )
	qs = rep( NA, length(cuts) )
	cnt = 1
	done = FALSE
	while( !done ) {
		
		s = CAST.calc.sample( Z, beta, stages, t, small.cut=cuts[[cnt]], ... )
		if ( s$q < 0 ) {
			done = TRUE
		} else {
			ns[cnt] = s$n
			qs[cnt] = s$q
		#	cat( "on ", cuts[[cnt]], " - ", ns[cnt], "\n" )
			cnt = cnt+1
		}
	}
	names(ns) = cuts
	ns = ns[1:(cnt-1)]
	qs = qs[1:(cnt-1)]
	cuts = cuts[1:(cnt-1)]
	
	if ( plot ) {
		plot( names(ns), ns, ylim=c(0,max(ns)), type="s", main="Sample Size by Small Cut", pch=19,
		ylab="sample size needed", xlab="cut-off for cutting out small precincts" )
		scale = max(qs) / max(ns) 
		qsU = unique( qs )
		points( names(ns), qs / scale, type="s", col="red" )
		axis( side=4, at=qsU / scale, labels=qsU )
	}
	cuts = data.frame( cut=cuts, n=ns, q=qs )
	v = which.min( cuts$n )
	cuts[v, ]
}


is.audit.plan = function( x ) {
	inherits( x, "audit.plan" )
}


## Print a nice pretty audit plan.
print.audit.plan = function( x, ... ) {
  P = x
  if ( x$as.taint ) {
  	met = paste( P$method, "(taint)", collapse="")
  } else {
  	met = P$method
  }
  cat(sprintf("Audit plan: beta=%.2f  stages=%d   beta1=%.2f   met=%s\n", 
  						P$beta, P$stages, P$beta1, met))
  if ( P$t < 1 ) {
  	cat(sprintf("\tt=%.3f\t q=%d\t N=%d / %d\t n=%d\n", 
  						P$t, P$q, P$N, P$Z$N, P$n ) )
  } else {
  	cat(sprintf("\tt=%d\t q=%d\t N=%d / %d\t n=%d\n", 
  						P$t, P$q, P$N, P$Z$N, P$n ) )
  }
  cat(sprintf("\tskipped=%d\t mar lost=%.0f%%\n",
  						P$skipped, 100*(1-P$threshold) ) )
  cat(sprintf("\tE[# pcts audited]=%.1f\t\t E[votes audited]=%.1f\n",
  						P$n, sum(P$E.votes) ) )
	if ( length(P$stratas) > 1 ) {
		tb =  rbind( P$stratas, P$ns, P$E.votes )
		rownames(tb) = c( "N", "n", "E.vts" )
		colnames(tb) = names(P$ns)
		print(tb)
	}
}
				

## Sample from the various strata according to the schedule set by 'ns'.
## Ignore all precincts that are known (i.e., have been previously audited).
##
## Param known: name of column of true/false of known vs. unknown counts
##
## Return: List of precincts to be audited.
CAST.sample = function( Z,  ns,  strata=NULL, seed=NULL, 
						print.trail=FALSE, known="known" ) {
	pt = function( ... ) { if ( print.trail ) { cat( ..., "\n" ) } }
		
	if ( is.audit.plan(ns) ) {
		strata = ns$strata
		ns = ns$ns
	}
	
	if ( !is.null( seed ) ) {
          set.seed( seed )
          pt( "Setting seed to ", seed )
	}
	
	if ( !is.null( Z$V[[known]] ) ) {
		V = Z$V[ !Z$V$known, ]
	} else {
		V = Z$V
	}

	if ( is.null( strata ) ) {
		stopifnot( length(ns) == 1 )
	
		ids = sample( nrow(V), ns )
		pt( "Full sample: ", ids )
		
		V$PID[ids]

	} else {
		strat = split( V, V[[strata]] )
	
		l = lapply( 1:length(ns), function( X ) {
			st = names(ns)[[X]]
			pids = sort( strat[[ st ]]$PID )
			ids = sample( length(pids), ns[[X]] )
			pt( "Strata", st, ": ", ids )
		
			pids[ids]
		} )
		unlist( l )
	}
}


## Given audit data, comptue p.values and all that.
CAST.audit = function( Z, audit=NULL, plan=NULL, ... ) {
	if ( is.null( plan ) ) {
		plan = CAST.calc.sample( Z, ... )
	}

	if ( is.null(audit) ) {
		plan$Z$audit = Z$audit
	} else {
		plan$Z$audit = audit
	}
		
	res = stark.test.Z( plan$Z, drop="skip", max_err=plan$bound.function,
			bound.col=Z$tot.votes.col,
			strat.col=plan$strata )
	res$certify = res$p.value < (1 - plan$beta1)
	
	cat( "Certify election: ", res$certify, "\n" )
	res
}



## Given the reported vote table, Z, and the actual truth (simulated)
## (a Z matrix with same precincts), and a list of precincts to audit,
## do the audit.  If audit.names
## is null and the ns is not null, it will sample from precincts via
## CAST.sample automatically.
##
## Return:   Overstatments for each candidate for each precinct.
do.audit = function( Z, truth, audit.names, ns=NULL ) {
	
	if ( is.null( audit.names ) ) {
		audit.names = CAST.sample( ns )
	}
	rownames(truth$V) = truth$V$PID
	truth$V = truth$V[ Z$V$PID, ]  # get in same order
	
	os = truth$V
	os[Z$C.names] =  Z$V[Z$C.names] - os[Z$C.names]
	os$clear = apply( os[Z$C.names], 1, function(X) { sum( abs(X) ) == 0 } )
	os = os[ os$PID %in% audit.names, ]
	
	os
}


## Single test
test.CAST = function() {
	Z = make.cartoon()
	
	samp.info = CAST.calc.sample( Z )
	samp.info
	
	samp = CAST.sample( Z, samp.info )
	
	list( samp.info, samp )
}



# testing the CAST system
# return: Number of stages audited (s+1=full recount) before stopping
sim.race = function( n=800, beta=0.75, stages=2, 
						truth.maker=make.truth.opt.bad,
						print.trail=FALSE) {

	Z = make.cartoon(n=n)
	Z$V$known = FALSE
	rownames(Z$V) = Z$V$PID
	Z
	
	truth = truth.maker(Z)
	rownames(truth$V) = truth$V$PID
	truth
			
	s = 0
	tots = rep( 0, stages )

	while ( s < stages && (s == 0 || t > samp.info$t) ) {
		s = s + 1
		samp.info = CAST.calc.sample( Z, stages=stages, beta=beta, drop="known" )
		
		tots[s] = sum( samp.info$ns )
		
		## add entropy and then use this function to generate table of audits
		audit.names = CAST.sample( Z, samp.info, 
							print.trail=print.trail, known="known" )
		aud = do.audit( Z, truth, audit.names )
		Z$audit = aud
		
		t = compute.stark.t( Z, "tot.votes" )
		
		if ( t > samp.info$t ) { # escalate!
			# replace count information with the 'true' audit
			# information for all audited precincts.
			Z$V[audit.names,] = truth$V[audit.names,]
			# mark these precincts as known so that future stages will
			# not audit again.
			Z$V[ audit.names, "known" ] = TRUE
			# rebuild Z to update margin totals, etc.
			Z = make.Z( Z$V, Z$C.names )  
		}
	}
	c( s, sum(tots), t < samp.info$t )
}



