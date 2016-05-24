### election functions

### Testing using simple random sample of precincts or stratified SRS of precincts.
### Miratrix Mar, 08
### From 4.1 - Stark

### The various weight functions, packaged in a single object
weight.function = function( name=c("no.weight","weight","weight.and.slop", "margin.weight","taint") ) {
  name = match.arg( name )
  switch( name,
     
         ## no weighting
         no.weight= c( function( x, b_p, M ) { x + 0 * b_p }, 
           function( x, b_p, M ) { x + 0 * b_p } ),
         
         ## weighted by size of precint
         weight = c( function( x, b_p, M ) { x / b_p }, 
           function( x, b_p, M ) { x * b_p } ),
         ## weight by size, after a slop of 2 votes has been taken off
         weight.and.slop = c( function( x, b_p, M ) { pmax( x - 2, 0 ) / b_p },
           function( t, b_p, M ) { t * b_p + 2 } ),
         ## for pairwise margin tests
         margin.weight = c( function( x, b_p, M ) { pmax(x - 2/M, 0) / b_p },
           function( x, b_p, M ) { x * b_p + 2/M } ),
         ## taint -- divide by maximum error in the precinct
         taint = c( function( x, b_p, M ) { x / b_p },
            function( x, b_p, M ) { x * b_p } )
           )
}


elec.data = function( V, C.names=names(V)[2:length(V)], f = 1, 
  audit=NULL, pool=TRUE, tot.votes.col="tot.votes", PID.col="PID" ) {
  ## Make the 'Z' matrix that holds all the vote totals, etc., as well as some
  ## precomputed information such as vote margins between candidates, the theoretical
  ## winners, and so on.
  ##
  ## make.Z does some cleaning and renaming of the passed data structure.  In particular
  ## it will rename the tot.votes column to "tot.votes" if it is not that name already.
  ##
  ## Param    V: Voter matrix OR 2-element list with Voter Matrix followed 
  ##             by Candidate names
  ##       pool: Combine small candidates into single pseudo-candidates to increase 
  ##             power
  ##
  ## Return:  A "elec.data" data structure. 
  ##         Note: Will _add_ PID (precinct ID) column (and generate PIDs)
  ##               if no PID provided.  Also, rownames
  ##               are always PIDs (so indexing by PID works)




  Z = list()
  class( Z ) = "elec.data"
  
  if ( length(V) == 2 ) {
    C.names = V[[2]]
    V = V[[1]]
  }
  
  Z$N = length(V[[1]])     # number of precincts
  Z$C = length(C.names)    # number of candidates
  Z$f = f                  # number of possible winners

  Z$V = V
  Z$C.names = C.names

  ##	Z$alpha = 0.10
  ##	Z$a_s = rep(0.5,10)^(1:10)

  stopifnot( all( Z$C.names %in% names(Z$V) ) )
  
  if ( is.numeric(tot.votes.col) ) {
  	tot.votes.col = names(V)[[tot.votes.col]]
  }
  if ( is.na(l <- match( tot.votes.col, names(V) ) ) ) {
    warning( "No tot.votes.col found in vote matrix--will be recounting" )
  } else {
  	names(Z$V)[l] = "tot.votes"
  }
  Z$tot.votes.col = "tot.votes"

  Z = countVotes(Z)

  if( !is.null( audit ) ) {
    stopifnot( all( Z$C.names %in% names(audit) ) )
    Z$audit = audit
  }

  ## Pooling!
  ln = length(Z$losers)
  if ( pool && 
      ln > 2 && ( (Z$totals[Z$losers[ln]] + Z$totals[Z$losers[ln-1]]) < Z$totals[Z$losers[1]] )) {
    
    grb = ln - 1
    while( grb > 1 &&  sum( Z$totals[Z$losers[grb:ln]] ) < Z$totals[Z$losers[1]] ) {
      grb = grb - 1
    }
    grb = grb + 1
    to.pool = Z$losers[grb:ln]
    
    Z$V[["pool"]] = apply(Z$V[to.pool],1,sum)
    Z$V = Z$V[setdiff(names(Z$V),to.pool)]
    if ( !is.null( audit )) {  # pool the audit info too.
      Z$audit[["pool"]] = apply(Z$audit[to.pool],1,sum)
      Z$audit = Z$audit[setdiff(names(Z$audit),to.pool)]
    }
    Z$C.names = c( setdiff(Z$C.names, to.pool ), "pool" )
    
    Z$C = length(Z$C.names)    # number of candidates
    Z = countVotes(Z)
    ln = ln - 1

  }
  
  
  if ( is.numeric(PID.col) ) {
  	names(Z$V)[PID.col] = "PID"
  } else if ( is.na(l <- match( PID.col, names(V) ) ) ) {
    warning( "No PID.col found in vote matrix--will be generated" )
  } else {
  	names(Z$V)[l] = "PID"
  }
  Z$PID.col = "PID"

  if ( is.null( Z$V$PID ) ) {
    PID = paste( "P", rownames(V), sep="-" )
    names(PID) = "PID"
    Z$V = cbind( PID, Z$V )
  }	
  Z$V$PID = as.character(Z$V$PID)
  stopifnot( sum(duplicated(Z$V$PID)) == 0 )
  rownames(Z$V) = Z$V$PID
  
  Z
}

                                        # this is what make.Z should be called.  Should fix.
make.Z = function( ... ) {
  elec.data( ... )
}

is.elec.data = function( x ) {
	inherits(x, "elec.data")
}

print.elec.data = function( x, n=4, ... ) {
                                        # "N"           "C"           "f"           "V"           "C.names"     "total.votes" "margin"      "margin.per"  "totals"      "winners"     "losers"      "audit"       "Ms"         
  Z = x
  cat( "Z frame:  N = ", Z$N, "\tC/f = ", Z$C, "/", Z$f, "\ttotal votes = ", 
      Z$total.votes, "    M = ", Z$margin,
      " (", round(100*Z$margin.per), "%)\n",
      "\t\tNames = ", paste(Z$C.names, " (", round(100*Z$totals/Z$total.votes), "%)", sep="", collapse=", "), "\t  Winners = ", 
      paste( Z$winners, collapse=", " ), "    Losers = ", 
      paste( Z$losers, collapse=", "), "\n", sep="")
  if ( Z$C > 2 ) { 
    cat( "Pairwise margins:\n" )
    print( Z$Ms )
  }
  
  cat( "Sample votes (", length( rownames(Z$V) ), " records)\n", sep="" )
  print( head( Z$V, n ) , fill=TRUE, labels=c("\t"))
  if ( !is.null(Z$audit ) ) {
    cat( "Sample Audits (", length( rownames(Z$audit) ), " records)\n", sep="" )
    print(	head( Z$audit, n ) )
    
                                        #fill=TRUE, labels=c("\t")
  } else {
    cat( "(No audit information)\n" )
  }
  
  invisible( Z )

}







maximumMarginBound = function( Z, votes=NULL ) {
## return the maximum margin reduction for each precint by computing
## all margin reductions between pairs of winners & losers and then
## scaling by that pair's total margin to get a proportion and then
## taking the max of all such proportions (usually will be the last 
## winner to the closest loser).
## Return:   Vector (of length of precincts) of maximum possible error for 
##           each precinct.
    
  if ( is.null( Z$winners ) ) {
    stop( "Need to count votes before computing bouns in maximumMarginBound()" )
  }
  if ( is.null( Z$Ms ) ) {
    stop( "Need to compute pairwise margins before computing error bounds." )
  }
  if ( is.null( Z$V[[Z$tot.votes.col]] ) ) {
    stop( "Vote total column, ", Z$tot.votes.col, " not found." )
  }
  
  
  err.bound = data.frame( row.names = row.names(Z$V) )
  
  for ( i in Z$winners ) {
    for ( j in Z$losers ) {	
      mrg = paste(i,j,sep="_")
      err.bound[[mrg]] = (Z$V[[Z$tot.votes.col]] + Z$V[[i]] - Z$V[[j]])/Z$Ms[[i,j]]		
    }
  }

  res = apply( err.bound, 1, max )
  
  if ( !is.null( votes ) ) {
  	stopifnot( !is.null( votes$PID ) && is.character( votes$PID ) )
  	res[ votes$PID ]
  } else {
  	res
  }
  
}


fractionOfVotesBound = function( Z, frac=0.4 ) {
  ## This is the 0.4b bound function.  It returns frac * total votes.
  ## Return:   Vector (of length of precincts) of maximum error for 
  ##           each precinct.
  
  Z$V[[Z$tot.votes.col]] * frac / Z$margin
}


countVotes = function( Z ) {
  ## Count the total votes for various candidates.
  ## Return:  Updated 'Z' matrix with the total votes as components
  ##          inside it.
  
  
    computeMargins = function( Z ) {
    ## Used by make.Z
    ##
    ## Return: the pairwise margins (as # of votes) as a data.frame with the rows being
    ## the winners and the columns being the losers.
    
    if ( is.null( Z$winners ) ) {
      stop( "Need to count votes before computing margins in computeMargins()" )
    }
    Ms = data.frame( row.names = Z$winners )
    
    for ( j in Z$losers ) {
      Ms[[j]] = rep(0, length(Z$winners))	
      
      for ( i in Z$winners ) {
        
        Ms[[ i, j ]] = Z$totals[[i]] - Z$totals[[j]]
        
      }
    }

    Ms
  }


  ## Recalculate N based on V
  Z$N = nrow( Z$V )
                                        # total votes in given precinct
  if ( is.null( Z$V[[Z$tot.votes.col]] ) ) {
    warning("Totalling votes from candidates in countVotes()\n" )
    Z$V[[Z$tot.votes.col]] = apply( Z$V[Z$C.names], 1, sum )
  }
  grand.tot = sum(Z$V[[Z$tot.votes.col]])
  
  tots = sapply(Z$V[Z$C.names],sum)
  nt = sort(tots, decreasing=TRUE)
  M = nt[Z$f] - nt[Z$f+1]   # smallest margin (last winner - best loser)
  

  winners = names(nt[1:Z$f])
  losers = setdiff( names(nt), winners )

  Z[c("total.votes", "margin", "margin.per", "totals","winners","losers")] = list(grand.tot, M, M/grand.tot, tots, winners, losers )
  
  ## Compute pairwise margins between all winners and losers
  Z$Ms = computeMargins( Z )

  Z
  #c( Z, list( margin=M, totals=tots, winners=winners, losers=losers ) )
}


## Functions to calculate error amount in audit sample ##

calc.overstatement.e_p = function( Z ) {
  ## One way of calculating the errors for a collection of audited precints.
  ## This one is the sum of all winner overcounts plus the sum of all 
  ## loser undercounts (for each precinct)
  ## Return: Vector (of length of audited precincts) of found errors by precinct. 
  
  apply( Z$audit, 1, 
        function(p) {sum( pmax(p[Z$winners],0) ) + sum( pmax(-1*p[Z$losers],0) )}  )
}


calc.pairwise.e_p = function( Z, audit=NULL, err.override=NULL ) {
  ## Calculate the error by finding the maximum margin reduction for each precint
  ## by computing all margin reductions between pairs of winners & losers
  ## (scaling by that pair's total margin to get a proportion) and then
  ## taking the max of all such proportions (usually will be the last 
  ## winner to the closest loser).
  ## Param err.override:  If non-null, use this as the found error in votes rather than
  ##                      the actual errors found in the audit.
  ## Return: Vector (of length of audited precincts) of found errors by precinct. 

 
  if ( is.null( audit ) ) {
 	audit = Z$audit
  }
  if ( is.null(audit) ) {
    stop( "No audit to calculate e_p values for." )
  }

  stopifnot( !is.null(audit$PID ) )
  rownames(audit) = audit$PID
  stopifnot( Z$C.names %in% names(audit) )
  
  # If there is an error override, we need to make sure not to
  # give more error to a precinct than it can hold. Thus we need
  # a bound on maximum error--hopefully it is passed.
  if ( !is.null( err.override ) ) {
  	if ( is.null( audit$e.max ) ) {
   		if ( is.null( audit[[Z$tot.votes.col]] ) ) {
		    warning( "No '", Z$tot.votes.col, "' or e.max in audit matrix in calc.pairwise.e_p" )
        	audit$e.max = Inf
        } else {
	        warning( "No e.max column in calc.pairwise.e_p with a set error override.  Using tot votes / M." )
    		audit$e.max = audit[[Z$tot.votes.col]] / Z$M
        }
     }
  }
	
  err.bound = data.frame( row.names = row.names(audit) )
  
  for ( i in Z$winners ) {
    for ( j in Z$losers ) {
      
      mrg = paste(i,j,sep="_")

      if ( is.null(err.override) ) {
        err.bound[[mrg]] =  (audit[[i]]-audit[[j]]) / Z$Ms[[i,j]]
      } else {
        err.bound[[mrg]] = pmin( audit$e.max, err.override/Z$Ms[[i,j]] )
      
      }
    }
  }

  apply( err.bound, 1, max )
}


## Calculate the measured error in each of the audited precicnts.
##
## Param bound.col:     This is the vector (in audit) containing the maximum number
##                      of votes (or error) possible in the various precincts.
## Param err.override:  If non-null, use this as the found error in votes rather than
##                      the actual errors found in the audit.
##
## Return:    Orig audit table from Z with two new columns, err and err.weighted, 
##            corresponding
##            to the errors found in each audited precinct before and after the 
##            weight function has been applied to them.
compute.audit.errors = function( Z, audit=NULL,
  calc.e_p=calc.pairwise.e_p,
  w_p = weight.function("no.weight"),
  bound.col="tot.votes",
  err.override = NULL ) {
  
  if ( !is.null(audit) ) {
  	Z$audit = audit
  }
                                        # calculate the errors
  Z$audit$err = calc.e_p( Z, err.override=err.override )
  
                                        # weight the errors
  Z$audit$err.weighted = w_p[[1]]( Z$audit$err, Z$audit[[bound.col]], Z$margin )

  Z$audit
}						

compute.stark.t = function( Z,
  bound.col,
  calc.e_p=calc.pairwise.e_p,
  w_p = weight.function("no.weight"),
  err.override = NULL,
  return.revised.audit = FALSE
  ) {
  ## Compute the error statistic given audit data and relevant functions
  ##
  ## Param Z              If it already has an audit table with err and err.weighted
  ##                      then it will use those errors, otherwise it will compute them
  ##                      with compute.stark.err
  ##
  ## Param bound.col:     This is the vector containing the maximum number of votes
  ##                      possible in the various precincts.
  ## Param err.override:  If non-null, use this as the found error in votes rather than
  ##                      the actual errors found in the audit.

  ## Return: the test statistic, i.e. the maximum found error in the audit 
  ##         sample, as computed by calc.e_p and weighted by w_p.
  
  if ( is.null( Z$audit[[bound.col]] ) ) {
    warning( "Assuming audit's rownames are unique and correspond to Z$V's row names" )
    Z$audit[[bound.col]] = Z$V[rownames(Z$audit), bound.col]
  }
  if ( is.null( Z$audit$err ) ) {
    Z$audit = compute.audit.errors( Z, bound.col=bound.col, 
                   calc.e_p=calc.e_p, w_p=w_p, err.override=err.override )
  }
  
  t.stat = max( Z$audit$err.weighted )
  
  if ( return.revised.audit ) {
    list( t.stat, Z$audit )
  } else {
    t.stat
  }
}

find.q = function( V, t.stat, bound.col, M, threshold=1.0,
			 w_p = weight.function("no.weight"),
 			 drop=NULL ) {
  ## Find q, the minimum number of precints with w_p's greater than given t.stat
  ## that can hold an entire election shift in them.

  ## I.e., find the number of precints that need to have "large taint" in order to
  ## flip the election.  This is, essentially, finding a collection of precints
  ## such that the max error (e.max) plus the background error (the w_p-inverse of the
  ## t.stat) for the rest of the precints is greater than the margin (or 1 if done
  ## by proportions).
  ##
  ##
  ## Param bound.col:  The name of the column in V to be used for the passed size 
  ##                   (max # votes, total votes, incl undervotes, etc.) to the error 
  ##                   function.
  ##       threshold:  The total amount of error to pack in the set of tainted precincts
  ##            drop:  Drop precincts with this column having a "true" value--they are
  ##                   previously audited or otherwise known, and thus can't hold error.
  ##                   Can also pass a logical T/F vector of the length of nrow(V)
  ##
  ## Return:  integer, number of badly tainted precints needed to hold 'threshold' error
  
  stopifnot( is.null( V$e.max ) == FALSE )
  
  if ( !is.null( drop ) ) {
    if ( length(drop) == nrow(V) ) {
      V = V[ !drop, ]
    } else {
      V = V[ !V[[drop]], ]
    }
  }
  
  sortErr = data.frame( e.max=V$e.max, wp.inv=w_p[[2]](t.stat, V[[bound.col]], M  ) )
  
  sortErr$bkg = pmin(sortErr$e.max, sortErr$wp.inv )

  sortErr$importance = sortErr$e.max - sortErr$bkg
  sortErr = sortErr[order( sortErr$importance, decreasing=TRUE ),]

  q = 0
  totErr = sum( sortErr$bkg )
  if ( sum( sortErr$importance ) + totErr < threshold ) {
    q = nrow(sortErr)
  } else {
    while( totErr < threshold ) {
      q = q + 1;
      totErr = totErr + sortErr$importance[q]
    }
  }
  q
}


find.stark.SRS.p = function( N, n, q ) {
                                        ## Find the p-value for a given q, n, and N
  
                                        ## N = total number of precints
                                        ## n = total number of audited precints (must be less than N)
                                        ## q = min number of precints that could hold taint to flip election
  
  if ( n >= N ) stop( "Audit size is equal or greater than population size." )
  
                                        #cat( "N=",N, "n=",n,"q=",q,"\n")
  if ( N-q < n ) {
    0
  } else {
                                        # q bad precints, N-q good ones, n draws -> chance of getting
                                        # no bad precints in audit
    phyper( 0, q, N-q, n )
  }
}



stark.test.Z = function( Z, 
  calc.e_p=calc.pairwise.e_p,
  w_p = weight.function("no.weight"),
  max_err = maximumMarginBound,
  bound.col = Z$tot.votes.col,
  strat.col = NULL,
  drop=NULL,
  strat.method = NULL,
  err.override = NULL,
  n = NULL, t = NULL, q = NULL
  ) {
  
  ## Param Z: The object holding all the voting information
  ##    In particular: Z$V:     The table of reported votes
  ##                   Z$audit: The table of audits as differences from 
  ##                            recorded votes
  ## Param calc.e_p: The Function used to calculate maximum error bounds
  ## Param      w_p: The function used to calculate weights of error
  ##                 (A list of two functions)
  ##   	    max_err: Function to compute _max_ error bounds for each precint
  ##      strat.col: Name of column that determines how to stratify
  ##          t,n,q: Elements of the test statistic.  Can pass to avoid computation
  ##                 if those values are already known (e.g., for a simulation)
  
   ##stopifnot( !is.null( strat.col) && is.null(strat.method) )
  stopifnot( is.null(drop) || !is.null(Z$V[[drop]]) )
  
  if ( is.null( t ) ) { 	
    res = compute.stark.t( Z, bound.col, calc.e_p, w_p, err.override=err.override,
      return.revised.audit = TRUE)
    Z$audit = res[[2]]
    t = res[[1]]
  }
  
                                        #cat( "weight function = " )
                                        #print( w_p )
  
  if ( is.null( q ) ) {
    Z$V$e.max = max_err( Z )
    q = find.q( Z$V, t, bound.col, Z$margin, w_p=w_p, drop=drop )
  }
  
  
                                        # computing sample size.
  if ( is.null( n ) ) {
    n = length(Z$audit[[1]])
    passed_n = NULL
  } else {
    ## compute scaling of n
    passed_n = n / length(Z$audit[[1]])
  }
  
  
  ## Once q has been computed, calculate final p-values.
  
  
  if ( is.null( strat.col ) ) {
	if ( !is.null(drop) ) {
 	 	eff.N = Z$N - sum(Z$V[[drop]])
 	} else {
  		eff.N = Z$N
  	}
    p.value = find.stark.SRS.p( eff.N, n, q )
    method = "Stark's Election Audit Test (SRS Audit)"
    DAT = paste( "# precincts = ", eff.N, " of ", Z$N, ", f/C=", Z$f, " of ", 
      Z$C, ", n=",n, sep="")
    
  } else {
	## Analysis as if the minimum sampling fraction applied
	##  everywhere.  Bound from sampling *with* replacement

	if ( !is.null(drop) ) {
		stop( "dropped column and stratification code not working right--no drop col allowed" )
	}
	
	## Find stratification levels ###
    strat = find.stratification( Z$V, Z$audit, strat.col )
        	
    ## The following takes out completely missed strata for which there
    ## is no audit data, and assigns maximum error to it
    not.audited = strat[[strat.col]][strat$sample<0.001]
    big.N = eff.N
    if ( length( not.audited ) > 0 ) {
      
      skipped = Z$V[[strat.col]] %in% not.audited
      cat( "# precincts missed", length(Z$V[skipped,]$e.max), " prec in ", length( not.audited ), "\n")
      M.red = sum( Z$V[skipped,]$e.max )
      
                                        #print( Z$V[skipped,] )
      
      qp = find.q( Z$V[!skipped, ], t, bound.col, Z$margin, threshold=(1-M.red), w_p=w_p )
                                        #               warning( "Some strata have no audit data--will subtract max error to get B=", 
                                        #               round( (1-M.red), digits=3), ", and q reduced to ", qp, " from ", q )
      
      big.N = big.N - (q - qp) 
      cat( " del q ", (q-qp), round( M.red, digits=3 ), "\n" )
      q = qp
    }
    
    ## find effective n by taking smallest audit percentage of strata (other than 0).
    n <- floor(big.N*min(strat$sample[strat$sample>=0.001]))
    if ( !is.null(passed_n) ) { 
      n = floor( passed_n *big.N*min(strat$sample[strat$sample>=0.01]))
    }
    p.value = pbinom(0,n,q/big.N)
    method = "Stark's Election Audit Test (Conservative Stratafication)"
    DAT = paste( "# precincts = ", Z$N, ", f/C=", Z$f, " of ", 
      Z$C, ", n=",length(Z$audit[[1]]), 
      " (", n, ")", sep="")
    if ( big.N != Z$N ) {
      DAT = paste( DAT, " N'=", big.N )
    }
  } # end stratification block
  
  names(t) = "max err";
  names(q) = "q";
  
  structure(list(statistic = t, p.value = p.value,
                 method=method, parameter=q,
                 data.name = DAT,
                 Z=Z, n=n), class = "htest")
}


stark.test = function( votes, audits, C.names=NULL, f=1, pool=TRUE, pairwise=FALSE, ... ) {
  ## Do the entire test. Basically a driver function that sets up 'Z' matrix and passes buck
  ## to the stark.test.Z
  ## 
  ## param       votes:  Table of reported votes.  Each row is precinct.
  ## param      audits:  Table of audits.  Each row is precinct.  Table reports overstatement by
  ##                     candidate.
  ## param     C.names:  Names of candidates (and names of cor columns in votes and audits tables.
  ## param           f:  The number of winners
  ## param    pairwise:  if TRUE then do a pairwise test for all pairs and return
  ##                     highest p-value
  ## param     C.names:  if NULL will derive from cols 2 on of votes
  ## param   strat.col:  if NULL will not stratify

  if ( pairwise ) {
    stark.pairwise.test( votes, audits, C.names=C.names, f=f, pool=pool, ... )
  } else {
    Z = make.Z(votes, C.names, f, audit=audits, pool=pool )

    T = stark.test.Z( Z, ... )
    
    T
  }
}


find.stratification = function( D, aud, strat.col ) {
  ## Finding how audit interacted with stratification levels for a table of votes and audits
  ## param     D:    Table of votes
  ##         aud:    Table of audit data
  ##     strat.col:  The column to use that identifies the stratification levels
  ## Return: Table of strata, for each strata (row) the name of the strata,
  ##         the number of precincts, the number of audited precincts
  ##         and percent of precincts audited is returned.
  if ( !(strat.col %in% names(D) ) ) {
    stop( "Cannot find stratification since strat.col '", strat.col, "' is not in vote data object." )
  }
  
  if ( !( strat.col %in% names(aud) ) ) {
  	stop( "Cannot find stratification since strat.col '", strat.col, "' is not in audit data object." )
  }
  
  nPrecincts <- table( D[[strat.col]] )
  audPrecincts = table( aud[[strat.col]] )
  tbl = merge( nPrecincts, audPrecincts, by="Var1", all.x=TRUE )

  names(tbl) = c(strat.col, "n", "audit" )
  tbl$audit[is.na(tbl$audit)] = 0
  tbl$sample = tbl$audit / tbl$n
  
  tbl
}

stark.pairwise.test = function( votes, audits, C.names=NULL, f=1, pool=TRUE, ... ) {
  ## Pairwise test: Look at all pairs of winners and losers, compute the p-value
  ## according to methods defined by passed parameters, and then return the
  ## worst p-value, using the bound that the chance that _all_ the statistics would be
  ## that small under their respective H0 is bounded by the chance that the hardest-to-
  ## detect statistic would be that small.
  ## param     C.names:  if NULL will derive from cols 2 on of votes
  ## param   strat.col:  if NULL will not stratify
  
  all.Z= make.Z(votes, C.names, f, audit=audits, pool=pool )
  
  max.pv = -1
  best.T = NULL

  for ( w in all.Z$winners ) {
    for ( l in all.Z$losers ) {
      ##cat( "Examining ", w, " beating ", l, "\n" )
      Z = make.Z( all.Z$V, C.names=c(w, l), 1, audit=all.Z$audit, pool=FALSE )
      T = stark.test.Z( Z, ... )
      ##print(T)
      if ( T$p.value > max.pv ) {
        max.pv = T$p.value
        best.T = T
      }
      ##print( "finished iteration" )
    }
  }

  best.T$method = paste( best.T$method, "(pairwise)" )
  best.T$data.name = paste( best.T$data.name, " (overall F/C=",all.Z$f, "/", all.Z$C, ")", sep="" )
  
  best.T
}


audit.totals.to.OS = function( Z, audit ) {
	stopifnot( !is.null( Z$V$PID ) )
	stopifnot( !is.null( audit$PID ) )
	stopifnot( all( Z$C.names %in% colnames(audit) ) )
	stopifnot( all( audit$PID %in% rownames(Z$V) ) )
	audit$PID = as.character(audit$PID)
	
	audit[ Z$C.names ] = Z$V[ audit$PID, Z$C.names ] - audit[ Z$C.names ]  
	if ( is.null( audit$tot.votes ) ) {
		audit$tot.votes = Z$V[ audit$PID, "tot.votes" ]
	}
	
	audit
}




