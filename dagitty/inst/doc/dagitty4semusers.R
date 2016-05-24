## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(comment = "")
library(dagitty)

## ------------------------------------------------------------------------
g1 <- dagitty( "dag {
	W1 -> Z1 -> X -> Y
	Z1 <- V -> Z2
	W2 -> Z2 -> Y
	X <-> W1 <-> W2 <-> Y
}")

g2 <- dagitty( "dag {
	Y <- X <- Z1 <- V -> Z2 -> Y
	Z1 <- W1 <-> W2 -> Z2
	X <- W1 -> Y
	X <- W2 -> Y
}")

plot(graphLayout(g1))

## ------------------------------------------------------------------------
print( impliedConditionalIndependencies( g1 ) )

## ------------------------------------------------------------------------
print( adjustmentSets( g1, "Z1", "X", effect="direct" ) )

## ------------------------------------------------------------------------
print( adjustmentSets( g2, "X", "Y", effect="direct" ) )

## ------------------------------------------------------------------------
for( n in names(g1) ){
	for( m in children(g1,n) ){
		a <- adjustmentSets( g1, n, m, effect="direct" )
		if( length(a) > 0 ){
			cat("The coefficient on ",n,"->",m,
				" is identifiable controlling for:\n",sep="")
			print( a, prefix=" * " )
		}
	}
}

## ------------------------------------------------------------------------
print( adjustmentSets( g1, "X", "Y" ) )

## ------------------------------------------------------------------------
print( adjustmentSets( g2, "X", "Y" ) )

## ------------------------------------------------------------------------
for( n in names(g1) ){
	for( m in setdiff( descendants( g1, n ), n ) ){
		a <- adjustmentSets( g1, n, m )
		if( length(a) > 0 ){
			cat("The total effect of ",n," on ",m,
				" is identifiable controlling for:\n",sep="")
			print( a, prefix=" * " )
		}
	}
}

## ------------------------------------------------------------------------
for( n in names(g1) ){
	for( m in children(g1,n) ){
		iv <- instrumentalVariables( g1, n, m )
		if( length( iv ) > 0 ){
			cat( n, m, "\n" )
			print( iv , prefix=" * " )
		}
	}
}

