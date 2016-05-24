library( "frontier" )
library( "plm" )
options( digits = 5 )

## data set of rice producers in the Philippines
data( riceProdPhil )
riceProdPhil$lPROD  <- log( riceProdPhil$PROD )
riceProdPhil$lAREA  <- log( riceProdPhil$AREA )
riceProdPhil$lLABOR <- log( riceProdPhil$LABOR )
riceProdPhil$lNPK   <- log( riceProdPhil$NPK )
riceProdPhil$farm <-
   paste( "F_", ifelse( riceProdPhil$FMERCODE > 9, "", "0" ),
   riceProdPhil$FMERCODE, sep = "" )
riceProdPhil$year <- riceProdPhil$YEARDUM + 1998
riceProdPhil <- plm.data( riceProdPhil, c( "farm", "year" ) )


########## cross-section data #############

## without mu / zIntercept
# Error Components Frontier (ECF)
sbb5ecf <- sfa( lPROD ~ lAREA + lLABOR + lNPK,
   data = as.data.frame( riceProdPhil ) )
bb5ecf <- frontier( data = as.data.frame( riceProdPhil ),
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ) )
all.equal( sbb5ecf[-42], bb5ecf[-42], tol = 1e-3 )

# Efficiency Effects Frontier (EEF)
sbb5eef <- sfa( lPROD ~ lAREA + lLABOR + lNPK | - 1,
   data = as.data.frame( riceProdPhil ) )
bb5eef <- frontier( data = as.data.frame( riceProdPhil ),
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   zNames = NA )
all.equal( sbb5eef[-42], bb5eef[-42], tol = 1e-3 )
tmp <- efficiencies( sbb5eef, margEff = TRUE )

# Comparisons
round( rbind( coef( bb5ecf ), coef( bb5eef ) ), 2 )
all.equal( coef( bb5ecf ), coef( bb5eef ), tol = 1e-3 )
all.equal( vcov( bb5ecf ), vcov( bb5eef ), tol = 1e-3 )
all.equal( efficiencies( bb5ecf ), efficiencies( bb5eef ), tol = 1e-3 )


## with mu / zIntercept
# Error Components Frontier (ECF)
sbb6ecf <- sfa( lPROD ~ lAREA + lLABOR + lNPK,
   data = as.data.frame( riceProdPhil ), truncNorm = TRUE, muBound = 0 )
bb6ecf <- frontier( data = as.data.frame( riceProdPhil ),
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   truncNorm = TRUE, muBound = 0 )
all.equal( sbb6ecf[-42], bb6ecf[-42], tol = 1e-3 )

# Efficiency Effects Frontier (EEF)
sbb6eef <- sfa( lPROD ~ lAREA + lLABOR + lNPK | 1,
   data = as.data.frame( riceProdPhil ) )
bb6eef <- frontier( data = as.data.frame( riceProdPhil ),
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   zIntercept = TRUE, zNames = NA )
all.equal( sbb6eef[-42], bb6eef[-42], tol = 1e-3 )
tmp <- efficiencies( sbb6eef, margEff = TRUE )

# Comparisons
round( rbind( coef( bb6ecf ), coef( bb6eef )[ c( 1:4, 6:7, 5 ) ] ), 2 )
all.equal( efficiencies( bb6ecf ), efficiencies( bb6eef ), tol = 1e-3 )


############ panel data ###############

## without mu / zIntercept
# Error Components Frontier (ECF)
b5ecf <- bb5ecf

# Efficiency Effects Frontier (EEF)
sb5eef <- sfa( lPROD ~ lAREA + lLABOR + lNPK | - 1,
   data = riceProdPhil )
b5eef <- frontier( data = riceProdPhil,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   zNames = NA )
all.equal( sb5eef[-42], b5eef[-42], tol = 1e-3 )
all.equal( b5eef[ -c( 4, 5, 20, 33, 34, 42 ) ], bb5eef[ -c( 4, 5, 20, 33, 34, 42 ) ], 
   tol = 1e-3 )
all.equal( c( t( residuals( b5eef ) ) ), c( residuals( bb5eef ) ), tol = 1e-3 )

# Comparisons
round( rbind( coef( b5ecf ), coef( b5eef ) ), 2 )
all.equal( coef( b5ecf ), coef( b5eef ), tol = 1e-3 )
all.equal( vcov( b5ecf ), vcov( b5eef ), tol = 1e-3 )
all.equal( c( efficiencies( b5ecf ) ), c( t( efficiencies( b5eef ) ) ), 
   tol = 1e-3 )


## without mu / zIntercept
# Error Components Frontier (ECF)
sb6ecf <- sfa( lPROD ~ lAREA + lLABOR + lNPK,
   data = as.data.frame( riceProdPhil ), truncNorm = TRUE, muBound = Inf )
b6ecf <- frontier( data = as.data.frame( riceProdPhil ),
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   truncNorm = TRUE, muBound = Inf )
all.equal( sb6ecf[-42], b6ecf[-42], tol = 1e-3 )
all.equal( b6ecf[-42], bb6ecf[-42], tol = 1e-3 )

# Efficiency Effects Frontier (EEF)
sb6eef <- sfa( lPROD ~ lAREA + lLABOR + lNPK | 1,
   data = riceProdPhil )
b6eef <- frontier( data = riceProdPhil,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   zIntercept = TRUE, zNames = NA )
all.equal( sb6eef[-42], b6eef[-42], tol = 1e-3 )
all.equal( b6eef[ -c( 4, 5, 20, 33, 34, 42 ) ], bb6eef[ -c( 4, 5, 20, 33, 34, 42 ) ], 
   tol = 1e-3 )
all.equal( c( efficiencies( b6ecf ) ), c( efficiencies( bb6eef ) ), tol = 1e-3 )
all.equal( c( residuals( b6ecf ) ), c( residuals( bb6eef ) ), tol = 5e-3 )

# Comparisons
round( rbind( coef( b6ecf ), coef( b6eef )[ c( 1:4, 6:7, 5 ) ] ), 2 )
all.equal( c( efficiencies( b6ecf ) ), c( t( efficiencies( b6eef ) ) ), 
   tol = 1e-3 )
all.equal( c( residuals( b6ecf ) ), c( t( residuals( b6eef ) ) ), tol = 5e-3 )

