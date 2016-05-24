## load miscTools package
library( "miscTools" )

## create a matrix
m <- matrix( 1:9, 3 )

# insert rows
print( insertRow( m, 1, 10:12 ) )
print( insertRow( m, 2, 10:12 ) )
print( insertRow( m, 3, 10:12 ) )
print( insertRow( m, 4, 10:12 ) )

# insert columns
print( insertCol( m, 1, 10:12 ) )
print( insertCol( m, 2, 10:12 ) )
print( insertCol( m, 3, 10:12 ) )
print( insertCol( m, 4, 10:12 ) )

# insert rows with name
print( insertRow( m, 1, 10:12, "R0" ) )
print( insertRow( m, 2, 10:12, "R1a" ) )
print( insertRow( m, 3, 10:12, "R2a" ) )
print( insertRow( m, 4, 10:12, "R4" ) )

# insert columns with name
print( insertCol( m, 1, 10:12, "C0" ) )
print( insertCol( m, 2, 10:12, "C1a" ) )
print( insertCol( m, 3, 10:12, "C2a" ) )
print( insertCol( m, 4, 10:12, "C4" ) )

## add row names and column names
rownames( m ) <- c( "R1", "R2", "R3" )
colnames( m ) <- c( "C1", "C2", "C3" )

# insert rows
print( insertRow( m, 1, 10:12 ) )
print( insertRow( m, 2, 10:12 ) )
print( insertRow( m, 3, 10:12 ) )
print( insertRow( m, 4, 10:12 ) )

# insert columns
print( insertCol( m, 1, 10:12 ) )
print( insertCol( m, 2, 10:12 ) )
print( insertCol( m, 3, 10:12 ) )
print( insertCol( m, 4, 10:12 ) )

# insert rows with name
print( insertRow( m, 1, 10:12, "R0" ) )
print( insertRow( m, 2, 10:12, "R1a" ) )
print( insertRow( m, 3, 10:12, "R2a" ) )
print( insertRow( m, 4, 10:12, "R4" ) )

# insert columns with name
print( insertCol( m, 1, 10:12, "C0" ) )
print( insertCol( m, 2, 10:12, "C1a" ) )
print( insertCol( m, 3, 10:12, "C2a" ) )
print( insertCol( m, 4, 10:12, "C4" ) )

# insert a row to a single-column matrix (example provided by Max Gordon)
insertRow( matrix( 1:3, ncol=1 ), 2, 4 )

# insert a column to a single-row matrix
insertCol( matrix( 1:3, nrow=1 ), 2, 4 )
