library( "miscTools" )


## matrix
m <- matrix( 1:24, nrow = 6, ncol = 4 )

cm1 <- colMedians( m )
print( cm1 )

rm1 <- rowMedians( m )
print( rm1 )

all.equal( cm1, rowMedians( t( m ) ) )
all.equal( rm1, colMedians( t( m ) ) )


## data.frame
data( "Electricity", package = "Ecdat" )
Electricity <- Electricity[ 1:20, ]

cm2 <- colMedians( Electricity )
print( cm2 )

rm2 <- rowMedians( Electricity )
print( rm2 )

all.equal( cm2, rowMedians( t( Electricity ) ) )
all.equal( rm2, colMedians( t( Electricity ) ) )

# array (3 dimensions)
a3 <- array( 1:24, dim = c(4,3,2),
   dimnames = list( c("a","b","c","d"), c("A","B","C"), c("x","y") ) )
colMedians( a3 )
all.equal( median( a3[ , "B", "y" ] ), colMedians( a3 )[ "B", "y" ] )

# array (4 dimensions)
a4 <- array( 1:120, dim = c(5,4,3,2),
   dimnames = list( c("a","b","c","d","e"), c("A","B","C","D"),
   c("x","y","z"), c("Y","Z") ) )
colMedians( a4 )
all.equal( median( a4[ , "B", "x", "Z" ] ), colMedians( a4 )[ "B", "x", "Z" ] )
